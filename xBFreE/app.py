# ##############################################################################
#                            GPLv3 LICENSE INFO                                #
#  Copyright (c) 2023.  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco  #
#                                                                              #
#  This program is free software; you can redistribute it and/or modify it     #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

import sys
from os.path import split
import logging
from pathlib import Path
import signal

try:
    from xBFreE.exceptions import GMXMMPBSA_ERROR, InputError, CommandlineError,MMPBSA_Error
    from xBFreE.mmpbsa.infofile import InfoFile
    from xBFreE.mmpbsa import main
    from xBFreE.mmpbsa.tester import run_test
    # from xBFreE.mmpbsa.commandlineparser import anaparser, testparser
    from xBFreE.mmpbsa.utils.misc import create_input_args
except ImportError:
    import os
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))


# Local methods
from mpi4py import MPI
if MPI.COMM_WORLD.size == 1:
    from xBFreE.mmpbsa.fake_mpi import MPI

_rank = MPI.COMM_WORLD.Get_rank()
_mpi_size = MPI.COMM_WORLD.Get_size()


def excepthook(exception_type, exception_value, tb):
    """
    Replaces sys.excepthook so fatal exceptions kill all MPI threads and we can
    control the printing of tracebacks. Those are helpful for debugging purposes,
    but may be unsightly to users. debug_printlevel set above controls this
    behavior
    """
    import traceback
    # global _stderr, _mpi_size, _rank
    if not isinstance(exception_type, MMPBSA_Error):
        traceback.print_tb(tb)
    _stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
    if _mpi_size > 1:
        _stderr.write('Error occurred on rank %d.' % _rank + os.linesep)
    _stderr.write('Exiting. All files have been retained.' + os.linesep)
    MPI.COMM_WORLD.Abort(1)


def interrupt_handler(signal, frame):
    """ Handles interrupt signals for a clean exit """
    global _stderr
    _stderr.write('\n%s interrupted! Program terminated. All files are kept.\n' %
                  os.path.split(sys.argv[0])[1])
    MPI.COMM_WORLD.Abort(1)


def setup_run():
    """
    Replace the uncaught exception handler to control traceback printing. Also
    add a signal handler for a SIGINT (Ctrl-C). However, we only want to do this
    if we're running gmx_MMPBSA -- for the API, we don't want to clobber the
    users' python environments like this.
    """
    sys.excepthook = excepthook
    signal.signal(signal.SIGINT, interrupt_handler)


def gmxmmpbsa():
    mmpbsa('gmx')
def ambermmpbsa():
    mmpbsa('amber')
def namdmmpbsa():
    mmpbsa('namd')

def mmpbsa(prog):
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(levelname)-7s] %(message)s")
    stream_handler.setFormatter(formatter)
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(levelname)-7s] %(message)s",
        handlers=[
            logging.FileHandler("gmx_MMPBSA.log", 'w'),
            stream_handler])
    # Just for compatibility as mpi4py works as serial when run without mpirun
    # (since v1.4.2)
    # from mpi4py import MPI
    # if MPI.COMM_WORLD.size == 1:
    #     from xBFreE.mmpbsa.fake_mpi import MPI
    # Set up error/signal handlers
    # main.setup_run()

    # Instantiate the main MMPBSA_App
    app = main.MMPBSA_App(MPI, prog)

    # Read the command-line arguments
    try:
        app.get_cl_args(sys.argv[1:])
    except CommandlineError as e:
        sys.stderr.write('%s: %s' % (type(e).__name__, e) + '\n')
        sys.exit(1)

    if app.FILES.createinput is not None:
        args_list = create_input_args(app.FILES.createinput)
        app.input_file.print_contents('mmpbsa.in', args_list)
        logging.info(f'Input file creation successful. Path: {Path("mmpbsa.in").absolute()}')
        sys.exit(0)

    # Perform our MMPBSA --clean now
    if app.FILES.clean:
        logging.info('Cleaning temporary files and quitting.\n')
        app.remove(-1)
        sys.exit(0)

    # See if we wanted to print out our input file options
    if app.FILES.infilehelp:
        app.input_file.print_contents(sys.stdout)
        sys.exit(0)

    # If we're not rewriting output do whole shebang, otherwise load info and parms
    # Throw up a barrier before and after running the actual calcs
    if not app.FILES.rewrite_output:
        try:
            app.read_input_file()
        except InputError as e:
            sys.stderr.write('%s: %s' % (type(e).__name__, e) + '\n')
            sys.stderr.write('  Enter `%s --help` for help\n' %
                             (split(sys.argv[0])[1]))
            sys.exit(1)
        app.process_input()
        app.check_for_bad_input()
        app.make_prmtops()
        app.loadcheck_prmtops()
        app.file_setup()
        app.run_mmpbsa()
    # If we are rewriting output, load the info and check prmtops
    else:
        info = InfoFile(app, True)
        info.read_info()
        app.loadcheck_prmtops()

    # Now we parse the output, print, and finish
    app.parse_output_files()
    app.write_final_outputs()
    app.finalize()


def gmxmmpbsa_ana():
    try:
        from PyQt6.QtWidgets import QApplication
        pyqt = True
    except:
        try:
            from PyQt5.QtWidgets import QApplication
            pyqt = True
        except:
            pyqt = False
    finally:
        if not pyqt:
            GMXMMPBSA_ERROR('Could not import PyQt5/PyQt6. gmx_MMPBSA_ana will be disabled until PyQt5/PyQt6 is '
                            'installed')

    from xBFreE.mmpbsa.analyzer.gui import GMX_MMPBSA_ANA
    from xBFreE.mmpbsa.analyzer.utils import get_files

    app = QApplication(sys.argv)
    app.setApplicationName('gmx_MMPBSA Analyzer (gmx_MMPBSA_ana)')
    try:
        parser = anaparser.parse_args(sys.argv[1:])
    except CommandlineError as e:
        GMXMMPBSA_ERROR('%s: %s' % (type(e).__name__, e))
        sys.exit(1)
    ifiles = get_files(parser)
    w = GMX_MMPBSA_ANA(ifiles)
    w.show()
    sys.exit(app.exec())


def gmxmmpbsa_test():
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(levelname)-7s] %(message)s",
        handlers=[
            logging.FileHandler("gmx_MMPBSA_test.log", 'w'),
            stream_handler])
    try:
        parser = testparser.parse_args(sys.argv[1:])
    except CommandlineError as e:
        GMXMMPBSA_ERROR('%s: %s' % (type(e).__name__, e))
        sys.exit(1)
    run_test(parser)

if __name__ == '__main__':


    logging.info('Finished')

    # gmxmmpbsa()
    gmxmmpbsa_ana()