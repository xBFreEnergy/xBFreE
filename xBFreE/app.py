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

from xBFreE.main import xBFreE_App
from xBFreE.mmpbsa.app import mmpbsa
from xBFreE.utils.misc import create_input_args
import os


# Local methods
from mpi4py import MPI
if MPI.COMM_WORLD.size == 1:
    from xBFreE.fake_mpi import MPI

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

setup_run()

def run_xbfree():
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(levelname)-7s] %(message)s")
    stream_handler.setFormatter(formatter)
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(levelname)-7s] %(message)s",
        handlers=[
            logging.FileHandler("xBFreE.log", 'w'),
            stream_handler])

    xbfree_app = xBFreE_App(MPI)
    files = xbfree_app.get_cmd_args(sys.argv[1:])

    method = "mmpbsa" if "mmpbsa" in xbfree_app.FILES.subparser.lower() else None

    # Remove all generated files in the directory
    if xbfree_app.FILES.clean:
        logging.info('Cleaning temporary files and quitting.\n')
        xbfree_app.remove(-1)
        sys.exit(0)

    if xbfree_app.FILES.createinput is not None:
        if method == "mmpbsa":
            args_list = create_input_args(xbfree_app.FILES.createinput, 'mmpbsa')
            from xBFreE.input.mmpbsa import input_file
            input_file.print_contents('mmpbsa.in', args_list)
            logging.info(f'Input file creation successful. Path: {Path("mmpbsa.in").absolute()}')
        # TODO: LIE method -- not implemented yet

        sys.exit(0)

    if method == 'mmpbsa':
        mmpbsa(files, MPI)


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