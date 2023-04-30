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

import logging
import os
import signal
import sys
from datetime import datetime
from pathlib import Path

from mpi4py import MPI

from xBFreE.exceptions import xBFreE_Error
from xBFreE.main import xBFreE_App
from xBFreE.mmpbsa.app import mmpbsa
from xBFreE.utils.misc import create_input_args

if MPI.COMM_WORLD.size == 1:
    from xBFreE.fake_mpi import MPI

_rank = MPI.COMM_WORLD.Get_rank()
_mpi_size = MPI.COMM_WORLD.Get_size()
_stdout = sys.stdout
_stderr = sys.stderr

def excepthook(exception_type, exception_value, tb):
    """
    Replaces sys.excepthook so fatal exceptions kill all MPI threads and we can
    control the printing of tracebacks. Those are helpful for debugging purposes,
    but may be unsightly to users. debug_printlevel set above controls this
    behavior
    """
    import traceback
    # global _stderr, _mpi_size, _rank
    if not isinstance(exception_type, xBFreE_Error):
        traceback.print_tb(tb)
    _stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
    if _mpi_size > 1:
        _stderr.write('Error occurred on rank %d.' % _rank + os.linesep)
    _stderr.write(f'Exiting. All files have been retained.{os.linesep}')
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
    if we're running xBFreE -- for the API, we don't want to clobber the
    users' python environments like this.
    """
    sys.excepthook = excepthook
    signal.signal(signal.SIGINT, interrupt_handler)


setup_run()


def run_xbfree():
    """
    Function tu launch xBFreE application
    - Create the logging handler
    - Get and process the args
    - Execute the sub-application according selected subcommand
    :return:
    """

    xbfree_log = Path("xBFreE.log")
    if _rank == 0 and xbfree_log.exists():
        xbfree_log.rename(f"xBFreE-{datetime.now().strftime('%Y-%m-%d-%H:%M:%S')}.log")

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
    # iterate over all posibles Path type variables in FILES to make them relative to wdir to chdir to
    # xBFreE_RESULTS/mmpbsa
    wdir = Path('.')
    for p in dir(files):
        if isinstance(getattr(files, p), list):
            templ = []
            for p1 in getattr(files, p):
                if isinstance(p1, Path):
                    if p1.is_relative_to(wdir):
                         templ.append(Path('../..').joinpath(p1).as_posix())
                    else:
                        templ.append(p1)
                else:
                    templ.append(p1)
            setattr(files, p, templ)
        else:
            arc = getattr(files, p)
            if isinstance(arc, Path) and arc.is_relative_to(wdir):
                setattr(files, p, Path('../..').joinpath(arc).as_posix())
    files.wdir = wdir.resolve()

    method = "mmpbsa" if "mmpbsa" in xbfree_app.FILES.subparser.lower() else None
    files.subwdir = Path('xBFreE_RESULTS', method).resolve()

    # Remove all generated files in the directory
    if xbfree_app.FILES.xbfree_clean:
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

    # See if we wanted to print out our input file options
    if xbfree_app.FILES.infilehelp:
        from xBFreE.input.mmpbsa import input_file
        input_file.print_contents(sys.stdout)
        sys.exit(0)

    if method == 'mmpbsa':
        mmpbsa(files, MPI)


# def gmxmmpbsa_test():
#     stream_handler = logging.StreamHandler()
#     stream_handler.setLevel(logging.INFO)
#     logging.basicConfig(
#         level=logging.DEBUG,
#         format="[%(levelname)-7s] %(message)s",
#         handlers=[
#             logging.FileHandler("gmx_MMPBSA_test.log", 'w'),
#             stream_handler])
#     try:
#         parser = testparser.parse_args(sys.argv[1:])
#     except CommandlineError as e:
#         GMXMMPBSA_ERROR('%s: %s' % (type(e).__name__, e))
#         sys.exit(1)
#     run_test(parser)


if __name__ == '__main__':
    logging.info('Finished')

    # gmxmmpbsa()
    # gmxmmpbsa_ana()
