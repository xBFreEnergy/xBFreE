# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
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
import os
import logging
from pathlib import Path

from xBFreE.utils import misc
from xBFreE.fake_mpi import MPI as FakeMPI
from xBFreE.utils.misc import check4dup_args
from xBFreE.utils.timer import Timer
from xBFreE.commandlineparser import cmdparser
from xBFreE import __version__

_unbuf_stdout = misc.Unbuffered(sys.stdout)  # unbuffered stdout
_unbuf_stderr = misc.Unbuffered(sys.stderr)  # unbuffered stderr
_stdout = sys.stdout
_stderr = sys.stderr
_mpi_size = 1
_rank = 0
_MPI = FakeMPI()


class xBFreE_App:
    def __init__(self, MPI, stdout=None, stderr=None, size=None):

        self.FILES = {}
        _MPI = self.MPI = MPI
        self.rfolder = Path('xBFreE_RESULTS')
        self.INPUT = {}
        if stdout is None:
            _stdout = self.stdout = _unbuf_stdout
        else:
            _stdout = self.stdout = stdout

        if stderr is None:
            _stderr = self.stderr = _unbuf_stderr
        else:
            _stderr = self.stderr = stderr

        # MPI-related variables. Squash output for non-master threads
        _rank = self.mpi_rank = self.MPI.COMM_WORLD.Get_rank()
        self.master = self.mpi_rank == 0
        _mpi_size = self.mpi_size = self.MPI.COMM_WORLD.Get_size()
        if not self.master:
            self.stdout = open(os.devnull, 'w')
        if self.master:
            logging.info(f'Starting xBFreE {__version__}')
            misc.get_sys_info()

            # create the Result folder if not exists
            self.rfolder.mkdir(exist_ok=True)

        # Set up timers
        timers = [Timer() for _ in range(self.mpi_size)]
        self.timer = timers[self.mpi_rank]

        # Support possible threading for those that don't use MPI. However, if
        # mpi_size is > 1, just use the MPI mechanism instead
        if size is not None and self.mpi_size == 1:
            self.mpi_size = size

    def get_cmd_args(self, args=None):
        """
        Get args and initialize cmd variables
        """
        if args is None:
            args = sys.argv
        if self.master:
            text_args = ' '.join(args)
            mpi_cl = f'  mpirun -np {self.mpi_size} ' if self.mpi_size > 1 else '  '
            # print cmd in the log file
            logging.info('Command-line\n' + mpi_cl + 'xbfree ' + text_args + '\n')
            # show usage and exit when no args was defined
            if not args:
                cmdparser.print_usage()
                sys.exit(1)

            # check if there is any duplicated argument
            check4dup_args(args)
            self.FILES = cmdparser.parse_args(args)
        else:
            self.FILES = object()
        # Broadcast the FILES
        self.FILES = self.MPI.COMM_WORLD.bcast(self.FILES)
        # return the parser/ FILES to use it in the selected method
        return self.FILES

    def remove(self, flag):
        """ Removes all/temporary files """
        if not self.master:
            return
        misc.remove()