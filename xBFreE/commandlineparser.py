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

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawTextHelpFormatter
from xBFreE.exceptions import GMXMMPBSA_ERROR
from xBFreE import __version__
import sys
import os
from pathlib import Path

from xBFreE.mmpbsa.commandlineparser import gmx_parser, amber_parser

if os.getenv('AMBERHOME'):
    # rism_xvv = os.path.join(os.getenv('AMBERHOME'), 'dat', 'mmpbsa', 'spc.xvv')
    path_to_file = os.path.join(os.getenv('AMBERHOME'), 'AmberTools', 'test', 'rism1d', 'tip3p-kh', 'tip3p.xvv.save')
    if os.path.exists(path_to_file):
        rism_xvv = path_to_file
    else:
        rism_xvv = os.path.join(Path(__file__).parent.joinpath('data'), 'xvv_files', 'tip3p.xvv')
else:
    rism_xvv = os.path.join(Path(__file__).parent.joinpath('data'), 'xvv_files', 'tip3p.xvv')


class GMXMMPBSA_ArgParser(ArgumentParser):
    """
    ArgumentParser subclass to redirect exit output to the xBFreE logging
    """

    def exit(self, status=0, message=None):
        if message:
            GMXMMPBSA_ERROR(message)
        sys.exit(status)


cmdparser = GMXMMPBSA_ArgParser(
    epilog='''xBFreEnergy is a tool to compute Binding Free Energy with different methods including MMPBSA and LIE''',
    description='''xBFreEnergy is a tool to compute Binding Free Energy with different methods including MMPBSA and 
    LIE''',
    formatter_class=ArgumentDefaultsHelpFormatter
)
cmdparser.add_argument('-v', '--version', action='version',
                    version='''%%(prog)s %s''' % __version__)


subparsers = cmdparser.add_subparsers(help='Methods to compute Binding Free Energy', dest='subparser')
gmxmmpbsa_parser = subparsers.add_parser(
    name='gmx_MMPBSA', parents=[gmx_parser],
    epilog='''xBFreEnergy is a tool to compute Binding Free Enrgy with different methods including MMPBSA and LIE''',
    description='''This is the core of gmx_MMPBSA and it will do all the calculations''',
    help='MMPBSA calculations for GROMACS'
)
ambermmpbsa_parser = subparsers.add_parser(
    name='amber_MMPBSA', parents=[amber_parser],
    epilog='''xBFreEnergy is a tool to compute Binding Free Enrgy with different methods including MMPBSA and LIE''',
    description='''This is the core of gmx_MMPBSA and it will do all the calculations''',
    help='MMPBSA calculations for AMBER'
)
# namdmmpbsa_parser = subparsers.add_parser(
#     name='namd_MMPBSA', parents=[gmx_parser],
#     epilog='''xBFreEnergy is a tool to compute Binding Free Enrgy with different methods including MMPBSA and LIE''',
#     description='''This is the core of gmx_MMPBSA and it will do all the calculations''',
#     help='MMPBSA calculations for NAMD',
#
# )
# charmmmmpbsa_parser = subparsers.add_parser(
#     name='charmm_MMPBSA', parents=[gmx_parser],
#     epilog='''xBFreEnergy is a tool to compute Binding Free Enrgy with different methods including MMPBSA and LIE''',
#     description='''This is the core of gmx_MMPBSA and it will do all the calculations''',
#     help='MMPBSA calculations for CHARMM'
# )


# parser.parse_args(['gmx_MMPBSA', '-h'])


