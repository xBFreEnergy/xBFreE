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

from argparse import ArgumentDefaultsHelpFormatter

from xBFreE import __version__
from xBFreE.mmpbsa.commandlineparser import gmx_parser, amber_parser, namd_parser, charmm_parser
from xBFreE.utils.misc import xBFreE_ArgParser

cmdparser = xBFreE_ArgParser(
    epilog='''xBFreEnergy is a tool to compute Binding Free Energy with different methods''',
    description='''xBFreEnergy is a tool to compute Binding Free Energy with different methods''',
    formatter_class=ArgumentDefaultsHelpFormatter
)
cmdparser.add_argument('-v', '--version', action='version',
                       version='''%%(prog)s %s''' % __version__)
cmdparser.add_argument('--clean', dest='xbfree_clean', help='Clean all files', action='store_true')

subparsers = cmdparser.add_subparsers(help='Methods to compute Binding Free Energy', dest='subparser')
gmxmmpbsa_parser = subparsers.add_parser(
    name='gmx_MMPBSA', parents=[gmx_parser],
    epilog='''xBFreE is a tool to compute Binding Free Enrgy with different methods''',
    description='''gmx_MMPBSA is the subcommand of xBFreE for performing MMPBSA calculations for GROMACS''',
    help='MMPBSA calculations for GROMACS',
    formatter_class=ArgumentDefaultsHelpFormatter
)
ambermmpbsa_parser = subparsers.add_parser(
    name='amber_MMPBSA', parents=[amber_parser],
    epilog='''xBFreE is a tool to compute Binding Free Enrgy with different methods''',
    description='''amber_MMPBSA is the subcommand of xBFreE for performing MMPBSA calculations for AMBER''',
    help='MMPBSA calculations for AMBER',
    formatter_class=ArgumentDefaultsHelpFormatter
)
namdmmpbsa_parser = subparsers.add_parser(
    name='namd_MMPBSA', parents=[namd_parser],
    epilog='''xBFreE is a tool to compute Binding Free Enrgy with different methods''',
    description='''namd_MMPBSA is the subcommand of xBFreE for performing MMPBSA calculations for NAMD''',
    help='MMPBSA calculations for NAMD',
    formatter_class=ArgumentDefaultsHelpFormatter
)
charmmmmpbsa_parser = subparsers.add_parser(
    name='charmm_MMPBSA', parents=[charmm_parser],
    epilog='''xBFreE is a tool to compute Binding Free Enrgy with different methods''',
    description='''charmm_MMPBSA is the subcommand of xBFreE for performing MMPBSA calculations for CHARMM''',
    help='MMPBSA calculations for CHARMM',
    formatter_class=ArgumentDefaultsHelpFormatter
)
