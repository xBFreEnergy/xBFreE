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

import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawTextHelpFormatter
from xBFreE import __version__, __mmpbsa_version__, __ambertools_version__
from types import SimpleNamespace

from pathlib import Path

from xBFreE.exceptions import GMXMMPBSA_ERROR


def check_arg(str_suffix, path=False):
    def decorator(f):
        def gmx_MMPBSA_file(*args):
            result = f(*args)
            if path and not Path(result).exists():
                GMXMMPBSA_ERROR(f'"{result}" do not exist or is inaccessible')
            if Path(result).suffix not in str_suffix:
                GMXMMPBSA_ERROR(f'"{result}" does not correspond to the required structure format {str_suffix}')
            return result

        return gmx_MMPBSA_file

    return decorator


@check_arg(['.tpr', '.pdb'], True)
def gmx_structure(arg):
    return arg


@check_arg(['.pdb'], True)
def pdb(arg):
    return arg


@check_arg(['.xtc', '.trr', '.pdb'], True)
def gmx_trajectory(arg):
    return arg


@check_arg(['.mdcrd', '.crd', '.nc'], True)
def amber_trajectory(arg):
    return arg


@check_arg(['.dcd', '.nc'], True)
def namd_trajectory(arg):
    return arg


@check_arg(['.top'], True)
def gmx_topology(arg):
    return arg


@check_arg(['.prmtop', '.top', '.parm7'], True)
def amber_topology(arg):
    return arg


@check_arg(['.psf'], True)
def namd_topology(arg):
    return arg


@check_arg(['.mol2'], True)
def mol2(arg):
    return arg


@check_arg(['.ndx'], True)
def index(arg):
    return arg


def mask(arg):
    if not arg.startswith(':'):
        arg = f":{arg}"
    return arg.replace(' ', '')


class OptionList(SimpleNamespace):
    """
    Just a container to hold the command-line options. Necessary when reading in
    a gmx_MMPBSA info file to have a container to load the results from the
    parser.
    """
    pass


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

    def exit(self, status=0, message=None):
        if message:
            GMXMMPBSA_ERROR(message)
        sys.exit(status)


description = (
    "gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy calculations"
    " with GROMACS files. This program is an adaptation of Amber's MMPBSA.py and essentially works as such. "
    "gmx_MMPBSA works with any GROMACS version. This program will calculate binding free energies using "
    "end-state free energy methods on an ensemble of snapshots using a variety of implicit solvent models.")

complex_group_des = ("Complex files and info that are needed to perform the calculation. If the receptor and/or the "
                     "ligand info is not defined, we generate them from that of the complex.")

receptor_group_des = ("Receptor files and info that are needed to perform the calculation. If the receptor info is "
                      "not defined, we generate it from that of the complex.")

ligand_group_des = '''Ligand files and info that are needed to perform the calculation. If the ligand are not defined, 
                       we generate it from that of the complex.'''


def mmpbsa_parser(epilog='', description=''):
    parser = GMXMMPBSA_ArgParser(
        epilog='''xBFreEnergy is a tool to compute Binding Free Enrgy with different methods including MMPBSA and LIE''',
        description='''This is the core of gmx_MMPBSA and it will do all the calculations''',
        add_help=False
        # formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-O', '--overwrite', default=False, action='store_true',
                        help='Allow output files to be overwritten', dest='overwrite')

    parser.add_argument('--create_input', dest='createinput', choices=['gb', 'pb', 'rism', 'ala', 'decomp', 'nmode',
                                                                       'gbnsr6', 'all'],
                        nargs='*', help='Create an new input file with selected calculation type.')
    parser.add_argument('--rewrite-output', dest='rewrite_output', default=False,
                        action='store_true', help='''Do not re-run any calculations,
                          just parse the output files from the previous calculation and
                          rewrite the output files.''')
    parser.add_argument('-s', '--stability', dest='stability', action='store_true', default=False,
                        help='''Perform stability calculation. Only the complex parameters are required. 
                        In any other case receptor and ligand parameters will be ignored''')
    parser.add_argument('-nogui', dest='gui', action='store_false', default=True,
                        help='No open gmx_MMPBSA_ana after all calculations finished')

    group = parser.add_argument_group('Input and Output Files', '''These options specify the input files and optional 
        output files.''')
    group.add_argument('-i', dest='input_file', metavar='FILE', help='MM/PBSA input file.')
    group.add_argument('-xvvfile', dest='xvvfile', help='XVV file for 3D-RISM.',
                       # default=rism_xvv
                       )
    group.add_argument('-o', dest='output_file', default='FINAL_RESULTS.dat', metavar='FILE',
                       help='Output file with selected method (MMPBSA or LIE) statistics.')
    group.add_argument('-do', dest='decompout', metavar='FILE', default='FINAL_DECOMP.dat',
                       help='Output file for decomposition statistics summary.')
    group.add_argument('-eo', dest='energyout', metavar='FILE',
                       help='''CSV-format output of all energy terms for every frame
                          in every calculation. File name forced to end in [.csv].
                          This file is only written when specified on the
                          command-line.''')
    group.add_argument('-deo', dest='dec_energies', metavar='FILE',
                       help='''CSV-format output of all energy terms for each printed
                          residue in decomposition calculations. File name forced to end
                          in [.csv]. This file is only written when specified on the
                          command-line.''')
    return parser


def add_miscellaneous_actions(parser):
    group = parser.add_argument_group('Miscellaneous Options')
    group.add_argument('-prefix', dest='prefix', default='_GMXMMPBSA_',
                       metavar='<file prefix>',
                       help='Prefix for intermediate files.')
    group.add_argument('--input-file-help', dest='infilehelp', action='store_true',
                       help='Print all available options in the input file.',
                       default=False)
    group.add_argument('--clean', dest='clean', action='store_true', default=False,
                       help='''Clean temporary files and quit.''')


# GUI parser
description = 'This program is part of gmx_MMPBSA and will show a workspace to analyze the gmx_MMPBSA results'
anaparser = ArgumentParser(epilog=f'gmx_MMPBSA is an effort to implement the GB/PB and others calculations in '
                                  f'GROMACS. \nBased on MMPBSA.py (version {__mmpbsa_version__}) and '
                                  f'AmberTools{__ambertools_version__}',
                           description=description,
                           formatter_class=ArgumentDefaultsHelpFormatter)
anaparser.add_argument('-v', '--version', action='version',
                       version='%%(prog)s %s based on MMPBSA version %s' % (__version__, __mmpbsa_version__))
group = anaparser.add_argument_group('Info file')
group.add_argument('-f', '--files', nargs='*', help='gmx_MMPBSA info files or container folder or list of them',
                   type=Path, default=None)
group.add_argument('-r', '--recursive', help='Search recursively in this folder at depth = 1', action='store_true',
                   default=False)

# tester parser
description = ('This program is part of gmx_MMPBSA and will allow you to run various gmx_MMPBSA examples easily.')
testparser = ArgumentParser(epilog=f'gmx_MMPBSA is an effort to implement the GB/PB and others calculations in '
                                   f'GROMACS. \nBased on MMPBSA.py (version {__mmpbsa_version__}) and '
                                   f'AmberTools{__ambertools_version__}',
                            description=description,
                            formatter_class=RawTextHelpFormatter)
testparser.add_argument('-v', '--version', action='version',
                        version='%%(prog)s %s based on MMPBSA version %s' % (__version__, __mmpbsa_version__))
group = testparser.add_argument_group('Test options')
group.add_argument('-t', dest='test', choices=list(range(19)) + [101], type=int, nargs='*', default=[2],
                   help='''\
The level the test is going to be run at. Multiple systems and analysis can be run at the same time.
      Nr. of Sys  
* 0      16     All -- Run all examples (Can take a long time!!!)
* 1      13     Minimal -- Does a minimal test with a set of systems and analyzes 
                that show that gmx_MMPBSA runs correctly. Only exclude 3drism, nmode
                protein-ligand MT because take a long time or are redundant
* 2       9     Fast -- Only the calculations that take a short time are run (Default)
[Systems]:
     Slow Frames
* 3    . | 10   Protein-Ligand (Single trajectory approximation)
* 4    . | 10   Protein-Protein
* 5    . | 10   Protein-DNA
* 6    x |  4   Protein-Membrane
* 7    . | 10   Protein-Glycan
* 8    x |  4   Metalloprotein-Peptide
* 9    . | 10   Protein-DNA-RNA-IONs-Ligand
* 10   x |  4   Protein-Ligand (CHARMM force field)
* 11   x |  4   Protein-ligand complex in membrane with CHARMMff 
[Analysis]:
     Slow Frames
* 12   . | 10   Alanine Scanning
* 13   . | 10   Stability calculation
* 14   . | 10   Decomposition Analysis
* 15   . | 16   Interaction Entropy approximation
* 16   . | 10   Protein-Ligand (Multiple trajectory approximation)
* 17   x |  4   Entropy calculation using Normal Mode approximation 
* 18   x |  4   Calculations using 3D-RISM approximation
''')
group.add_argument('-f', '--folder', help='Defines the folder to store all data', type=Path, default='.')
group.add_argument('-r', '--reuse', help='Defines the existing test forlder will be reuse', action='store_true')
group.add_argument('-ng', '--nogui', help='No open gmx_MMPBSA_ana after all calculations finished',
                   action='store_true', default=False)
group.add_argument('-n', '--num_processors', type=int, default=4,
                   help='Defines the number of processor cores you want to use with MPI per calculation. If the number '
                        'of frames is less than the number of cpus defined, the calculation will be performed with'
                        'the number of processors = number of frames')
