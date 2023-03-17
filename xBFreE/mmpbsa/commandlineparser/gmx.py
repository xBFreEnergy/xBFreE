from xBFreE.mmpbsa.commandlineparser.parsers import (mmpbsa_parser,
                                                 gmx_topology, gmx_structure, index, gmx_trajectory, pdb, mol2,
                                                 complex_group_des, receptor_group_des, ligand_group_des,
                                                 add_miscellaneous_actions)
from xBFreE.exceptions import GMXMMPBSA_ERROR
from argparse import ArgumentParser
import sys


class GMXMMPBSA_ArgParser(ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def exit(self, status=0, message=None):
        if message:
            GMXMMPBSA_ERROR(message)
        sys.exit(status)

gmx_parser = GMXMMPBSA_ArgParser(
    epilog='''xBFreEnergy is a tool to compute Binding Free Enrgy with different methods including MMPBSA and LIE''',
    description='''This is the core of gmx_MMPBSA and it will do all the calculations''',
    add_help=False
    # formatter_class=ArgumentDefaultsHelpFormatter
)

gmx_parser.add_argument('-O', '--overwrite', default=False, action='store_true',
                   help='Allow output files to be overwritten', dest='overwrite')

gmx_parser.add_argument('--create_input', dest='createinput', choices=['gb', 'pb', 'rism', 'ala', 'decomp', 'nmode',
                                                                   'gbnsr6', 'all'],
                    nargs='*', help='Create an new input file with selected calculation type.')
gmx_parser.add_argument('--rewrite-output', dest='rewrite_output', default=False,
                   action='store_true', help='''Do not re-run any calculations,
                  just parse the output files from the previous calculation and
                  rewrite the output files.''')
gmx_parser.add_argument('-s', '--stability', dest='stability', action='store_true', default=False,
                   help='''Perform stability calculation. Only the complex parameters are required. If
                         ligand is non-Protein (small molecule) type, then ligand *.mol2 file is 
                         required. In any other case receptor and ligand parameters will be ignored.
                         See description bellow''')
gmx_parser.add_argument('-nogui', dest='gui', action='store_false', default=True,
                   help='No open gmx_MMPBSA_ana after all calculations finished')

group = gmx_parser.add_argument_group('Input and Output Files', '''These options specify the input files and optional 
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


group = gmx_parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cs', dest='complex_tpr', metavar='<Structure File>', default=None, type=gmx_structure,
                   help='''Structure file of the complex. Allowed formats: *.tpr (recommended), *.pdb''')
group.add_argument('-ci', dest='complex_index', metavar='<Index File>', default=None, type=index,
                   help='Index file of the bound complex.')
group.add_argument('-cg', dest='complex_groups', metavar='index', nargs=2, default=None, type=int,
                   help='Groups of receptor and ligand in complex index file. The notation is as follows: "-cg '
                        '<Receptor group> <Ligand group>", ie. -cg 1 13')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', type=gmx_trajectory,
                   help='''Complex trajectories. Make sure the trajectory is fitted and
                         pbc have been removed. Allowed formats: *.xtc (recommended), *.trr, *.pdb
                         (specify as many as you'd like).''')
group.add_argument('-cp', dest='complex_top', metavar='<Topology>', default=None, type=gmx_topology,
                   help='''The complex Topology file. When it is defined -lm option is not needed''')
group.add_argument('-cr', dest='reference_structure', metavar='<PDB File>', default=None, type=pdb,
                   help='''Complex Reference Structure file. This option is optional but recommended 
                         (Use the PDB file used to generate the topology in GROMACS). If not defined,
                         the chains ID assignment (if the structure used in -cs does not have chain
                         IDs) will be done automatically according to the structure (can generate
                         wrong mapping).''')

group = gmx_parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rs', dest='receptor_tpr', metavar='<Structure File>', default=None, type=gmx_structure,
                   help='''Structure file of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.tpr (recommended), *.pdb''')
group.add_argument('-ri', dest='receptor_index', metavar='<Index File>', default=None, type=index,
                   help='Index file of the unbound receptor.')
group.add_argument('-rg', dest='receptor_group', metavar='index', default=None, type=int,
                   help='''Receptor group in receptor index file. Notation: "-lg <Receptor group>", 
                         e.g. -rg 1''')
group.add_argument('-rt', dest='receptor_trajs', nargs='*', metavar='TRJ', type=gmx_trajectory,
                   help='''Input trajectories of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.xtc (recommended), *.trr, *.pdb, *.gro (specify as many
                         as you'd like).''')
group.add_argument('-rp', dest='receptor_top', metavar='<Topology>', default=None, type=gmx_topology,
                   help='''Topology file of the receptor.''')

group = gmx_parser.add_argument_group('Ligand', ligand_group_des)
group.add_argument('-ls', dest='ligand_tpr', metavar='<Structure File>', default=None, type=gmx_structure,
                   help='''Structure file of the unbound ligand for multiple trajectory approach.
                         Allowed formats: *.tpr (recommended), *.pdb''')
group.add_argument('-li', dest='ligand_index', metavar='<Index File>', type=index,
                   default=None, help='Index file of the unbound ligand. Only if tpr file was define in -ls.')
group.add_argument('-lg', dest='ligand_group', metavar='index', default=None, type=int,
                   help='''Ligand group in ligand index file. Notation: "-lg <Ligand group>", 
                         e.g. -lg 13''')
group.add_argument('-lt', dest='ligand_trajs', nargs='*', metavar='TRJ', type=gmx_trajectory,
                   help='''Input trajectories of the unbound ligand for multiple trajectory approach. 
                         Allowed formats: *.xtc (recommended), *.trr, *.pdb, *.gro (specify as many
                         as you'd like).''')
group.add_argument('-lp', dest='ligand_top', metavar='<Topology>', default=None, type=gmx_topology,
                   help='''Topology file of the ligand.''')

group = gmx_parser.add_argument_group('Miscellaneous Options')
group.add_argument('-prefix', dest='prefix', default='_GMXMMPBSA_',
                   metavar='<file prefix>',
                   help='Prefix for intermediate files.')
group.add_argument('--input-file-help', dest='infilehelp', action='store_true',
                    help='Print all available options in the input file.',
                    default=False)
group.add_argument('--clean', dest='clean', action='store_true', default=False,
                   help='''Clean temporary files and quit.''')

# add actions to parser
# add_miscellaneous_actions(gmx_parser)

# gmx_parser.parse_args(['-h'])
