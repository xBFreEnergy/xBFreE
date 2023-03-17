from xBFreE.mmpbsa.commandlineparser.parsers import (mmpbsa_parser,
                                                 namd_trajectory, namd_topology, pdb,
                                                 complex_group_des, receptor_group_des, ligand_group_des,
                                                 add_miscellaneous_actions)

namd_parser = mmpbsa_parser()

group = namd_parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cp', dest='complex_prmtop', metavar='<Topology File>', type=namd_topology, default=None,
                   help='''The complex Topology file''')
group.add_argument('-cg', dest='complex_groups', metavar='aa selection', nargs=2, default=None, type=str,
                   help='Selection for receptor and ligand in complex topology file. Both Amber mask '
                        'and residue selection are allowed. The notation is as follows: "-cg '
                        '<Receptor group> <Ligand group>", ie. -cg :1-100 :101-120 (Amber mask) or A/1-100 B/101-120 '
                        '(Residue selection)')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', type=namd_trajectory,
                   help='''Complex trajectories. Make sure the trajectory is fitted and
                         pbc have been removed. Allowed formats: *.dcd (recommended), *.pdb,
                         (specify as many as you'd like).''')
group.add_argument('-cr', dest='reference_structure', metavar='<PDB File>', default=None, type=pdb,
                   help='''Complex Reference Structure file. This option is optional but recommended
                         (Use the PDB file used to generate the topology in NAMD). If not defined,
                         the chains ID assignment will be done automatically according to the
                         structure (can generate wrong mapping).''')
group = namd_parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rp', dest='receptor_prmtop', metavar='<Topology File>', type=namd_topology,
                   help='''Topology file of the unbound receptor.''')
group.add_argument('-rg', dest='receptor_group', metavar='index', default=None, type=str,
                   help='''Receptor group in receptor topology file. Both Amber mask '
                        'and residue selection are allowed. The notation is as follows: "-rg '
                        '<Receptor group>", ie. -rg :1-100 (Amber mask) or A/1-100 (Residue selection)''')
group.add_argument('-rt', dest='receptor_trajs', nargs='*', metavar='TRJ', type=namd_trajectory,
                   help='''Input trajectories of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.dcd (recommended), *.pdb, *.nc (specify as many
                         as you'd like).''')
group = namd_parser.add_argument_group('Ligand', ligand_group_des)
group.add_argument('-lp', dest='ligand_prmtop', metavar='<Topology File>', type=namd_topology,
                   help='''Topology file of the unbound ligand.''')
group.add_argument('-lg', dest='ligand_group', metavar='index', default=None, type=str,
                   help='''Ligand group in ligand topology file. Both Amber mask '
                        'and residue selection are allowed. The notation is as follows: "-lg '
                        '<Ligand group>", ie. -lg :101-120 (Amber mask) or B/101-120 (Residue selection)''')
group.add_argument('-lt', dest='ligand_trajs', nargs='*', metavar='TRJ', type=namd_trajectory,
                   help='''Input trajectories of the unbound ligand for multiple trajectory approach.
                         Allowed formats: *.dcd (recommended), *.pdb, *.nc (specify as many
                         as you'd like).''')

add_miscellaneous_actions(namd_parser)

# namd_parser.parse_args(['-h'])
