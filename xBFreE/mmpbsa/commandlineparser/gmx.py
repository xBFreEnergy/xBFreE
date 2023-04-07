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

"""
This submodule contains the GROMACS parser
"""


from xBFreE.mmpbsa.commandlineparser.parsers import (mmpbsa_parser,
                                                     gmx_topology, gmx_structure, index, gmx_trajectory, pdb,
                                                     complex_group_des, receptor_group_des, ligand_group_des,
                                                     add_miscellaneous_actions)

gmx_parser = mmpbsa_parser()

group = gmx_parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cs', dest='complex_structure', metavar='<Structure File>', default=None, type=gmx_structure,
                   help="Structure file of the complex. Allowed formats: *.tpr (recommended), *.pdb")
group.add_argument('-ci', dest='complex_index', metavar='<Index File>', default=None, type=index,
                   help='Index file of the bound complex.')
group.add_argument('-cg', dest='complex_groups', metavar='index', nargs=2, default=None, type=int,
                   help='Groups of receptor and ligand in complex index file. The notation is as follows: "-cg '
                        '<Receptor group> <Ligand group>", ie. -cg 1 13')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', type=gmx_trajectory,
                   help="Complex trajectories. Make sure the trajectory is fitted and pbc have been removed. "
                        "Allowed formats: *.xtc, *.trr, *.pdb (specify as many as you'd like).")
group.add_argument('-cp', dest='complex_top', metavar='<Topology>', default=None, type=gmx_topology,
                   help="The complex Topology file. When it is defined -lm option is not needed")
group.add_argument('-cr', dest='reference_structure', metavar='<PDB File>', default=None, type=pdb,
                   help="Complex Reference Structure file. This option is optional but recommended (Use the PDB file "
                        "used to generate the topology in GROMACS). If not defined, the chains ID assignment (if the "
                        "structure used in -cs does not have chain IDs) will be done automatically according to the "
                        "structure (can generate wrong mapping).")

group = gmx_parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rs', dest='receptor_tpr', metavar='<Structure File>', default=None, type=gmx_structure,
                   help="Structure file of the unbound receptor for multiple trajectory approach. Allowed formats: "
                        "*.tpr (recommended), *.pdb")
group.add_argument('-ri', dest='receptor_index', metavar='<Index File>', default=None, type=index,
                   help='Index file of the unbound receptor.')
group.add_argument('-rg', dest='receptor_group', metavar='index', default=None, type=int,
                   help='Receptor group in receptor index file. Notation: "-lg <Receptor group>", e.g. -rg 1"')
group.add_argument('-rt', dest='receptor_trajs', nargs='*', metavar='TRJ', type=gmx_trajectory,
                   help="Input trajectories of the unbound receptor for multiple trajectory approach. "
                        "Allowed formats: *.xtc, *.trr, *.pdb, *.gro (specify as many as you'd like).")
group.add_argument('-rp', dest='receptor_top', metavar='<Topology>', default=None, type=gmx_topology,
                   help="Topology file of the receptor.")

group = gmx_parser.add_argument_group('Ligand', ligand_group_des)
group.add_argument('-ls', dest='ligand_tpr', metavar='<Structure File>', default=None, type=gmx_structure,
                   help="Structure file of the unbound ligand for multiple trajectory approach. "
                        "Allowed formats: *.tpr (recommended), *.pdb")
group.add_argument('-li', dest='ligand_index', metavar='<Index File>', type=index,
                   default=None, help='Index file of the unbound ligand. Only if tpr file was define in -ls.')
group.add_argument('-lg', dest='ligand_group', metavar='index', default=None, type=int,
                   help='Ligand group in ligand index file. Notation: "-lg <Ligand group>", e.g. -lg 13"')
group.add_argument('-lt', dest='ligand_trajs', nargs='*', metavar='TRJ', type=gmx_trajectory,
                   help="Input trajectories of the unbound ligand for multiple trajectory approach. "
                        "Allowed formats: *.xtc (recommended), *.trr, *.pdb, *.gro (specify as many as you'd like).")
group.add_argument('-lp', dest='ligand_top', metavar='<Topology>', default=None, type=gmx_topology,
                   help="Topology file of the ligand.")

add_miscellaneous_actions(gmx_parser)
