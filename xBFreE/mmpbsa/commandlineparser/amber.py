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

"""
This module contains classes and such that are responsible for parsing
command-line arguments for gmx_MMPBSA.  All the files specified for use
in gmx_MMPBSA will be assigned as attributes to the returned class.

"""

from xBFreE.mmpbsa.commandlineparser.parsers import (mmpbsa_parser,
                                                       amber_topology, amber_trajectory, pdb,
                                                       complex_group_des, receptor_group_des, ligand_group_des,
                                                       add_miscellaneous_actions)

amber_parser = mmpbsa_parser()

group = amber_parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cp', dest='complex_prmtop', metavar='<Topology File>', type=amber_topology, default=None,
                   help='''The complex Topology file''')
group.add_argument('-cg', dest='complex_groups', metavar='aa selection', nargs=2, default=None, type=str,
                   help='Selection for receptor and ligand in complex topology file. Both amber mask and residue '
                        'selection are allowed. The notation is as follows: '
                        '"-cg <Receptor group> <Ligand group>", ie. -cg :1-100 :101-120 (Amber mask) or '
                        'A/1-100 B/101-120 (Residue selection)')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', type=amber_trajectory,
                   help='''Complex trajectories. Make sure the trajectory is fitted and pbc have been removed. 
                   Allowed formats: *.mdcrd (recommended), *.crd, *.pdb, *.nc (specify as many as you'd like).''')
group.add_argument('-cr', dest='reference_structure', metavar='<PDB File>', default=None, type=pdb,
                   help='''Complex Reference Structure file. This option is optional but recommended (Use the PDB 
                   file used to generate the topology in Amber). If not defined, the chains ID assignment will be 
                   done automatically according to the structure (can generate wrong mapping).''')
group = amber_parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rp', dest='receptor_prmtop', metavar='<Topology File>', type=amber_topology,
                   help='''Topology file of the unbound receptor.''')
group.add_argument('-rg', dest='receptor_group', metavar='index', default=None, type=str,
                   help='''Receptor group in receptor topology file. Both amber mask  and residue selection are 
                   allowed. The notation is as follows: 
                   "-rg <Receptor group>", ie. -rg :1-100 (Amber mask) or A/1-100 (Residue selection)''')
group.add_argument('-rt', dest='receptor_trajs', nargs='*', metavar='TRJ', type=amber_trajectory,
                   help='''Input trajectories of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.mdcrd (recommended), *.crd, *.pdb, *.nc (specify as many
                         as you'd like).''')
group = amber_parser.add_argument_group('Ligand', ligand_group_des)
group.add_argument('-lp', dest='ligand_prmtop', metavar='<Topology File>', type=amber_topology,
                   help='''Topology file of the unbound ligand.''')
group.add_argument('-lg', dest='ligand_group', metavar='index', default=None, type=str,
                   help='''Ligand group in ligand topology file. Both amber mask '
                        'and residue selection are allowed. The notation is as follows: "-lg '
                        '<Ligand group>", ie. -lg :101-120 (Amber mask) or B/101-120 (Residue selection)''')
group.add_argument('-lt', dest='ligand_trajs', nargs='*', metavar='TRJ', type=amber_trajectory,
                   help='''Input trajectories of the unbound ligand for multiple trajectory approach.
                         Allowed formats: *.mdcrd (recommended), *.crd, *.pdb, *.nc (specify as many
                         as you'd like).''')

add_miscellaneous_actions(amber_parser)
#
# from icecream import ic
# ic('test----')
# amber_parser.parse_args(['-h'])
