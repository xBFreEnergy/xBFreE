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
This submodule contains the AMBER parser
"""

from xBFreE.mmpbsa.commandlineparser.parsers import (mmpbsa_parser,
                                                     amber_topology, amber_trajectory, pdb, mask,
                                                     complex_group_des, receptor_group_des, ligand_group_des,
                                                     add_miscellaneous_actions)

amber_parser = mmpbsa_parser()
group = amber_parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cs', dest='complex_structure', metavar='<Structure File>', default=None, type=pdb,
                   help='Structure file of the complex. Allowed formats: *.pdb')
group.add_argument('-cp', dest='complex_top', metavar='<Topology File>', type=amber_topology, default=None,
                   help='The complex Topology file')
group.add_argument('-cg', dest='complex_groups', metavar=('<receptor_mask>', '<ligand_mask>'), nargs=2, default=None,
                   type=mask,
                   help='Selection for receptor and ligand in complex topology file using amber mask. '
                        'The notation is as follows: '
                        '"-cg <receptor_mask> <ligand_mask>", ie. -cg :1-100,121-150 :101-120')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', type=amber_trajectory,
                   help="Complex trajectories. Make sure the trajectory is fitted and pbc have been removed. "
                        "Allowed formats:  *.nc, *.mdcrd, *.crd, *.pdb. (specify as many as you'd like)")

group = amber_parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rp', dest='receptor_top', metavar='<Topology File>', type=amber_topology,
                   help='Topology file of the unbound receptor.')
group.add_argument('-rg', dest='receptor_group', metavar='<receptor_mask>', default=None, type=mask,
                   help='Receptor group in receptor topology file using amber mask. '
                        'The notation is as follows: '
                        '"-rg <receptor_mask>", ie. -rg :1-100')
group.add_argument('-rt', dest='receptor_trajs', nargs='*', metavar='TRJ', type=amber_trajectory,
                   help='Input trajectories of the unbound receptor for multiple trajectory approach. '
                        "Allowed formats: *.nc, *.mdcrd, *.crd, *.pdb. (specify as many as you'd like).")

group = amber_parser.add_argument_group('Ligand', ligand_group_des)
group.add_argument('-lp', dest='ligand_top', metavar='<Topology File>', type=amber_topology,
                   help='Topology file of the unbound ligand.')
group.add_argument('-lg', dest='ligand_group', metavar='<ligand_mask>', default=None, type=mask,
                   help='Ligand group in ligand topology file using amber mask. '
                        'The notation is as follows: '
                        '"-lg <ligand_mask>", ie. -lg :1-20')
group.add_argument('-lt', dest='ligand_trajs', nargs='*', metavar='TRJ', type=amber_trajectory,
                   help='Input trajectories of the unbound ligand for multiple trajectory approach. '
                        "Allowed formats: *.nc, *.mdcrd, *.crd, *.pdb. (specify as many as you'd like).")

add_miscellaneous_actions(amber_parser)
