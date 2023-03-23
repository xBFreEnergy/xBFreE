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
import logging

from parmed.amber import AmberMask
from parmed.tools.changeradii import ChRad

from .core import BuildTop
import parmed
from xBFreE.exceptions import GMXMMPBSA_ERROR
from xBFreE.utils.molecule import list2range, get_indexes_from_str, res2map, check_str, eq_strs
from ..utils.changeradii import LoadRadii


class BuildTopAmber(BuildTop):
    def __init__(self, FILES, INPUT, external_programs):
        super().__init__(FILES, INPUT, external_programs)
        self.use_temp = False
        self.com_mut_index = None
        self.radii = LoadRadii(self.INPUT['general']['PBRadii'], self.INPUT['general']['radii_path'])

    def buildTopology(self):
        """
        :return: complex, receptor, ligand topologies and their mutants
        """

        self.str2pdb()
        tops = self.prmtop2prmtop()

        # FIXME: Is this step necessary? When the trajs are created with cpptraj a cleanup is made. However,
        #  this step get the complex trajectory in multi-component system
        self.cleanup_trajs()
        return tops

    def str2pdb(self):

        # create complex, receptor and ligand PDBs
        com_rec_group, com_lig_group = self.FILES.complex_groups
        # check pdb consistency
        check_str(self.FILES.complex_structure)

        # create PDBs. Here we use the complex_structure for get the receptor and ligand PBDs as well
        com_pdb = parmed.read_PDB(self.FILES.complex_structure)
        com_pdb.strip(f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}")
        rec_pdb = com_pdb.__copy__()
        rec_pdb.strip(f"!:{com_rec_group.strip(':')}")
        lig_pdb = com_pdb.__copy__()
        lig_pdb.strip(f"!:{com_lig_group.strip(':')}")

        # get amber selection based on AmberMask
        com_ind = AmberMask(com_pdb, f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}").Selection()
        com_rec_ind = AmberMask(com_pdb, f"!:{com_rec_group.strip(':')}").Selection()
        com_lig_ind = AmberMask(com_pdb, f"!:{com_lig_group.strip(':')}").Selection()

        rec_ind = (AmberMask(rec_pdb, f"!{self.FILES.receptor_group}").Selection() if self.FILES.receptor_top else None)
        lig_ind = (AmberMask(lig_pdb, f"!{self.FILES.ligand_group}").Selection() if self.FILES.ligand_top else None)

        # # initialize receptor and ligand structures. Needed to get residues map
        self.complex_str = self.molstr(com_pdb)
        self.complex_str.save(self.complex_str_file)
        self.receptor_str = self.molstr(rec_pdb)
        self.receptor_str.save(self.receptor_str_file)
        self.ligand_str = self.molstr(lig_pdb)
        self.ligand_str.save(self.ligand_str_file)

        # # FIXME: remove, since now the topology is processed with the index file
        # # self.check4water()
        self.indexes = get_indexes_from_str(com_ind={'com': com_ind, 'rec': com_rec_ind, 'lig': com_lig_ind},
                                            rec_ind=rec_ind, lig_ind=lig_ind)
        self.resi, self.resl, self.orderl = res2map(self.indexes, self.complex_str)

    def prmtop2prmtop(self):

        # create complex, receptor and ligand PDBs
        com_rec_group, com_lig_group = self.FILES.complex_groups
        # create com topology
        com_amb_prm = parmed.amber.AmberParm(self.FILES.complex_top)
        com_amb_prm.strip(f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}")

        # create rec topology
        if self.FILES.receptor_top:
            rec_amb_prm = parmed.amber.AmberParm(self.FILES.receptor_top)
            rec_amb_prm.strip(f"!{self.FILES.receptor_group}")
        else:
            rec_amb_prm = com_amb_prm.__copy__()
            rec_amb_prm.strip(f"!:{com_rec_group.strip(':')}")

        # create lig topology
        if self.FILES.receptor_top:
            lig_amb_prm = parmed.amber.AmberParm(self.FILES.ligand_top)
            lig_amb_prm.strip(f"!{self.FILES.ligand_group}")
        else:
            lig_amb_prm = com_amb_prm.__copy__()
            lig_amb_prm.strip(f"!:{com_lig_group.strip(':')}")

        self.check_consistency([self.complex_str, com_amb_prm], [self.receptor_str, rec_amb_prm],
                               [self.ligand_str, lig_amb_prm])


        com_amb_prm.coordinates = self.complex_str.coordinates
        com_amb_prm.save(f"{self.FILES.prefix}COM.inpcrd", format='rst7', overwrite=True)
        # IMPORTANT: make_trajs ends in error if the box is defined
        com_amb_prm.box = None

        rec_amb_prm.coordinates = self.receptor_str.coordinates
        rec_amb_prm.save(f"{self.FILES.prefix}REC.inpcrd", format='rst7', overwrite=True)
        rec_amb_prm.box = None

        lig_amb_prm.coordinates = self.ligand_str.coordinates
        lig_amb_prm.save(f"{self.FILES.prefix}LIG.inpcrd", format='rst7', overwrite=True)
        lig_amb_prm.box = None

        logging.info(f"Assigning PBRadii {self.INPUT['general']['PBRadii']} to Normal topologies...")
        self.radii.assign_radii(com_amb_prm)
        com_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        self.radii.assign_radii(rec_amb_prm)
        rec_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        self.radii.assign_radii(lig_amb_prm)
        lig_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text

        logging.info("Saving Normal Topology files...")
        com_amb_prm.write_parm(self.complex_pmrtop)
        rec_amb_prm.write_parm(self.receptor_pmrtop)
        lig_amb_prm.write_parm(self.ligand_pmrtop)

        # make the mutant
        if self.INPUT['ala']['alarun']:
            logging.info('Building Mutant Complex Topology...')
            # get mutation index in complex
            # FIXME: change to com_top from gromacs
            self.com_mut_index, self.part_mut, self.part_index = self.getMutationInfo()
            mut_com_amb_prm = self.makeMutTop(com_amb_prm, self.com_mut_index)
            mut_com_amb_prm.save(f"{self.FILES.prefix}MUT_COM.inpcrd", format='rst7', overwrite=True)
            mut_com_amb_prm.box = None

            if self.part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor topology...')
                out_prmtop = self.mutant_receptor_pmrtop
                self.mutant_ligand_pmrtop = None
                mut_amb_prm = self.makeMutTop(rec_amb_prm, self.part_index)
            else:
                logging.info('Detecting mutation in Ligand. Building Mutant Ligand topology...')
                out_prmtop = self.mutant_ligand_pmrtop
                self.mutant_receptor_pmrtop = None
                mut_amb_prm = self.makeMutTop(lig_amb_prm, self.part_index)

            mut_amb_prm.box = None

            logging.info(f"Assigning PBRadii {self.INPUT['general']['PBRadii']} to Mutant topologies...")
            self.radii.assign_radii(mut_com_amb_prm)
            mut_com_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
            mut_com_amb_prm.write_parm(self.mutant_complex_pmrtop)
            self.radii.assign_radii(mut_amb_prm)
            mut_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
            mut_amb_prm.write_parm(out_prmtop)
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)


    def cleantop(self, top_file, masks, id='complex'):
        """
        Create a new top instance
        - get defined molecules
            - remove ions
            - remove water (check is explicit_water was defined)
            - remove non-selected molecules

        :param top_file: User-defined topology file
        :param ndx: atoms index
        :return: new and clean top instance
        """

        # read the topology with parmed
        top = parmed.amber.AmberParm(top_file)

        # # get water residues index (starting for 1 like amber) when explicit_water is defined
        # water_index = []
        # if self.INPUT['general']['explicit_water']:
        #     for c, r in enumerate(top.residues):
        #         if r.name in ['TIP3P', 'TIP3', 'TP3', 'TIPS3P', 'TIP3o', 'TIP4P', 'TIP4PEW', 'T4E', 'TIP4PD', 'TIP5P',
        #                       'SPC', 'SPC/E', 'SPCE', 'WAT', 'OPC', 'HOH']:
        #             if c > self.INPUT['general']['explicit_water']:
        #                 break
        #             water_index.append(c)

        # get the residues in the top from the com_ndx
        res_list = []

        for i in masks:
            try:
                idx = top.atoms[i - 1].residue.idx + 1
                if idx not in res_list:
                    res_list.append(top.atoms[i - 1].residue.number + 1)
            except IndexError:
                GMXMMPBSA_ERROR(f'The atom {i} in the {id} index is not found in the topology file. Please check that '
                                'the files are consistent.')

        ranges = list2range(res_list)
        top.strip(f"!:{','.join(ranges['string'])}")

        return top


    def cleanup_trajs(self):
        """
        Generate dry/clean trajectories removing anything not selected as part of the system
        :return:
        """
        # clear trajectory
        if not self.INPUT['general']['solvated_trajectory']:
            return
        logging.info('Cleaning normal trajectories...')
        from xBFreE.mmpbsa.make_trajs import Trajectory

        com_rec_group, com_lig_group = self.FILES.complex_groups

        # FIXME: include explicit water

        new_trajs = []
        for i in range(len(self.FILES.complex_trajs)):
            trajsystem = Trajectory(self.FILES.complex_top, self.FILES.complex_trajs[i])
            trajsystem.Setup()
            trajsystem.Image()
            trajsystem.Strip(f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}")
            trajsystem.Outtraj(f'COM_traj_{i}.nc', filetype='nc')
            trajsystem.Run(f'{self.FILES.prefix}cleanup_trajs.out')
            new_trajs.append(f'COM_traj_{i}.nc')
        self.FILES.complex_trajs = new_trajs

        # clear trajectory
        if self.FILES.receptor_trajs:
            new_trajs = []
            for i in range(len(self.FILES.receptor_trajs)):
                trajsystem = Trajectory(self.FILES.complex_top, self.FILES.complex_trajs[i])
                trajsystem.Setup()
                trajsystem.Image()
                trajsystem.Strip(f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}")
                trajsystem.Outtraj(f'REC_traj_{i}.nc', filetype='nc')
                trajsystem.Run(f'{self.FILES.prefix}cleanup_trajs.out')
                new_trajs.append(f'REC_traj_{i}.nc')
            self.FILES.receptor_trajs = new_trajs

        if self.FILES.ligand_trajs:
            new_trajs = []
            for i in range(len(self.FILES.receptor_trajs)):
                trajsystem = Trajectory(self.FILES.complex_top, self.FILES.complex_trajs[i])
                trajsystem.Setup()
                trajsystem.Image()
                trajsystem.Strip(f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}")
                trajsystem.Outtraj(f'LIG_traj_{i}.nc', filetype='nc')
                trajsystem.Run(f'{self.FILES.prefix}cleanup_trajs.out')
                new_trajs.append(f'LIG_traj_{i}.nc')
            self.FILES.ligand_trajs = new_trajs

    def check_consistency(self, com_data, rec_data=None, lig_data=None):
        logging.info('Checking the structures consistency...')
        com_str, com_top = com_data
        rec_str, rec_top = rec_data or [None, None]
        lig_str, lig_top = lig_data or [None, None]
        check_str(com_str)
        if len(com_top.residues) != len(com_str.residues):
            raise

        eq_strs(com_str, com_top, molid='complex')

        if rec_data:
            check_str(rec_str, skip=True)
            eq_strs(rec_str, rec_top, molid='receptor')
            if len(rec_top.residues) != len(rec_str.residues):
                raise
        if lig_data:
            check_str(lig_str, skip=True)
            eq_strs(lig_str, lig_top, molid='ligand')
            if len(lig_top.residues) != len(lig_str.residues):
                raise

        # Save fixed complex structure for analysis and set it in FILES to save in info file
        com_str.save(f'{self.FILES.prefix}COM_FIXED.pdb', 'pdb', True, renumber=False)
        logging.info('')
