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
from pathlib import Path

import parmed
from parmed.amber import AmberMask
from xBFreE.exceptions import xBFreEErrorLogging
from xBFreE.utils.molecule import get_indexes_from_str, res2map, check_str, eq_strs
from .core import BuildTop
from ..utils.changeradii import LoadRadii


class BuildTopCHARMM(BuildTop):
    def __init__(self, FILES, INPUT, external_programs):
        super().__init__(FILES, INPUT, external_programs)
        self.use_temp = False
        self.com_mut_index = None
        self.radii = LoadRadii(self.INPUT['general']['PBRadii'], self.INPUT['general']['radii_path'])
        self.checkFiles()

        if (not self.FILES.complex_top and not self.FILES.complex_structure or
                not self.FILES.complex_trajs or not self.FILES.complex_groups):
            xBFreEErrorLogging('You must define the topology, structure and trajectories files, as well as the groups!')

    def buildTopology(self):
        """
        :return: complex, receptor, ligand topologies and their mutants
        """

        self.str2pdb()
        tops = self.psf2prmtop()
        # check if decomp or qmmm for residue selection
        self.decomp_qmmm_ressel()
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

    def psf2prmtop(self):

        # create complex, receptor and ligand PDBs
        com_rec_group, com_lig_group = self.FILES.complex_groups
        # create com topology
        namd_com_top = self.read_psf(self.FILES.complex_top)
        namd_com_top.strip(f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}")

        if namd_com_top.impropers or namd_com_top.urey_bradleys:
            top_class = parmed.amber.ChamberParm
        else:
            top_class = parmed.amber.AmberParm

        rec_indexes_string = ','.join(self.resi['REC']['string'])

        # create rec topology
        if self.FILES.receptor_top:
            namd_rec_top = self.read_psf(self.FILES.receptor_top)
            namd_rec_top.strip(f"!{self.FILES.receptor_group}")
        else:
            namd_rec_top = namd_com_top.__copy__()
            namd_rec_top.strip(f"!:{com_rec_group.strip(':')}")

        # create lig topology
        if self.FILES.receptor_top:
            namd_lig_top = self.read_psf(self.FILES.ligand_top)
            namd_lig_top.strip(f"!{self.FILES.ligand_group}")
        else:
            namd_lig_top = namd_com_top.__copy__()
            namd_lig_top.strip(f"!:{com_lig_group.strip(':')}")

        self.check_consistency([self.complex_str, namd_com_top], [self.receptor_str, namd_rec_top],
                               [self.ligand_str, namd_lig_top])


        namd_com_top.coordinates = self.complex_str.coordinates
        namd_com_top.save(f"{self.FILES.prefix}COM.inpcrd", format='rst7', overwrite=True)
        # IMPORTANT: make_trajs ends in error if the box is defined
        namd_com_top.box = None

        namd_rec_top.coordinates = self.receptor_str.coordinates
        namd_rec_top.save(f"{self.FILES.prefix}REC.inpcrd", format='rst7', overwrite=True)
        namd_rec_top.box = None

        namd_lig_top.coordinates = self.ligand_str.coordinates
        namd_lig_top.save(f"{self.FILES.prefix}LIG.inpcrd", format='rst7', overwrite=True)
        namd_lig_top.box = None

        logging.info(f"Assigning PBRadii {self.INPUT['general']['PBRadii']} to Normal topologies...")
        self.radii.assign_radii(namd_com_top)
        self.radii.assign_radii(namd_rec_top)
        self.radii.assign_radii(namd_lig_top)

        logging.info("Saving Normal Topology files...")
        com_amb_prm = top_class.from_structure(namd_com_top)
        self.fixparm2amber(com_amb_prm)
        # IMPORTANT: In this case, we need to assign RADIUS_SET manually since GromacsTopologyFile don't contain it
        com_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        com_amb_prm.write_parm(self.complex_pmrtop)
        rec_amb_prm = top_class.from_structure(namd_rec_top)
        rec_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        self.fixparm2amber(rec_amb_prm)
        rec_amb_prm.write_parm(self.receptor_pmrtop)
        lig_amb_prm = top_class.from_structure(namd_lig_top)
        lig_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        self.fixparm2amber(lig_amb_prm)
        lig_amb_prm.write_parm(self.ligand_pmrtop)

        # make the mutant
        if self.INPUT['ala']['alarun']:
            logging.info('Building Mutant Complex Topology...')
            # get mutation index in complex
            self.com_mut_index, self.part_mut, self.part_index = self.getMutationInfo()
            namd_mut_com_top = self.makeMutTop(namd_com_top, self.com_mut_index)
            namd_mut_com_top.save(f"{self.FILES.prefix}MUT_COM.inpcrd", format='rst7', overwrite=True)
            namd_mut_com_top.box = None

            if self.part_mut == 'REC':
                logging.debug('Detecting mutation in Receptor. Building Mutant Receptor topology...')
                out_prmtop = self.mutant_receptor_pmrtop
                self.mutant_ligand_pmrtop = None
                if self.FILES.receptor_top:
                    mut_namd_top = self.makeMutTop(namd_rec_top, self.part_index)
                else:
                    mut_namd_top = namd_com_top.__copy__()
                    mut_namd_top.strip(f'!:{rec_indexes_string}')
            else:
                logging.debug('Detecting mutation in Ligand. Building Mutant Ligand topology...')
                out_prmtop = self.mutant_ligand_pmrtop
                self.mutant_receptor_pmrtop = None
                if self.FILES.ligand_top:
                    mut_namd_top = self.makeMutTop(namd_lig_top, self.part_index)
                else:
                    mut_namd_top = namd_com_top.__copy__()
                    mut_namd_top.strip(f':{rec_indexes_string}')
            mut_namd_top.box = None

            logging.info(f"Assigning PBRadii {self.INPUT['general']['PBRadii']} to Mutant topologies...")
            self.radii.assign_radii(namd_mut_com_top)
            self.radii.assign_radii(mut_namd_top)

            logging.info("Saving Mutant Topology files...")
            mut_com_amb_prm = top_class.from_structure(namd_mut_com_top)
            self.fixparm2amber(mut_com_amb_prm)
            mut_com_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
            mut_com_amb_prm.write_parm(self.mutant_complex_pmrtop)
            mut_amb_prm = top_class.from_structure(mut_namd_top)
            self.fixparm2amber(mut_amb_prm)
            mut_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
            mut_amb_prm.write_parm(out_prmtop)
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)

    def read_psf(self, top_file, molid='complex'):
        # load topology
        namd_top = parmed.charmm.CharmmPsfFile(top_file)
        if molid == 'complex':
            toppar_name = ['complex_toppar', 'com_toppar', 'toppar']
        elif molid == 'receptor':
            toppar_name = ['receptor_toppar', 'rec_toppar']
        else:
            toppar_name = ['ligand_toppar', 'lig_toppar']

        file_list = []
        exist = False
        for name in toppar_name:
            toppar_folder = Path(top_file).parent.joinpath(name)
            if not toppar_folder.exists():
                continue
            file_list = [x.as_posix() for x in toppar_folder.glob('*')]
            exist = bool(len(file_list))
            break

        if not exist:
            xBFreEErrorLogging(f"The {top_file} associated toppar folder not exist. Please, define a valid one. Note "
                            f"that the toppar folder name can be defined as follow:\n"
                            f"COM: complex_toppar, com_toppar, toppar (for ST and MT)\n"
                            f"REC: receptor_toppar, rec_toppar (for MT only)\n"
                            f"LIG: ligand_toppar, lig_toppar (for MT only)\n")

        params = parmed.charmm.CharmmParameterSet(*file_list)
        namd_top.load_parameters(params)

        return namd_top

    def fixparm2amber(self, structure):
        his = ['HIS', 'HIE', 'HID', 'HIP']
        cys_name = ['CYS', 'CYX', 'CYM']

        for c, residue in enumerate(structure.residues, start=1):
            # change atoms name from GROMACS to AMBER
            for atom in residue.atoms:
                if atom.name == 'OC1':
                    atom.name = 'O'
                elif atom.name == 'OC2':
                    atom.name = 'OXT'
                    residue.ter = True  # parmed terminal
            # change residues name according to AMBER
            if residue.name == 'ILE':
                for atom in residue.atoms:
                    if atom.name == 'CD':
                        atom.name = 'CD1'
                        break
            elif residue.name == 'LYS':
                atoms = [atom.name for atom in residue.atoms]
                if 'HZ3' not in atoms:
                    residue.name = 'LYN'
            elif residue.name == 'ASP':
                atoms = [atom.name for atom in residue.atoms]
                if 'HD2' in atoms:
                    residue.name = 'ASH'
            elif residue.name == 'GLU':
                atoms = [atom.name for atom in residue.atoms]
                if 'HE2' in atoms:
                    residue.name = 'GLH'
            elif residue.name in his:
                atoms = [atom.name for atom in residue.atoms if atom.atomic_number == 1]
                if 'HD1' in atoms and 'HE2' in atoms:
                    residue.name = 'HIP'
                elif 'HD1' in atoms:
                    residue.name = 'HID'
                elif 'HE2' in atoms:
                    residue.name = 'HIE'
            elif residue.name in cys_name:
                for atom in residue.atoms:
                    if 'SG' in atom.name:
                        for bondedatm in atom.bond_partners:
                            if bondedatm.name == 'SG':
                                if residue.name == 'CYX' and bondedatm.residue.name == 'CYX':
                                    continue
                                residue.name = 'CYX'
                                bondedatm.residue.name = 'CYX'
                        break

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

        namd_com_top = self.read_psf(self.FILES.complex_top)
        if namd_com_top.impropers or namd_com_top.urey_bradleys:
            top_class = parmed.amber.ChamberParm
        else:
            top_class = parmed.amber.AmberParm
        com_amb_prm = top_class.from_structure(namd_com_top)
        new_trajs = []
        for i in range(len(self.FILES.complex_trajs)):
            trajsystem = Trajectory(self.FILES.complex_top, self.FILES.complex_trajs[i])
            trajsystem.Setup()
            # trajsystem.Image()
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
                # trajsystem.Image()
                trajsystem.Strip(f"!:{com_rec_group.strip(':')},{com_lig_group.strip(':')}")
                trajsystem.Outtraj(f'REC_traj_{i}.nc', filetype='nc')
                trajsystem.Run(f'{self.FILES.prefix}cleanup_trajs.out')
                new_trajs.append(f'REC_traj_{i}.nc')
            self.FILES.receptor_trajs = new_trajs

        if self.FILES.ligand_trajs:
            new_trajs = []
            for i in range(len(self.FILES.receptor_trajs)):
                trajsystem = Trajectory(com_amb_prm, self.FILES.complex_trajs[i])
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
