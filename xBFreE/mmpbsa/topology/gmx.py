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

import platform
import string
import subprocess

import parmed
from xBFreE.exceptions import *
from xBFreE.utils.misc import log_subprocess_output
from xBFreE.utils.molecule import (list2range, res2map, get_indexes, check_str, eq_strs, get_index_groups)
from .core import BuildTop
from ..utils.changeradii import LoadRadii

echo_command = ['echo'] if platform.system() == "Darwin" else ['echo', '-e']

chains_letters = list(string.ascii_uppercase)


positive_aa = ['LYS', 'ARG', 'HIP']
negative_aa = ['GLU', 'ASP']
nonpolar_aa = ['PHE', 'TRP', 'VAL', 'ILE', 'LEU', 'MET', 'PRO', 'CYX', 'ALA', 'GLY']
polar_aa = ['TYR', 'SER', 'THR', 'CYM', 'CYS', 'HIE', 'HID', 'GLN', 'ASN', 'ASH', 'GLH', 'LYN']



class BuildTopGromacs(BuildTop):
    def __init__(self, FILES, INPUT, external_programs):
        super().__init__(FILES, INPUT, external_programs)
        self.use_temp = False
        self.com_mut_index = None

        # Define Gromacs executable
        self.make_ndx = [self.external_progs['gmx'], 'make_ndx']
        self.trjconv = [self.external_progs['gmx'], 'trjconv']
        self.editconf = [self.external_progs['gmx'], 'editconf']
        self.radii = LoadRadii(self.INPUT['general']['PBRadii'], self.INPUT['general']['radii_path'])
        self.checkFiles()

    def checkFiles(self):
        if (not self.FILES.complex_top and not self.FILES.complex_structure or not self.FILES.complex_index and
                not self.FILES.complex_trajs or not self.FILES.complex_groups):
            xBFreEErrorLogging('You must define the topology, structure, index and trajectories files, as well as the '
                               'groups!')

    def buildTopology(self):
        """
        :return: complex, receptor, ligand topologies and their mutants
        """

        self.gmx2pdb()
        tops = self.gmxtop2prmtop()
        # check if decomp or qmmm for residue selection
        self.decomp_qmmm_ressel()
        # FIXME: Is this step necessary? When the trajs are created with cpptraj a cleanup is made. However,
        #  this step get the complex trajectory in multi-component system
        self.cleanup_trajs()
        return tops

    def gmx2pdb(self):
        """
        Generate PDB file to generate topology
        # TODO: make a function to get water residues index for explicit water calculations
        :return:
        """

        logging.info('Get PDB files from GROMACS structures files...')

        # wt complex
        # make index for extract pdb structure
        com_rec_group, com_lig_group = self.FILES.complex_groups
        if com_rec_group == com_lig_group:
            xBFreEErrorLogging('The receptor and ligand groups have to be different')
        print(self.FILES.complex_index, com_rec_group)
        num_com_rec_group, str_com_rec_group = get_index_groups(self.FILES.complex_index, com_rec_group)
        num_com_lig_group, str_com_lig_group = get_index_groups(self.FILES.complex_index, com_lig_group)

        logging.info('Making gmx_MMPBSA index for complex...')
        # merge both (rec and lig) groups into complex group, modify index and create a copy
        # 1-rename groups, 2-merge
        make_ndx_echo_args = echo_command + ['name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | '
                                             '{l}\n q\n'.format(r=num_com_rec_group, l=num_com_lig_group)]
        c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

        com_ndx = 'COM_index.ndx'
        make_ndx_args = self.make_ndx + ['-n', self.FILES.complex_index, '-o', com_ndx, '-f',
                                         self.FILES.complex_structure]
        logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                      (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                      ' '.join(make_ndx_args))
        c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_subprocess_output(c2)
        if c2.wait():  # if it quits with return code != 0
            xBFreEErrorLogging('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.complex_index))
        self.FILES.complex_index = com_ndx

        logging.info(f'Normal Complex: Saving group {str_com_rec_group}_{str_com_lig_group} '
                     f'({num_com_rec_group}_{num_com_lig_group}) in {self.FILES.complex_index} file as '
                     f'{self.complex_str_file}')
        # avoid PBC and not chain ID problems
        pdbcom_echo_args = echo_command + ['GMXMMPBSA_REC_GMXMMPBSA_LIG']
        c3 = subprocess.Popen(pdbcom_echo_args, stdout=subprocess.PIPE)

        str_format = 'tpr' if self.FILES.complex_structure[-3:] == 'tpr' else 'pdb'
        if str_format == 'tpr':
            comprog = self.trjconv
            # we extract the pdb from the first frame of trajs to make amber topology
            pdbcom_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_structure, '-o',
                                          self.complex_str_file, '-n', self.FILES.complex_index, '-dump', '0']
        else:
            comprog = self.editconf
            pdbcom_args = self.editconf + ['-f', self.FILES.complex_structure, '-n', self.FILES.complex_index, '-o',
                                           self.complex_str_file]
        logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                      (' '.join(pdbcom_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                      '| ' + ' '.join(pdbcom_args))
        c4 = subprocess.Popen(pdbcom_args, stdin=c3.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_subprocess_output(c4)
        if c4.wait():  # if it quits with return code != 0
            xBFreEErrorLogging('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))

        # make a temp receptor pdb (even when stability) if decomp to get correct receptor residues from complex. This
        # avoids get multiples molecules from complex.split()
        if self.INPUT['decomp']['decomprun'] and self.FILES.stability:
            self.use_temp = True
            logging.warning('When &decomp is defined, we generate a receptor file in order to extract interface '
                            'residues')
            rec_echo_args = echo_command + ['{}'.format(num_com_rec_group)]
            cp1 = subprocess.Popen(rec_echo_args, stdout=subprocess.PIPE)
            if str_format == 'tpr':
                # we extract the pdb from the first frame of trajs to make amber topology
                pdbrec_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_structure, '-o',
                                              'rec_temp.pdb', '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdbrec_args = self.editconf + ['-f', self.FILES.complex_structure, '-n', self.FILES.complex_index, '-o',
                                               'rec_temp.pdb']
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(rec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=cp1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # check if stability
        if self.FILES.stability and (
                (self.FILES.receptor_tpr or self.FILES.ligand_tpr)
        ):
            logging.warning('When Stability calculation mode is selected, receptor and ligand files are not '
                            'needed...')
        # wt receptor
        if self.FILES.receptor_tpr:
            num_rec_group, str_rec_group = get_index_groups(self.FILES.receptor_index, self.FILES.receptor_group)

            logging.info('Making gmx_MMPBSA index for receptor...')
            make_ndx_echo_args = echo_command + ['name {r} GMXMMPBSA_REC\n q\n'.format(r=num_rec_group)]
            c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

            rec_ndx = 'REC_index.ndx'
            make_ndx_args = self.make_ndx + ['-n', self.FILES.receptor_index, '-o', rec_ndx, '-f',
                                             self.FILES.receptor_tpr]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                          ' '.join(make_ndx_args))
            c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c2)
            if c2.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.receptor_index))
            self.FILES.receptor_index = rec_ndx

            logging.info(f'Normal Receptor: Saving group {str_rec_group} ({num_rec_group}) in '
                         f'{self.FILES.receptor_index} file as {self.receptor_str_file}')
            pdbrec_echo_args = echo_command + ['{}'.format(num_rec_group)]
            p1 = subprocess.Popen(pdbrec_echo_args, stdout=subprocess.PIPE)
            str_format = 'tpr' if self.FILES.receptor_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                prog = self.trjconv
                # we extract a pdb from structure file to make amber topology
                pdbrec_args = self.trjconv + ['-f', self.FILES.receptor_trajs[0], '-s', self.FILES.receptor_tpr, '-o',
                                              self.receptor_str_file, '-n', self.FILES.receptor_index, '-dump', '0']
            else:
                prog = self.editconf
                pdbrec_args = self.editconf + ['-f', self.FILES.receptor_tpr, '-n', self.FILES.receptor_index, '-o',
                                               self.receptor_str_file]

            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdbrec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(prog), self.FILES.receptor_trajs[0]))
        else:
            logging.info(f'Normal Receptor: Saving group {str_com_rec_group} ({num_com_rec_group}) in '
                         f'{self.FILES.complex_index} file as {self.receptor_str_file}')
            pdbrec_echo_args = echo_command + ['{}'.format(num_com_rec_group)]
            cp1 = subprocess.Popen(pdbrec_echo_args, stdout=subprocess.PIPE)
            str_format = 'tpr' if self.FILES.complex_structure[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                # we extract a pdb from structure file to make amber topology
                pdbrec_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_structure, '-o',
                                              self.receptor_str_file, '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdbrec_args = self.editconf + ['-f', self.FILES.complex_structure, '-n', self.FILES.complex_index, '-o',
                                               self.receptor_str_file]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdbrec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=cp1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # ligand
        # # check consistence
        if self.FILES.ligand_tpr:  # ligand is protein
            logging.info('A ligand structure file was defined. Using MT approach...')
            num_lig_group, str_lig_group = get_index_groups(self.FILES.ligand_index, self.FILES.ligand_group)

            logging.info('Making gmx_MMPBSA index for ligand...')
            make_ndx_echo_args = echo_command + ['name {l} GMXMMPBSA_LIG\n q\n'.format(l=num_lig_group)]
            c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

            lig_ndx = self.mm + 'LIG_index.ndx'
            make_ndx_args = self.make_ndx + ['-n', self.FILES.ligand_index, '-o', lig_ndx, '-f', self.FILES.ligand_tpr]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                          ' '.join(make_ndx_args))
            c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c2)
            if c2.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.ligand_index))
            self.FILES.ligand_index = lig_ndx

            logging.info(f'Normal Ligand: Saving group {str_lig_group} ({num_lig_group}) in {self.FILES.ligand_index}'
                         f' file as {self.ligand_str_file}')
            # wt ligand
            pdblig_echo_args = echo_command + ['{}'.format(num_lig_group)]
            l1 = subprocess.Popen(pdblig_echo_args, stdout=subprocess.PIPE)
            str_format = 'tpr' if self.FILES.ligand_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                prog = self.trjconv
                # we extract a pdb from structure file to make amber topology
                pdblig_args = self.trjconv + ['-f', self.FILES.ligand_trajs[0], '-s', self.FILES.ligand_tpr, '-o',
                                              self.ligand_str_file, '-n', self.FILES.ligand_index, '-dump', '0']
            else:
                prog = self.editconf
                pdblig_args = self.editconf + ['-f', self.FILES.ligand_tpr, '-n', self.FILES.ligand_index, '-o',
                                               self.ligand_str_file]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdblig_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdblig_args))
            l2 = subprocess.Popen(pdblig_args, stdin=l1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(l2)
            if l2.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(prog), self.FILES.ligand_trajs[0]))
        else:
            # wt complex ligand
            logging.info(f'Normal Ligand: Saving group {str_com_lig_group} ({num_com_lig_group}) in '
                         f'{self.FILES.complex_index} file as {self.ligand_str_file}')
            pdblig_echo_args = echo_command + ['{}'.format(num_com_lig_group)]
            l1 = subprocess.Popen(pdblig_echo_args, stdout=subprocess.PIPE)

            str_format = 'tpr' if self.FILES.complex_structure[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                # we extract a pdb from structure file to make amber topology
                pdblig_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_structure, '-o',
                                              self.ligand_str_file, '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdblig_args = self.editconf + ['-f', self.FILES.complex_structure, '-n', self.FILES.complex_index, '-o',
                                               self.ligand_str_file]

            # we extract a pdb from structure file to make amber topology
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdblig_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdblig_args))
            l2 = subprocess.Popen(pdblig_args, stdin=l1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(l2)
            if l2.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # check for IE variable
        # FIXME: change it to ERROR?
        if (self.FILES.receptor_tpr or self.FILES.ligand_tpr) and (
                self.INPUT['general']['interaction_entropy'] or self.INPUT['general']['c2_entropy']
        ):
            logging.warning("The IE or C2 entropy method don't support the MT approach...")
            self.INPUT['general']['interaction_entropy'] = self.INPUT['general']['c2_entropy'] = 0

        # initialize receptor and ligand structures. Needed to get residues map
        self.complex_str = self.molstr(self.complex_str_file)
        self.receptor_str = self.molstr(self.receptor_str_file)
        self.ligand_str = self.molstr(self.ligand_str_file)
        if self.FILES.reference_structure:
            self.ref_str = check_str(self.FILES.reference_structure, ref=True)
        # FIXME: remove, since now the topology is processed with the index file
        # self.check4water()
        self.indexes = get_indexes(com_ndx=self.FILES.complex_index,
                                   rec_ndx=self.FILES.receptor_index,
                                   lig_ndx=self.FILES.ligand_index)
        self.resi, self.resl, self.orderl = res2map(self.indexes, self.complex_str)
        self.check_structures(self.complex_str, self.receptor_str, self.ligand_str)

    def gmxtop2prmtop(self):
        # FIXME: use Tan&Luo radii to avoid using radiopt variable (end in error when its not amber protein)
        # self.INPUT['pb']['radiopt'] = 0
        logging.info('Building Normal Amber topologies...')

        gmx_com_top = self.cleantop(self.FILES.complex_top, self.indexes['COM']['COM'])

        eq_strs(gmx_com_top, self.complex_str, molid='complex')

        gmx_com_top.coordinates = self.complex_str.coordinates
        gmx_com_top.save("COM.inpcrd", format='rst7', overwrite=True)
        # IMPORTANT: make_trajs ends in error if the box is defined
        gmx_com_top.box = None

        if gmx_com_top.impropers or gmx_com_top.urey_bradleys:
            top_class = parmed.amber.ChamberParm
        else:
            top_class = parmed.amber.AmberParm

        rec_indexes_string = ','.join(self.resi['REC']['string'])

        if self.FILES.receptor_top:
            gmx_rec_top = self.cleantop(self.FILES.receptor_top, self.indexes['REC'], 'receptor')
            eq_strs(gmx_rec_top, self.receptor_str, molid='receptor')
            gmx_rec_top.coordinates = self.receptor_str.coordinates
        else:
            # we make a copy for receptor topology
            gmx_rec_top = self.molstr(gmx_com_top)
            gmx_rec_top.strip(f'!:{rec_indexes_string}')

        gmx_rec_top.save("REC.inpcrd", format='rst7', overwrite=True)
        gmx_rec_top.box = None

        if self.FILES.ligand_top:
            gmx_lig_top = self.cleantop(self.FILES.ligand_top, self.indexes['LIG'], 'ligand')

            eq_strs(gmx_lig_top, self.ligand_str, molid='ligand')

            gmx_lig_top.coordinates = self.ligand_str.coordinates
        else:
            # we make a copy for ligand topology
            gmx_lig_top = self.molstr(gmx_com_top)
            gmx_lig_top.strip(f':{rec_indexes_string}')

        gmx_lig_top.save("LIG.inpcrd", format='rst7', overwrite=True)
        gmx_lig_top.box = None

        logging.info(f"Assigning PBRadii {self.INPUT['general']['PBRadii']} to Normal topologies...")
        self.radii.assign_radii(gmx_com_top)
        self.radii.assign_radii(gmx_rec_top)
        self.radii.assign_radii(gmx_lig_top)

        logging.info("Saving Normal Topology files...")
        com_amb_prm = top_class.from_structure(gmx_com_top)
        self.fixparm2amber(com_amb_prm)
        # IMPORTANT: In this case, we need to assign RADIUS_SET manually since GromacsTopologyFile don't contain it
        com_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        com_amb_prm.write_parm(self.complex_prmtop)
        rec_amb_prm = top_class.from_structure(gmx_rec_top)
        rec_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        self.fixparm2amber(rec_amb_prm)
        rec_amb_prm.write_parm(self.receptor_prmtop)
        lig_amb_prm = top_class.from_structure(gmx_lig_top)
        lig_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
        self.fixparm2amber(lig_amb_prm)
        lig_amb_prm.write_parm(self.ligand_prmtop)

        if self.INPUT['ala']['alarun']:
            logging.debug('Building Mutant Amber Topology...')
            # get mutation index in complex
            self.com_mut_index, self.part_mut, self.part_index = self.getMutationInfo()
            gmx_mut_com_top = self.makeMutTop(gmx_com_top, self.com_mut_index)
            gmx_mut_com_top.save("MUT_COM.inpcrd", format='rst7', overwrite=True)
            gmx_mut_com_top.box = None

            if self.part_mut == 'REC':
                logging.debug('Detecting mutation in Receptor. Building Mutant Receptor topology...')
                out_prmtop = self.mutant_receptor_prmtop
                self.mutant_ligand_prmtop = None
                if self.FILES.receptor_top:
                    mut_gmx_top = self.makeMutTop(gmx_rec_top, self.part_index)
                else:
                    mut_gmx_top = gmx_com_top.__copy__()
                    mut_gmx_top.strip(f'!:{rec_indexes_string}')
            else:
                logging.debug('Detecting mutation in Ligand. Building Mutant Ligand topology...')
                out_prmtop = self.mutant_ligand_prmtop
                self.mutant_receptor_prmtop = None
                if self.FILES.ligand_top:
                    mut_gmx_top = self.makeMutTop(gmx_lig_top, self.part_index)
                else:
                    mut_gmx_top = gmx_com_top.__copy__()
                    mut_gmx_top.strip(f':{rec_indexes_string}')
            mut_gmx_top.box = None

            logging.info(f"Assigning PBRadii {self.INPUT['general']['PBRadii']} to Mutant topologies...")
            self.radii.assign_radii(gmx_mut_com_top)
            self.radii.assign_radii(mut_gmx_top)

            logging.info("Saving Mutant Topology files...")
            mut_com_amb_prm = top_class.from_structure(gmx_mut_com_top)
            self.fixparm2amber(mut_com_amb_prm)
            mut_com_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
            mut_com_amb_prm.write_parm(self.mutant_complex_prmtop)
            # save receptor or ligand mutant
            mut_amb_prm = top_class.from_structure(mut_gmx_top)
            self.fixparm2amber(mut_amb_prm)
            mut_amb_prm.parm_data['RADIUS_SET'][0] = self.radii.radius_set_text
            mut_amb_prm.write_parm(out_prmtop)
        else:
            self.mutant_complex_prmtop = None

        return (self.complex_prmtop, self.receptor_prmtop, self.ligand_prmtop, self.mutant_complex_prmtop,
                self.mutant_receptor_prmtop, self.mutant_ligand_prmtop)

    def read_top(self, top_file, id='complex'):
        """
        Read the top file. In this function it will be extracted explicit water molecules. Not implemented yet.
        :param top_file:
        :param id:
        :return:
        """
        top = parmed.gromacs.GromacsTopologyFile(top_file)

        return top

    def cleantop(self, top_file, ndx, id='complex'):
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
        # top_file = Path(top_file)
        # molsect = False
        #
        # ttp_file = top_file.parent.joinpath('_temp_top.top')
        # temp_top = ttp_file.open(mode='w')
        # # temp_top.write('; Modified by gmx_MMPBSA\n')
        # # TODO: keep solvent when n-wat is implemented
        # with open(top_file) as topf:
        #     for line in topf:
        #         if '[ molecules ]' in line:
        #             molsect = True
        #         if molsect:
        #             # not copy ions and solvent
        #             sol_ion = [
        #                 # standard gmx form
        #                 'NA', 'CL', 'SOL', 'K'
        #                 # charmm-GUI form ??
        #                 'SOD', 'Na+', 'CLA', 'Cl-', 'POT', 'K+',
        #                 'TIP3P', 'TIP3', 'TP3', 'TIPS3P', 'TIP3o',
        #                 'TIP4P', 'TIP4PEW', 'T4E', 'TIP4PD',
        #                 'TIP5P',
        #                 'SPC', 'SPC/E', 'SPCE',
        #                 'WAT',
        #                 'OPC']
        #             if not line.split():
        #                 continue
        #             if line.split()[0].strip() in sol_ion:
        #                 continue
        #         temp_top.write(line)
        # temp_top.close()

        # read the topology with parmed
        top = parmed.gromacs.GromacsTopologyFile(top_file)

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

        for i in ndx:
            try:
                idx = top.atoms[i - 1].residue.idx + 1
                if idx not in res_list:
                    res_list.append(top.atoms[i - 1].residue.number + 1)
            except IndexError:
                xBFreEErrorLogging(f'The atom {i} in the {id} index is not found in the topology file. Please check that '
                                'the files are consistent.')

        ranges = list2range(res_list)
        top.strip(f"!:{','.join(ranges['string'])}")

        return top

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
        # clear trajectory
        if not self.INPUT['general']['solvated_trajectory']:
            return
        logging.info('Cleaning normal complex trajectories...')
        new_trajs = []
        for i in range(len(self.FILES.complex_trajs)):
            trjconv_echo_args = echo_command + ['GMXMMPBSA_REC_GMXMMPBSA_LIG']
            c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file and make amber topology for complex
            trjconv_args = self.trjconv + ['-f', self.FILES.complex_trajs[i], '-s', self.FILES.complex_structure, '-o',
                                           f'COM_traj_{i}.xtc', '-n', self.FILES.complex_index]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(trjconv_args))
            c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c6)
            if c6.wait():  # if it quits with return code != 0
                xBFreEErrorLogging('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.complex_trajs[i]))
            new_trajs.append(f'COM_traj_{i}.xtc')
        self.FILES.complex_trajs = new_trajs

        # clear trajectory
        if self.FILES.receptor_tpr:
            logging.info('Cleaning normal receptor trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.receptor_trajs)):
                trjconv_echo_args = echo_command + ['GMXMMPBSA_REC']
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = self.trjconv + ['-f', self.FILES.receptor_trajs[i], '-s', self.FILES.receptor_tpr,
                                               '-o', f'REC_traj_{i}.xtc', '-n', self.FILES.receptor_index]
                logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                              (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                              '| ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_subprocess_output(c6)
                if c6.wait():  # if it quits with return code != 0
                    xBFreEErrorLogging(
                        '%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.receptor_trajs[i]))
                new_trajs.append(f'REC_traj_{i}.xtc')
            self.FILES.receptor_trajs = new_trajs

        if self.FILES.ligand_tpr:
            logging.info('Cleaning normal ligand trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.ligand_trajs)):
                trjconv_echo_args = echo_command + ['GMXMMPBSA_LIG']
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = self.trjconv + ['-f', self.FILES.ligand_trajs[i], '-s', self.FILES.ligand_tpr, '-o',
                                               f'LIG_traj_{i}.xtc', '-n', self.FILES.ligand_index]
                logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                              (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                              '| ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_subprocess_output(c6)
                if c6.wait():  # if it quits with return code != 0
                    xBFreEErrorLogging('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.ligand_trajs[i]))
                new_trajs.append(f'LIG_traj_{i}.xtc')
            self.FILES.ligand_trajs = new_trajs

    def check_structures(self, com_str, rec_str=None, lig_str=None):
        logging.info('Checking the structures consistency...')
        check_str(com_str)
        check_str(rec_str, skip=True)
        check_str(lig_str, skip=True)

        if self.FILES.reference_structure:
            logging.info('Assigning chain ID to structures files according to the reference structure...')
            ref_str = check_str(self.FILES.reference_structure)
            if len(ref_str.residues) != len(com_str.residues):
                xBFreEErrorLogging(f'The number of residues of the complex ({len(com_str.residues)}) and of the '
                                f'reference structure ({len(ref_str.residues)}) are different. Please check that the '
                                f'reference structure is correct')
            for c, res in enumerate(ref_str.residues):
                if com_str.residues[c].number != res.number or com_str.residues[c].name != res.name:
                    xBFreEErrorLogging('There is no match between the complex and the reference structure used. An '
                                    f'attempt was made to assign the chain ID to "{com_str.residues[c].name}'
                                    f':{com_str.residues[c].number}:{com_str.residues[c].insertion_code}" in the '
                                    f'complex, but "{res.name}:{res.number}:{res.insertion_code}" was expected '
                                    'based on the reference structure. Please check that the reference structure is '
                                    'correct')
                com_str.residues[c].chain = res.chain
                i = self.resl[c].id_index - 1
                if self.resl[c].is_receptor():
                    rec_str.residues[i].chain = res.chain
                else:
                    lig_str.residues[i].chain = res.chain
        else:
            assign = False
            if self.INPUT['general']['assign_chainID'] == 1:
                assign = not com_str.residues[0].chain  # pretty simple
                if assign:
                    logging.info('Chains ID not found. Assigning chains IDs...')
                else:
                    logging.info('Chains ID found. Ignoring chains ID assignation...')
            elif self.INPUT['general']['assign_chainID'] == 2:
                assign = True
                if com_str.residues[0].chain:
                    logging.warning('Assigning chains ID...')
                else:
                    logging.warning('Already have chain ID. Re-assigning ID...')
            elif self.INPUT['general']['assign_chainID'] == 0 and not com_str.residues[0].chain:
                assign = True
                logging.warning('No reference structure was found and the complex structure not contain any chain ID. '
                                'Assigning chains ID automatically...')
            if assign:
                self._assign_chains_IDs(com_str, rec_str, lig_str)
        # Save fixed complex structure for analysis and set it in FILES to save in info file
        com_str.save('COM_FIXED.pdb', 'pdb', True, renumber=False)
        self.FILES.complex_fixed = 'COM_FIXED.pdb'
        logging.info('')

    def _assign_chains_IDs(self, com_str, rec_str, lig_str):
        chains_ids = []
        chain_by_num = False
        chain_by_ter = False
        previous_res_number = 0
        curr_chain_id = 'A'
        has_nucl = 0
        for c, res in enumerate(com_str.residues):
            if res.chain:
                if res.chain != curr_chain_id:
                    res.chain = curr_chain_id
                    i = self.resl[c].id_index - 1
                    if self.resl[c].is_receptor():
                        rec_str.residues[i].chain = res.chain
                    else:
                        lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            else:
                res.chain = curr_chain_id

                i = self.resl[c].id_index - 1
                if self.resl[c].is_receptor():
                    rec_str.residues[i].chain = res.chain
                else:
                    lig_str.residues[i].chain = res.chain
                if curr_chain_id not in chains_ids:
                    chains_ids.append(curr_chain_id)
                    # see if it is the end of chain
            if res.number != previous_res_number + 1 and previous_res_number != 0:
                chain_by_num = True
            if chain_by_num and chain_by_ter:
                chain_by_num = False
                chain_by_ter = False
                curr_chain_id = chains_letters[chains_letters.index(chains_ids[-1]) + 1]
                res.chain = curr_chain_id

                i = self.resl[c].id_index - 1
                if self.resl[c].is_receptor():
                    rec_str.residues[i].chain = res.chain
                else:
                    lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            elif chain_by_ter:
                chain_by_ter = False
            elif chain_by_num:
                chain_by_num = False
                curr_chain_id = chains_letters[chains_letters.index(chains_ids[-1]) + 1]
                res.chain = curr_chain_id
                i = self.resl[c].id_index - 1
                if self.resl[c + 1].is_receptor():
                    rec_str.residues[i].chain = res.chain
                else:
                    lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            for atm in res.atoms:
                if atm.name == 'OXT':  # only for protein
                    res.ter = True
                    chain_by_ter = True
            if parmed.residue.RNAResidue.has(res.name) or parmed.residue.DNAResidue.has(res.name):
                has_nucl += 1

            previous_res_number = res.number
        if has_nucl == 1:
            logging.warning('This structure contains nucleotides. We recommend that you use the reference structure')

