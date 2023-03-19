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
import textwrap

import parmed
from xBFreE.exceptions import *
from xBFreE.mmpbsa.utils.molecule import (list2range, res2map, get_indexes, check_str, eq_strs, get_index_groups)
from xBFreE.mmpbsa.utils.misc import log_subprocess_output
from xBFreE.mmpbsa.alamdcrd import _scaledistance
import subprocess
import logging
import string
from parmed.tools.changeradii import ChRad


if platform.system() == "Darwin":
    echo_command = ['echo']
else:
    echo_command = ['echo', '-e']

chains_letters = list(string.ascii_uppercase)
his = ['HIS', 'HIE', 'HID', 'HIP']
cys_name = ['CYS', 'CYX', 'CYM']
lys = ['LYS', 'LYN']
asp = ['ASP', 'ASH']
glu = ['GLU', 'GLH']
positive_aa = ['LYS', 'ARG', 'HIP']
negative_aa = ['GLU', 'ASP']
nonpolar_aa = ['PHE', 'TRP', 'VAL', 'ILE', 'LEU', 'MET', 'PRO', 'CYX', 'ALA', 'GLY']
polar_aa = ['TYR', 'SER', 'THR', 'CYM', 'CYS', 'HIE', 'HID', 'GLN', 'ASN', 'ASH', 'GLH', 'LYN']

PBRadii = {1: 'bondi', 2: 'mbondi', 3: 'mbondi2', 4: 'mbondi3', 5: 'mbondi_pb2', 6: 'mbondi_pb3', 7: 'charmm_radii'}

ions_para_files = {1: 'frcmod.ions234lm_126_tip3p', 2: 'frcmod.ions234lm_126_spce', 3: 'frcmod.ions234lm_126_tip4pew',
                   4: 'frcmod.ions234lm_hfe_tip3p', 5: 'frcmod.ions234lm_hfe_spce', 6: 'frcmod.ions234lm_hfe_tip4pew',
                   7: 'frcmod.ions234lm_iod_tip3p', 8: 'frcmod.ions234lm_iod_spce', 9: 'frcmod.ions234lm_iod_tip4pew',
                   10: 'frcmod.ionslm_126_opc', 11: 'frcmod.ionslm_hfe_opc', 12: 'frcmod.ionslm_iod_opc',
                   13: 'frcmod.ions1lm_126_tip3p', 14: 'frcmod.ions1lm_126_spce', 15: 'frcmod.ions1lm_126_tip4pew',
                   16: 'frcmod.ions1lm_iod'}

ions = ["AG", "AL", "Ag", "BA", "BR", "Be", "CA", "CD", "CE", "CL", "CO", "CR", "CS", "CU", "CU1", "Ce", "Cl-", "Cr",
        "Dy", "EU", "EU3", "Er", "F", "FE", "FE2", "GD3", "H3O+", "HE+", "HG", "HZ+", "Hf", "IN", "IOD", "K", "K+",
        "LA", "LI", "LU", "MG", "MN", "NA", "NH4", "NI", "Na+", "Nd", "PB", "PD", "PR", "PT", "Pu", "RB", "Ra", "SM",
        "SR", "Sm", "Sn", "TB", "TL", "Th", "Tl", "Tm", "U4+", "V2+", "Y", "YB2", "ZN", "Zr"]


class BuildAmberTopology:
    def __init__(self, prog, FILES, INPUT, external_programs):
        self.md_prog = prog
        self.FILES = FILES
        self.INPUT = INPUT
        self.external_progs = external_programs
        self.use_temp = False
        self.com_mut_index = None

        # Define Gromacs executable
        self.make_ndx = self.external_progs['make_ndx']
        self.trjconv = self.external_progs['trjconv']
        self.editconf = self.external_progs['editconf']

        self.ref_str = None

        # create the * prmtop variables for compatibility with the original code
        self.complex_pmrtop = 'COM.prmtop'
        self.receptor_pmrtop = 'REC.prmtop'
        self.ligand_pmrtop = 'LIG.prmtop'

        self.mutant_complex_pmrtop = 'MUT_COM.prmtop'
        self.mutant_receptor_pmrtop = 'MUT_REC.prmtop'
        self.mutant_ligand_pmrtop = 'MUT_LIG.prmtop'

        self.complex_str_file = f'{self.FILES.prefix}COM.pdb'
        self.receptor_str_file = f'{self.FILES.prefix}REC.pdb'
        self.ligand_str_file = f'{self.FILES.prefix}LIG.pdb'

        self.checkFiles()

    def checkFiles(self):
        if (not self.FILES.complex_tpr or not self.FILES.complex_index or
                not self.FILES.complex_trajs or not self.FILES.complex_groups):
            GMXMMPBSA_ERROR('You must define the structure, topology and index files, as well as the groups!')

    def build_from_gmx(self):
        pass

    def build_from_amber(self):
        pass

    def build_from_namd(self):
        pass

    def build_from_charmm(self):
        pass


    def buildTopology(self):
        """
        :return: complex, receptor, ligand topologies and their mutants
        """
        self.gmx2pdb()
        tops = self.gmxtop2prmtop()


        # FIXME: move on to a new function to be applied for all tops
        if self.INPUT['decomp']['decomprun']:
            decomp_res = self.get_selected_residues(self.INPUT['decomp']['print_res'])
            if 'within' in self.INPUT['decomp']['print_res']:
                if len(decomp_res) < 2:
                    logging.warning(f"Number of decomp residues to print using "
                                    f"print_res = '{self.INPUT['decomp']['print_res']}' < 2")
                    logging.info(
                        'Increasing cutoff value by 0.1 until number of decomp residues to print >= 2'
                    )
                    cutoff = float(self.INPUT['decomp']['print_res'].split()[1])
                    it = 0
                    while len(decomp_res) < 2:
                        cutoff = round(cutoff, 1) + 0.25
                        decomp_res = self.get_selected_residues(f'within {cutoff}')
                        if it == 20:
                            # probably not needed, but...
                            GMXMMPBSA_ERROR('The maximum number of iterations to select interaction residues was '
                                            'reached. Please set print_res with a valid selection.')
                        it += 1

                    logging.info(f"Selecting residues by distance ({round(cutoff, 1)} Å) between "
                                 f"receptor and ligand for decomposition analysis...")
                else:
                    logging.info(f"Selecting residues by distance ({self.INPUT['decomp']['print_res'].split()[1]} Å) between "
                                 f"receptor and ligand for decomposition analysis...")
            elif self.INPUT['decomp']['print_res'] == 'all':
                logging.info('Selecting all residues for decomposition analysis...')
            else:
                logging.info('User-selected residues for decomposition analysis...')

            textwraped = textwrap.wrap('\t'.join(x.string for x in decomp_res), tabsize=4, width=120)
            logging.info(f'Selected {len(decomp_res)} residues:\n' + '\n'.join(textwraped) + '\n')

            if self.INPUT['decomp']['idecomp'] in [3, 4]:
                if self.INPUT['decomp']['dec_verbose'] == 0:
                    mol_terms = 1
                elif self.INPUT['decomp']['dec_verbose'] == 1:
                    mol_terms = 3
                elif self.INPUT['decomp']['dec_verbose'] == 2:
                    mol_terms = 4
                else:
                    mol_terms = 12
                energy_terms = 6
                num_res = len(decomp_res)
                total_items = energy_terms * mol_terms * num_res ** 2
                if total_items > 250:
                    logging.warning(f"Using idecomp = {self.INPUT['decomp']['idecomp']} and dec_verbose ="
                                    f" {self.INPUT['decomp']['dec_verbose']} will generate approximately {total_items} items. "
                                    f"Large print selections demand a large amount of memory and take a "
                                    f"significant amount of time to print!")

            self.INPUT['decomp']['print_res'] = ','.join(list2range(decomp_res)['string'])
        if self.INPUT['gb']['ifqnt']:
            qm_residues, (rec_charge, lig_charge) = self.get_selected_residues(self.INPUT['gb']['qm_residues'], True)

            if 'within' in self.INPUT['gb']['qm_residues']:
                if len(qm_residues) == 0:
                    logging.warning(f"Number of qm_residues using print_res = '{self.INPUT['gb']['qm_residues']}' = 0")
                    logging.info(
                        'Increasing cutoff value by 0.1 until number of qm_residues > 0'
                    )
                    cutoff = float(self.INPUT['gb']['qm_residues'].split()[1])
                    it = 0
                    while len(qm_residues) == 0:
                        cutoff = round(cutoff, 1) + 0.25
                        qm_residues, (rec_charge, lig_charge) = self.get_selected_residues(f'within {cutoff}', True)
                        if it == 20:
                            # probably not needed, but...
                            GMXMMPBSA_ERROR('The maximum number of iterations to select interaction residues was '
                                            'reached. Please set print_res with a valid selection.')
                        it += 1

                    logging.info(f"Selecting residues by distance ({round(cutoff, 1)} Å) between "
                                 f"receptor and ligand for QM/MM calculation...")
                else:
                    logging.info(f"Selecting residues by distance ({self.INPUT['gb']['qm_residues'].split()[1]} Å) between "
                                 f"receptor and ligand for QM calculation...")
            elif self.INPUT['gb']['qm_residues'] == 'all':
                logging.info('Selecting all residues for QM calculation...')
            else:
                logging.info('User-selected residues for QM calculation...')

            textwraped = textwrap.wrap('\t'.join(x.string for x in qm_residues), tabsize=4, width=120)
            logging.info(f'Selected {len(qm_residues)} residues:\n' + '\n'.join(textwraped) + '\n')
            self.INPUT['gb']['qm_residues'] = ','.join(list2range(qm_residues)['string'])

            if self.INPUT['gb']['qmcharge_com'] != rec_charge + lig_charge:
                logging.warning('System specified with odd number of electrons. Most likely the charge of QM region '
                                '(qmcharge_com) have been set incorrectly.')
                self.INPUT['gb']['qmcharge_com'] = rec_charge + lig_charge
                logging.warning(f'Setting qmcharge_com = {rec_charge + lig_charge}')

            if self.INPUT['gb']['qmcharge_rec'] != rec_charge:
                logging.warning(f'Setting qmcharge_rec = {rec_charge}')
                self.INPUT['gb']['qmcharge_rec'] = rec_charge
            if self.INPUT['gb']['qmcharge_lig'] != lig_charge:
                logging.warning(f'Setting qmcharge_lig = {lig_charge}')
                self.INPUT['gb']['qmcharge_lig'] = lig_charge

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
            GMXMMPBSA_ERROR('The receptor and ligand groups have to be different')
        num_com_rec_group, str_com_rec_group = get_index_groups(self.FILES.complex_index, com_rec_group)
        num_com_lig_group, str_com_lig_group = get_index_groups(self.FILES.complex_index, com_lig_group)

        logging.info('Making gmx_MMPBSA index for complex...')
        # merge both (rec and lig) groups into complex group, modify index and create a copy
        # 1-rename groups, 2-merge
        make_ndx_echo_args = echo_command + ['name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | '
                                             '{l}\n q\n'.format(r=num_com_rec_group, l=num_com_lig_group)]
        c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

        com_ndx = self.FILES.prefix + 'COM_index.ndx'
        make_ndx_args = self.make_ndx + ['-n', self.FILES.complex_index, '-o', com_ndx, '-f', self.FILES.complex_tpr]
        logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                      (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                      ' '.join(make_ndx_args))
        c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_subprocess_output(c2)
        if c2.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.complex_index))
        self.FILES.complex_index = com_ndx

        logging.info(f'Normal Complex: Saving group {str_com_rec_group}_{str_com_lig_group} '
                     f'({num_com_rec_group}_{num_com_lig_group}) in {self.FILES.complex_index} file as '
                     f'{self.complex_str_file}')
        # avoid PBC and not chain ID problems
        pdbcom_echo_args = echo_command + ['GMXMMPBSA_REC_GMXMMPBSA_LIG']
        c3 = subprocess.Popen(pdbcom_echo_args, stdout=subprocess.PIPE)

        str_format = 'tpr' if self.FILES.complex_tpr[-3:] == 'tpr' else 'pdb'
        if str_format == 'tpr':
            comprog = self.trjconv
            # we extract the pdb from the first frame of trajs to make amber topology
            pdbcom_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                          self.complex_str_file, '-n', self.FILES.complex_index, '-dump', '0']
        else:
            comprog = self.editconf
            pdbcom_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                           self.complex_str_file]
        logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                      (' '.join(pdbcom_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                      '| ' + ' '.join(pdbcom_args))
        c4 = subprocess.Popen(pdbcom_args, stdin=c3.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_subprocess_output(c4)
        if c4.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))

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
                pdbrec_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                              'rec_temp.pdb', '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdbrec_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                               'rec_temp.pdb']
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(rec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=cp1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # check if stability
        if self.FILES.stability and (
                (self.FILES.receptor_tpr or self.FILES.ligand_tpr)
        ):
            logging.warning('When Stability calculation mode is selected, receptor and ligand files are not '
                            'needed...')
        # wt receptor
        if self.FILES.receptor_tpr:
            logging.info('A receptor structure file was defined. Using MT approach...')
            num_rec_group, str_rec_group = get_index_groups(self.FILES.receptor_index, self.FILES.receptor_group)

            logging.info('Making gmx_MMPBSA index for receptor...')
            make_ndx_echo_args = echo_command + ['name {r} GMXMMPBSA_REC\n q\n'.format(r=num_rec_group)]
            c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

            rec_ndx = self.FILES.prefix + 'REC_index.ndx'
            make_ndx_args = self.make_ndx + ['-n', self.FILES.receptor_index, '-o', rec_ndx, '-f',
                                             self.FILES.receptor_tpr]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                          ' '.join(make_ndx_args))
            c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c2)
            if c2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.receptor_index))
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
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(prog), self.FILES.receptor_trajs[0]))
        else:
            logging.info('No receptor structure file was defined. Using ST approach...')
            logging.info('Using receptor structure from complex to generate AMBER topology')
            logging.info(f'Normal Receptor: Saving group {str_com_rec_group} ({num_com_rec_group}) in '
                         f'{self.FILES.complex_index} file as {self.receptor_str_file}')
            pdbrec_echo_args = echo_command + ['{}'.format(num_com_rec_group)]
            cp1 = subprocess.Popen(pdbrec_echo_args, stdout=subprocess.PIPE)
            str_format = 'tpr' if self.FILES.complex_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                # we extract a pdb from structure file to make amber topology
                pdbrec_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                              self.receptor_str_file, '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdbrec_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                               self.receptor_str_file]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdbrec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=cp1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # ligand
        # # check consistence
        if self.FILES.ligand_tpr:  # ligand is protein
            logging.info('A ligand structure file was defined. Using MT approach...')
            num_lig_group, str_lig_group = get_index_groups(self.FILES.ligand_index, self.FILES.ligand_group)

            logging.info('Making gmx_MMPBSA index for ligand...')
            make_ndx_echo_args = echo_command + ['name {l} GMXMMPBSA_LIG\n q\n'.format(l=num_lig_group)]
            c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

            lig_ndx = self.FILES.prefix + 'LIG_index.ndx'
            make_ndx_args = self.make_ndx + ['-n', self.FILES.ligand_index, '-o', lig_ndx, '-f', self.FILES.ligand_tpr]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                          ' '.join(make_ndx_args))
            c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c2)
            if c2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.ligand_index))
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
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(prog), self.FILES.ligand_trajs[0]))
        else:
            # wt complex ligand
            logging.info('No ligand structure file was defined. Using ST approach...')
            logging.info('Using ligand structure from complex to generate AMBER topology')
            logging.info(f'Normal Ligand: Saving group {str_com_lig_group} ({num_com_lig_group}) in '
                         f'{self.FILES.complex_index} file as {self.ligand_str_file}')
            pdblig_echo_args = echo_command + ['{}'.format(num_com_lig_group)]
            l1 = subprocess.Popen(pdblig_echo_args, stdout=subprocess.PIPE)

            str_format = 'tpr' if self.FILES.complex_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                # we extract a pdb from structure file to make amber topology
                pdblig_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                              self.ligand_str_file, '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdblig_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                               self.ligand_str_file]

            # we extract a pdb from structure file to make amber topology
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdblig_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdblig_args))
            l2 = subprocess.Popen(pdblig_args, stdin=l1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(l2)
            if l2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # check for IE variable
        # FIXME: change it to ERROR?
        if (self.FILES.receptor_tpr or self.FILES.ligand_tpr) and (
                self.INPUT['general']['interaction_entropy'] or self.INPUT['general']['c2_entropy']
        ):
            logging.warning("The IE or C2 entropy method don't support the MTP approach...")
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

    def check4water(self):
        if counter := sum(
                res.name
                in [
                    'SOD', 'Na+', 'NA', 'Na', 'CLA', 'Cl-', 'CL', 'Cl', 'POT', 'K+', 'K',
                    'SOL', 'WAT',
                    'TIP3P', 'TIP3', 'TP3', 'TIPS3P', 'TIP3o',
                    'TIP3P', 'TIP3', 'TP3', 'TIPS3P', 'TIP3o',
                    'TIP4P', 'TIP4PEW', 'T4E', 'TIP4PD',
                    'TIP5P',
                    'SPC', 'SPCE',
                    'OPC'
                ]
                for res in self.complex_str.residues
        ):
            GMXMMPBSA_ERROR(f'gmx_MMPBSA does not support water/ions molecules in any structure, but we found'
                            f' {counter} molecules in the complex.')

    def _check_periodicity(self, parm, system):
        """
        check for periodicities == 0 and change them to 1. This is required especially for nmode calculations
        """
        invalid_per = 0

        for dt in parm.dihedral_types:
            if dt.per == 0:
                invalid_per += 1
                dt.per = 1

        if invalid_per:
            logging.warning(f'{invalid_per} invalid DIHEDRAL_PERIODICITY = 0 found in {system.capitalize()} '
                            f'topology... Setting DIHEDRAL_PERIODICITY = 1')

        return parm

    def gmxtop2prmtop(self):
        # FIXME: use Tan&Luo radii to avoid using radiopt variable (end in error when its not amber protein)
        logging.info('Using topology conversion. Setting radiopt = 0...')
        self.INPUT['pb']['radiopt'] = 0
        logging.info('Building Normal Complex Amber topology...')
        com_top = self.cleantop(self.FILES.complex_top, self.indexes['COM']['COM'])
        if error_info := eq_strs(com_top, self.complex_str):
            if error_info[0] == 'atoms':
                GMXMMPBSA_ERROR(f"The number of atoms in the topology ({error_info[1]}) and the complex structure "
                                f"({error_info[2]}) are different. Please check these files and verify that they are "
                                f"correct. Otherwise report the error...")
            else:
                GMXMMPBSA_ERROR(f"The number of residues in the topology ({error_info[1]}) and the complex structure "
                                f"({error_info[2]}) are different. Please check these files and verify that they are "
                                f"correct. Otherwise report the error...")

        com_top.coordinates = self.complex_str.coordinates
        com_top.save(f"{self.FILES.prefix}COM.inpcrd", format='rst7', overwrite=True)
        # try:
        if com_top.impropers or com_top.urey_bradleys:
            com_amb_prm = parmed.amber.ChamberParm.from_structure(com_top)
            com_top_parm = 'chamber'
            logging.info('Detected CHARMM force field topology format...')
        else:
            com_amb_prm = parmed.amber.AmberParm.from_structure(com_top)
            com_top_parm = 'amber'
            logging.info('Detected Amber/OPLS force field topology format...')

        # IMPORTANT: make_trajs ends in error if the box is defined
        com_amb_prm.box = None

        # check periodicity
        com_amb_prm = self._check_periodicity(com_amb_prm, 'complex')

        self.fixparm2amber(com_amb_prm)

        logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Complex...")
        if com_top_parm == 'amber' and self.INPUT['general']['PBRadii'] == 7:
            GMXMMPBSA_ERROR(f"The PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} is not compatible with Amber/OPLS "
                            f"topologies...")
        action = ChRad(com_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
        logging.info('Writing Normal Complex AMBER topology...')
        com_amb_prm.write_parm(self.complex_pmrtop)

        rec_indexes_string = ','.join(self.resi['REC']['string'])

        rec_hastop = True
        if self.FILES.receptor_top:
            logging.info('A Receptor topology file was defined. Using MT approach...')
            logging.info('Building AMBER Receptor Topology from GROMACS Receptor Topology...')
            rec_top = self.cleantop(self.FILES.receptor_top, self.indexes['REC'], 'receptor')

            if error_info := eq_strs(rec_top, self.receptor_str):
                if error_info[0] == 'atoms':
                    GMXMMPBSA_ERROR(f"The number of atoms in the topology ({error_info[1]}) and the receptor "
                                    f"structure ({error_info[2]}) are different. Please check this files and verify "
                                    f"that they are correct. Otherwise report the error...")
                else:
                    GMXMMPBSA_ERROR(f"The number of residues in the topology ({error_info[1]}) and the receptor "
                                    f"structure ({error_info[2]}) are different. Please check this files and verify "
                                    f"that they are correct. Otherwise report the error...")

            rec_top.coordinates = self.receptor_str.coordinates
            # rec_top.save(f"{self.FILES.prefix}REC.inpcrd", format='rst7', overwrite=True)
            if rec_top.impropers or rec_top.urey_bradleys:
                if com_top_parm == 'amber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is Amber/OPLS type while the '
                                    'Receptor is CHAMBER type!')
                rec_amb_prm = parmed.amber.ChamberParm.from_structure(rec_top)
            else:
                if com_top_parm == 'chamber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is CHAMBER type while the '
                                    'Receptor is Amber/OPLS type!')
                rec_amb_prm = parmed.amber.AmberParm.from_structure(rec_top)
            logging.info('Changing the Receptor residues name format from GROMACS to AMBER...')

            # check periodicity
            rec_amb_prm = self._check_periodicity(rec_amb_prm, 'receptor')

            self.fixparm2amber(rec_amb_prm)
        else:
            logging.info('No Receptor topology file was defined. Using ST approach...')
            logging.info('Building AMBER Receptor topology from Complex...')
            # we make a copy for receptor topology
            rec_amb_prm = self.molstr(com_amb_prm)
            rec_amb_prm.strip(f'!:{rec_indexes_string}')
            rec_hastop = False

        logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Receptor...")
        action = ChRad(rec_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
        logging.info('Writing Normal Receptor AMBER topology...')
        rec_amb_prm.write_parm(self.receptor_pmrtop)
        rec_amb_prm.save(f"{self.FILES.prefix}REC.inpcrd", format='rst7', overwrite=True)

        lig_hastop = True
        if self.FILES.ligand_top:
            logging.info('A Ligand Topology file was defined. Using MT approach...')
            logging.info('Building AMBER Ligand Topology from GROMACS Ligand Topology...')
            lig_top = self.cleantop(self.FILES.ligand_top, self.indexes['LIG'], 'ligand')

            if error_info := eq_strs(lig_top, self.ligand_str):
                if error_info[0] == 'atoms':
                    GMXMMPBSA_ERROR(f"The number of atoms in the topology ({error_info[1]}) and the ligand "
                                    f"structure ({error_info[2]}) are different. Please check this files and verify "
                                    f"that they are correct. Otherwise report the error...")
                else:
                    GMXMMPBSA_ERROR(f"The number of residues in the topology ({error_info[1]}) and the ligand "
                                    f"structure ({error_info[2]}) are different. Please check this files and verify "
                                    f"that they are correct. Otherwise report the error...")

            lig_top.coordinates = self.ligand_str.coordinates
            # lig_top.save(f"{self.FILES.prefix}LIG.inpcrd", format='rst7', overwrite=True)
            if lig_top.impropers or lig_top.urey_bradleys:
                if com_top_parm == 'amber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is Amber/OPLS type while the '
                                    'Ligand is CHAMBER type!')
                lig_amb_prm = parmed.amber.ChamberParm.from_structure(lig_top)
            else:
                if com_top_parm == 'chamber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is CHAMBER type while the '
                                    'Ligand is Amber/OPLS type!')
                lig_amb_prm = parmed.amber.AmberParm.from_structure(lig_top)
            logging.info('Changing the Ligand residues name format from GROMACS to AMBER...')

            # check periodicity
            lig_amb_prm = self._check_periodicity(lig_amb_prm, 'ligand')

            self.fixparm2amber(lig_amb_prm)
        else:
            logging.info('No Ligand topology file was defined. Using ST approach...')
            logging.info('Building AMBER Ligand topology from Complex...')
            # we make a copy for ligand topology
            lig_amb_prm = self.molstr(com_amb_prm)
            lig_amb_prm.strip(f':{rec_indexes_string}')
            lig_hastop = False
        logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Ligand...")
        action = ChRad(lig_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
        logging.info('Writing Normal Ligand AMBER topology...')
        lig_amb_prm.write_parm(self.ligand_pmrtop)
        lig_amb_prm.save(f"{self.FILES.prefix}LIG.inpcrd", format='rst7', overwrite=True)

        if self.INPUT['ala']['alarun']:
            logging.info('Building Mutant Complex Topology...')
            # get mutation index in complex
            self.com_mut_index, self.part_mut, self.part_index = self.getMutationInfo()
            mut_com_amb_prm = self.makeMutTop(com_amb_prm, self.com_mut_index)
            logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Mutant Complex...")
            action = ChRad(mut_com_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
            logging.info('Writing Mutant Complex AMBER topology...')
            mut_com_amb_prm.write_parm(self.mutant_complex_pmrtop)

            if self.part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor topology...')
                out_prmtop = self.mutant_receptor_pmrtop
                self.mutant_ligand_pmrtop = None
                if rec_hastop:
                    mtop = self.makeMutTop(rec_amb_prm, self.part_index)
                else:
                    mut_com_amb_prm.strip(f'!:{rec_indexes_string}')
                    mtop = mut_com_amb_prm
            else:
                logging.info('Detecting mutation in Ligand. Building Mutant Ligand topology...')
                out_prmtop = self.mutant_ligand_pmrtop
                self.mutant_receptor_pmrtop = None
                if lig_hastop:
                    mtop = self.makeMutTop(lig_amb_prm, self.part_index)
                else:
                    mut_com_amb_prm.strip(f':{rec_indexes_string}')
                    mtop = mut_com_amb_prm

            if com_top_parm == 'chamber':
                mut_prot_amb_prm = parmed.amber.ChamberParm.from_structure(mtop)
            else:
                mut_prot_amb_prm = parmed.amber.AmberParm.from_structure(mtop)
            logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Mutant "
                         f"{'Receptor' if self.part_mut == 'REC' else 'Ligand'}...")
            action = ChRad(mut_prot_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
            logging.info(f"Writing Mutant {'Receptor' if self.part_mut == 'REC' else 'Ligand'} AMBER topology...")
            mut_prot_amb_prm.write_parm(out_prmtop)
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)

    def _split_str(self, start, r, c, basename, struct, mut_index=0):
        end = start + (r[1] - r[0])
        mask = f'!:{start}-{end}'
        str_ = self.molstr(struct)
        if mut_index:
            str_ = self.makeMutTop(str_, mut_index, True)
        str_.strip(mask)
        str_file = f'{self.FILES.prefix}{basename}_F{c}.pdb'
        str_.save(str_file, 'pdb', True, renumber=False)
        return end, str_file


    @staticmethod
    def cleantop(top_file, ndx, id='complex'):
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
                GMXMMPBSA_ERROR(f'The atom {i} in the {id} index is not found in the topology file. Please check that '
                                'the files are consistent.')

        ranges = list2range(res_list)
        top.strip(f"!:{','.join(ranges['string'])}")

        return top

    def get_masks(self):
        rec_mask = ':' + ','.join(self.resi['REC']['string'])
        lig_mask = ':' + ','.join(self.resi['LIG']['string'])

        if self.INPUT['ala']['alarun']:
            self.resl[self.com_mut_index].set_mut(self.INPUT['ala']['mutant'])
        return rec_mask, lig_mask, self.resl

    def get_selected_residues(self, select, qm_sele=False):
        """
        Convert string selection format to amber index list
        """
        if qm_sele:
            com_top = parmed.load_file(self.complex_pmrtop)

        dist, res_selection = selector(select)
        residues_selection = {'rec': [], 'lig': []}
        rec_charge = 0
        lig_charge = 0
        if dist:
            for rres in self.resl:
                if rres.is_ligand():
                    continue
                for lres in self.resl:
                    if lres.is_receptor():
                        continue
                    for rat in self.complex_str.residues[rres - 1].atoms:
                        rat_coor = [rat.xx, rat.xy, rat.xz]
                        for lat in self.complex_str.residues[lres - 1].atoms:
                            lat_coor = [lat.xx, lat.xy, lat.xz]
                            if get_dist(rat_coor, lat_coor) <= dist:
                                if rres not in residues_selection['rec']:
                                    residues_selection['rec'].append(rres)
                                    if qm_sele:
                                        rec_charge += round(
                                            sum(atm.charge for atm in com_top.residues[rres - 1].atoms), 0)
                                if lres not in residues_selection['lig']:
                                    residues_selection['lig'].append(lres)
                                    if qm_sele:
                                        lig_charge += round(
                                            sum(atm.charge for atm in com_top.residues[lres - 1].atoms), 0)
                                break
        elif res_selection:
            for i in self.resl:
                rres = self.complex_str.residues[i - 1]
                if [rres.chain, rres.number, rres.insertion_code] in res_selection:
                    if i.is_ligand():
                        residues_selection['lig'].append(i)
                        if qm_sele:
                            rec_charge += round(sum(atm.charge for atm in com_top.residues[i - 1].atoms), 0)
                    else:
                        residues_selection['rec'].append(i)
                        if qm_sele:
                            lig_charge += round(sum(atm.charge for atm in com_top.residues[i - 1].atoms), 0)
                    res_selection.remove([rres.chain, rres.number, rres.insertion_code])
            for res in res_selection:
                logging.warning("We couldn't find this residue CHAIN:{} RES_NUM:{} ICODE: {}".format(*res))
            # check if residues in receptor and ligand was defined
            if not residues_selection['rec'] or not residues_selection['lig']:
                if not self.INPUT['ala']['alarun']:
                    GMXMMPBSA_ERROR('For decomposition analysis, you most define residues for both receptor and ligand!')
        else:
            for i in self.resl:
                if i.is_ligand():
                    residues_selection['lig'].append(i)
                    if qm_sele:
                        rec_charge += round(sum(atm.charge for atm in com_top.residues[i - 1].atoms), 0)
                else:
                    residues_selection['rec'].append(i)
                    if qm_sele:
                        lig_charge += round(sum(atm.charge for atm in com_top.residues[i - 1].atoms), 0)
        sele_res = sorted([r for m in residues_selection.values() for r in m], key=lambda x: x.index)
        return (sele_res, (rec_charge, lig_charge)) if qm_sele else sele_res

    def fixparm2amber(self, structure, str_name=None):

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

            # FIXME: remove Gromacs 4.x support
            # GROMACS 4.x save the pdb without atom element column, so parmed does not recognize some H atoms.
            # Parmed assigns 0 to the atomic number of these atoms. In order to correctly eliminate hydrogens,
            # it is necessary to assign the atomic number.
            if len(self.make_ndx) == 2:
                for atom in residue.atoms:
                    if 'H' in atom.name and atom.atomic_number == 0:
                        atom.atomic_number = 1
                # Remove H atoms. Only when using the pdb files with tleap to build the topologies
        if str_name:
            structure.strip('@/H')

    def getMutationInfo(self):
        if not self.INPUT['ala']['mutant_res']:
            GMXMMPBSA_ERROR("No residue for mutation was defined")
        # dict = { resind: [chain, resnum, icode]
        sele_res_dict = self.get_selected_residues(self.INPUT['ala']['mutant_res'])
        if len(sele_res_dict) != 1:
            GMXMMPBSA_ERROR('Only ONE mutant residue is allowed.')
        r = sele_res_dict[0]
        res = self.complex_str.residues[r - 1]
        icode = f':{res.insertion_code}' if res.insertion_code else ''
        if (not parmed.residue.AminoAcidResidue.has(res.name) or res.name in ['CYX', 'PRO', 'GLY'] or
                res.name == 'ALA' and self.INPUT['ala']['mutant'] == 'ALA'):
            GMXMMPBSA_ERROR(f"Selecting residue {res.chain}:{res.name}:{res.number}{icode} can't be mutated. Please, "
                            f"define a valid residue...")

        if r.is_receptor():
            part_index = r.id_index - 1
            part_mut = 'REC'
        elif r.is_ligand():
            part_index = r.id_index - 1
            part_mut = 'LIG'
        else:
            part_index = None
            part_mut = None
            if icode:
                GMXMMPBSA_ERROR(f'Residue {res.chain}:{res.number}:{res.insertion_code} not found')
            else:
                GMXMMPBSA_ERROR(f'Residue {res.chain}:{res.number} not found')

        # return r - 1 since r is the complex mutant index from amber selection format. Needed for top mutation only
        return r - 1, part_mut, part_index

    def makeMutTop(self, wt_top, mut_index, pdb=False):
        """

        :param wt_top: Amber parm from GROMACS topology
        :param mut_index: index of mutation in structure
        :return: Mutant AmberParm
        """
        mut_top = self.molstr(wt_top)
        mut_aa = self.INPUT['ala']['mutant']

        bb_atoms = 'N,H,CA,HA,C,O,HN'
        nterm_atoms = 'H1,H2,H3'
        cterm_atoms = 'OXT'
        sc_cb_atom = 'CB'
        sc_ala_atoms = ('HB,' +  # VAL, ILE, THR
                        'HB1,HB2,' + # charmm -> HB1, HB2
                        'HB3,' + # amber -> HB2, HB3
                        'CG1,CG2,OG1,' +  # VAL, ILE, THR
                        'OG,' +  # SER
                        'SG,' +  # CYS
                        'CG')

        if mut_aa in ['GLY', 'G']:
            # FIXME: allow terminal residues to mutate?
            strip_mask = f":{mut_index + 1} &!@{','.join([bb_atoms, nterm_atoms, cterm_atoms])}"
            if not pdb:
                strip_mask += f",{sc_cb_atom}"
        else:
            # FIXME: allow terminal residues to mutate?
            strip_mask = f":{mut_index + 1} &!@{','.join([bb_atoms, sc_cb_atom, nterm_atoms, cterm_atoms])}"
            if not pdb:
                strip_mask += f",{sc_ala_atoms}"
        mut_top.strip(strip_mask)

        h_atoms_prop = {}
        # get an example HB atom if not PDB
        if not pdb:
            for res in mut_top.residues:
                if res.name == mut_aa:
                    for at in res.atoms:
                        if (
                                mut_aa == 'GLY'
                                and at.name in ['HA2']
                                or mut_aa != 'GLY'
                                and at.name in ['HB2']
                        ):
                            h_atoms_prop['mass'] = at.mass
                            h_atoms_prop['element'] = at.element
                            h_atoms_prop['atomic_number'] = at.atomic_number
                            h_atoms_prop['charge'] = at.charge
                            h_atoms_prop['atom_type'] = at.atom_type
                            h_atoms_prop['type'] = at.type
                            break
                    break
        cb_atom = None
        ca_atom = None
        logging.info(
            f"Mutating {self.complex_str.residues[mut_index].chain}/{self.complex_str.residues[mut_index].number} "
            f"{self.complex_str.residues[mut_index].name} by {mut_aa}")

        mutant_resname = mut_top.residues[mut_index].name

        mut_top.residues[mut_index].name = mut_aa

        for at in mut_top.residues[mut_index].atoms:
            if mut_aa == 'GlY':
                if at.name == 'CA':
                    ca_atom = at
                if at.name in ['CB']:
                    at.name = 'HA2'
                    ca_atom.xx, ca_atom.xy, ca_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [ca_atom.xx, ca_atom.xy,
                         ca_atom.xz, at.xx, at.xy,
                         at.xz], 1.09)
                    for ref_at in h_atoms_prop:
                        setattr(at, ref_at, h_atoms_prop[ref_at])
                elif at.name in ['HA']:
                    at.name = 'HA1'
                    at.type = h_atoms_prop['type']
                    at.atom_type = h_atoms_prop['atom_type']
            else:
                if at.name == 'CB':
                    cb_atom = at
                    continue
                if at.name == 'CG2':  # VAL, LEU and THR
                    at.name = 'HB2'
                    cb_atom.xx, cb_atom.xy, cb_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [cb_atom.xx, cb_atom.xy,
                         cb_atom.xz, at.xx, at.xy,
                         at.xz], 1.09)
                    for ref_at in h_atoms_prop:
                        setattr(at, ref_at, h_atoms_prop[ref_at])
                elif at.name in ['HB']:  # VAL, LEU and THR
                    at.name = 'HB1'
                    at.type = h_atoms_prop['type']
                    at.atom_type = h_atoms_prop['atom_type']
                elif at.name in ['CG', 'OG', 'SG',  # LEU, PHE, TRP, MET, TYR, ARG, LYS, ASN, GLN, ASP, GLU, HIS,
                                 # PRO (EXCLUDED), CYS (EXCLUDED IF S-S), SER
                                 'CG1', 'OG1'  # VAL, LEU and THR
                                 ]:
                    at.name = 'HB3'
                    cb_atom.xx, cb_atom.xy, cb_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [cb_atom.xx, cb_atom.xy,
                         cb_atom.xz, at.xx, at.xy,
                         at.xz], 1.09)
                    for ref_at in h_atoms_prop:
                        setattr(at, ref_at, h_atoms_prop[ref_at])
                elif at.name in ['HB1', 'HB2', 'HB3']:
                    at.type = h_atoms_prop['type']
                    at.atom_type = h_atoms_prop['atom_type']
        # change intdiel if cas_intdiel was defined before end the mutation process
        if self.INPUT['ala']['cas_intdiel']:
            if self.INPUT['gb']['gbrun']:
                if self.INPUT['gb']['intdiel'] != 1.0:
                    logging.warning('Both cas_intdiel and intdiel were defined. The dielectric constants associated '
                                    'with cas_intdiel will be ignored and intdiel will be used instead')
                elif mutant_resname in polar_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_polar']
                    logging.info(f"Setting intdiel = intdiel_polar = {self.INPUT['ala']['intdiel_polar']} for "
                                 f"Alanine scanning")
                elif mutant_resname in nonpolar_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_nonpolar']
                    logging.info(f"Setting intdiel = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                                 f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(f"Setting intdiel = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for Alanine "
                                 f"scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(f"Setting intdiel = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
                                 f"scanning")
                else:
                    logging.warning(f"Unclassified mutant residue {mutant_resname}. The default "
                                    f"intdiel will be used")
            if self.INPUT['gbnsr6']['gbnsr6run']:
                if self.INPUT['gbnsr6']['epsin'] != 1.0:
                    logging.warning('Both cas_intdiel and epsin were defined. The dielectric constants associated '
                                    'with cas_intdiel will be ignored and epsin will be used instead')
                elif mutant_resname in polar_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_polar']
                    logging.info(f"Setting epsin = intdiel_polar = {self.INPUT['ala']['intdiel_polar']} for "
                                 f"Alanine scanning")
                elif mutant_resname in nonpolar_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_nonpolar']
                    logging.info(f"Setting epsin = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                                 f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(f"Setting epsin = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for Alanine "
                                 f"scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(f"Setting epsin = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
                                 f"scanning")
                else:
                    logging.warning(f"Unclassified mutant residue {mutant_resname}. The default "
                                    f"intdiel will be used")

            if self.INPUT['pb']['pbrun']:
                if self.INPUT['pb']['indi'] != 1.0:
                    logging.warning('Both cas_intdiel and indi were defined. The dielectric constants associated with '
                                    'cas_intdiel will be ignored and indi will be used instead')
                elif mutant_resname in polar_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_polar']
                    logging.info(f"Setting indi = intdiel_polar = {self.INPUT['ala']['intdiel_polar']} for Alanine scanning")
                elif mutant_resname in nonpolar_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_nonpolar']
                    logging.info(f"Setting indi = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                                 f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(f"Setting intdiel = indi = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for "
                                 f"Alanine scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(f"Setting indi = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
                                 f"scanning")
                else:
                    logging.warning(f"Unclassified mutant residue {mutant_resname}. The default indi will be used")
        return mut_top

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
            trjconv_args = self.trjconv + ['-f', self.FILES.complex_trajs[i], '-s', self.FILES.complex_tpr, '-o',
                                           f'COM_traj_{i}.xtc', '-n', self.FILES.complex_index]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(trjconv_args))
            c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c6)
            if c6.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.complex_trajs[i]))
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
                                               '-o', 'REC_traj_{}.xtc'.format(i), '-n', self.FILES.receptor_index]
                logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                              (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                              '| ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_subprocess_output(c6)
                if c6.wait():  # if it quits with return code != 0
                    GMXMMPBSA_ERROR(
                        '%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.receptor_trajs[i]))
                new_trajs.append('REC_traj_{}.xtc'.format(i))
            self.FILES.receptor_trajs = new_trajs

        if self.FILES.ligand_tpr:
            logging.info('Cleaning normal ligand trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.ligand_trajs)):
                trjconv_echo_args = echo_command + ['GMXMMPBSA_LIG']
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = self.trjconv + ['-f', self.FILES.ligand_trajs[i], '-s', self.FILES.ligand_tpr, '-o',
                                               'LIG_traj_{}.xtc'.format(i), '-n', self.FILES.ligand_index]
                logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                              (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                              '| ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_subprocess_output(c6)
                if c6.wait():  # if it quits with return code != 0
                    GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.ligand_trajs[i]))
                new_trajs.append('LIG_traj_{}.xtc'.format(i))
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
                GMXMMPBSA_ERROR(f'The number of residues of the complex ({len(com_str.residues)}) and of the '
                                f'reference structure ({len(ref_str.residues)}) are different. Please check that the '
                                f'reference structure is correct')
            for c, res in enumerate(ref_str.residues):
                if com_str.residues[c].number != res.number or com_str.residues[c].name != res.name:
                    GMXMMPBSA_ERROR('There is no match between the complex and the reference structure used. An '
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
        com_str.save(f'{self.FILES.prefix}COM_FIXED.pdb', 'pdb', True, renumber=False)
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

    @staticmethod
    def molstr(data):
        if type(data) == str:
            # data is a pdb file
            pdb_file = data
            try:
                new_str = []
                with open(pdb_file) as fo:
                    fo = fo.readlines()
                    for line in fo:
                        if 'MODEL' in line or 'ENDMDL' in line:
                            continue
                        # check new charmm-gui format for Amber ff19SB (with N- and C- terminals)
                        if 'ATOM' in line:
                            resn = line[17:21].strip()
                            if len(resn) == 4 and resn.startswith(('N', 'C')):
                                line = f'{line[:17]}{resn[1:]} {line[21:]}'
                        new_str.append(line)
                with open(pdb_file, 'w') as fw:
                    for x in new_str:
                        fw.write(x)
            except IOError as e:
                GMXMMPBSA_ERROR('', str(e))

            structure = parmed.read_PDB(pdb_file)
        else:
            # data is Structure, AmberParm, ChamberParm or GromacsTopologyFile. This make a copy
            structure = data.__copy__()
            for c, at in enumerate(structure.atoms, start=1):
                at.number = c
        return structure