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

import string
import textwrap
from pathlib import Path

import parmed
from xBFreE.exceptions import *
from xBFreE.mmpbsa.alamdcrd import _scaledistance
from xBFreE.utils.molecule import get_dist, selector, list2range

positive_aa = ['LYS', 'ARG', 'HIP']
negative_aa = ['GLU', 'ASP']
nonpolar_aa = ['PHE', 'TRP', 'VAL', 'ILE', 'LEU', 'MET', 'PRO', 'CYX', 'ALA', 'GLY']
polar_aa = ['TYR', 'SER', 'THR', 'CYM', 'CYS', 'HIE', 'HID', 'GLN', 'ASN', 'ASH', 'GLH', 'LYN']
chains_letters = list(string.ascii_uppercase)


class BuildTop:
    def __init__(self, FILES, INPUT, external_programs):
        self.FILES = FILES
        self.INPUT = INPUT
        self.external_progs = external_programs

        # create the * prmtop variables for compatibility with the original code
        self.complex_prmtop = 'COM.prmtop'
        self.receptor_prmtop = 'REC.prmtop'
        self.ligand_prmtop = 'LIG.prmtop'

        self.ref_str = None
        self.complex_str = None
        self.receptor_str = None
        self.ligand_str = None

        self.mutant_complex_prmtop = 'MUT_COM.prmtop'
        self.mutant_receptor_prmtop = 'MUT_REC.prmtop'
        self.mutant_ligand_prmtop = 'MUT_LIG.prmtop'

        self.complex_str_file = 'COM.pdb'
        self.receptor_str_file = 'REC.pdb'
        self.ligand_str_file = 'LIG.pdb'


        # FIXME: adapt this function to implemnted md programs
        # self.checkFiles()

    def checkFiles(self):
        pass

    def buildTopology(self):
        pass

    def decomp_qmmm_ressel(self):
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
                            xBFreEErrorLogging('The maximum number of iterations to select interaction residues was '
                                               'reached. Please set print_res with a valid selection.')
                        it += 1

                    logging.info(f"Selecting residues by distance ({round(cutoff, 1)} Å) between "
                                 f"receptor and ligand for decomposition analysis...")
                else:
                    logging.info(
                        f"Selecting residues by distance ({self.INPUT['decomp']['print_res'].split()[1]} Å) between "
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
                            xBFreEErrorLogging('The maximum number of iterations to select interaction residues was '
                                               'reached. Please set print_res with a valid selection.')
                        it += 1

                    logging.info(f"Selecting residues by distance ({round(cutoff, 1)} Å) between "
                                 f"receptor and ligand for QM/MM calculation...")
                else:
                    logging.info(
                        f"Selecting residues by distance ({self.INPUT['gb']['qm_residues'].split()[1]} Å) between "
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

    def get_selected_residues(self, select, qm_sele=False):
        """
        Convert string selection format to amber index list
        """
        if qm_sele:
            com_top = parmed.load_file(self.complex_prmtop)

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
                    xBFreEErrorLogging(
                        'For decomposition analysis, you most define residues for both receptor and ligand!')
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

    def get_masks(self):
        rec_mask = ':' + ','.join(self.resi['REC']['string'])
        lig_mask = ':' + ','.join(self.resi['LIG']['string'])

        if self.INPUT['ala']['alarun']:
            self.resl[self.com_mut_index].set_mut(self.INPUT['ala']['mutant'])
        return rec_mask, lig_mask, self.resl

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
                        'HB1,HB2,' +  # charmm -> HB1, HB2
                        'HB3,' +  # amber -> HB2, HB3
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
                    logging.info(
                        f"Setting intdiel = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                        f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(
                        f"Setting intdiel = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for Alanine "
                        f"scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(
                        f"Setting intdiel = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
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
                    logging.info(
                        f"Setting epsin = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                        f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(
                        f"Setting epsin = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for Alanine "
                        f"scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(
                        f"Setting epsin = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
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
                    logging.info(
                        f"Setting indi = intdiel_polar = {self.INPUT['ala']['intdiel_polar']} for Alanine scanning")
                elif mutant_resname in nonpolar_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_nonpolar']
                    logging.info(
                        f"Setting indi = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                        f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(
                        f"Setting intdiel = indi = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for "
                        f"Alanine scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(
                        f"Setting indi = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
                        f"scanning")
                else:
                    logging.warning(f"Unclassified mutant residue {mutant_resname}. The default indi will be used")
        return mut_top

    def getMutationInfo(self):
        if not self.INPUT['ala']['mutant_res']:
            xBFreEErrorLogging("No residue for mutation was defined")
        # dict = { resind: [chain, resnum, icode]
        sele_res_dict = self.get_selected_residues(self.INPUT['ala']['mutant_res'])
        if len(sele_res_dict) != 1:
            xBFreEErrorLogging('Only ONE mutant residue is allowed.')
        r = sele_res_dict[0]
        res = self.complex_str.residues[r - 1]
        icode = f':{res.insertion_code}' if res.insertion_code else ''
        if (not parmed.residue.AminoAcidResidue.has(res.name) or res.name in ['CYX', 'PRO', 'GLY'] or
                res.name == 'ALA' and self.INPUT['ala']['mutant'] == 'ALA'):
            xBFreEErrorLogging(
                f"Selecting residue {res.chain}:{res.name}:{res.number}{icode} can't be mutated. Please, "
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
                xBFreEErrorLogging(f'Residue {res.chain}:{res.number}:{res.insertion_code} not found')
            else:
                xBFreEErrorLogging(f'Residue {res.chain}:{res.number} not found')

        # return r - 1 since r is the complex mutant index from amber selection format. Needed for top mutation only
        return r - 1, part_mut, part_index

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
                xBFreEErrorLogging('', str(e))

            structure = parmed.read_PDB(pdb_file)
        else:
            # data is Structure, AmberParm, ChamberParm or GromacsTopologyFile. This make a copy
            structure = data.__copy__()
            for c, at in enumerate(structure.atoms, start=1):
                at.number = c
        return structure
