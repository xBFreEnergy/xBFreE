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

import parmed
from xBFreE.exceptions import *
from xBFreE.mmpbsa.utils.molecule import check_str, get_dist
from xBFreE.mmpbsa.utils.misc import selector
from xBFreE.mmpbsa.alamdcrd import _scaledistance
import logging
import string

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
        self.complex_pmrtop = 'COM.prmtop'
        self.receptor_pmrtop = 'REC.prmtop'
        self.ligand_pmrtop = 'LIG.prmtop'

        self.ref_str = None
        self.complex_str = None
        self.receptor_str = None
        self.ligand_str = None

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

    def buildTopology(self):
        pass

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