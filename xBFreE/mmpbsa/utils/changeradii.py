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


import json
import re
from pathlib import Path

import parmed
from parmed.residue import AminoAcidResidue, RNAResidue, DNAResidue
from parmed.tools.changeradii import _screen1
from xBFreE.mmpbsa import data as radii_path

"""
bondi           -> | atom_number | atom_type |                    | mass | radii |  (General)
amber6          -> | atom_number | atom_type | bonded_atom_number | mass | radii |  (General)
mbondi          -> | atom_number | atom_type | bonded_atom_number | mass | radii |  (General)
mbondi2         -> | atom_number | atom_type | bonded_atom_number | mass | radii |  (General)
mbondi3         -> | atom_number | atom_type | bonded_atom_number | mass | res_name | atom_name | radii |   (General)
mbondi_pb2      -> | atom_number | atom_type | bonded_atom_number | mass | radii |  (General) 
mbondi_pb3      -> | atom_number | atom_type | bonded_atom_number | mass | radii |  (General)
tyl06           -> | res_name    | atom_name |                    |      | radii |  (Protein, NA, Ions and Water) 
yamagishi       -> | res_name    | atom_name |                    |      | radii |  (Protein)
"""


class LoadRadii:
    def __init__(self, radiinames, searchpath=None):

        self.radiinames = radiinames
        self.searchpath = searchpath
        self.radiinames_list = self.radiinames.split('+')
        self.radius_set_text = ''
        self._check_radii_names()
        self.mapper = self._get_mapper()

    def reset_radii(self, parm):
        """
        define each atom radii as None to check if was assigned a new radii or not
        """
        for at in parm.atoms:
            at.solvent_radius = None

    def check_radii(self, parm):

        atm_wo_radii = []
        for r in parm.residues:
            # if r.name in ['SOD', 'CLA', 'TIP3']:
            #     continue
            atm_wo_radii.extend(str(at) for at in r.atoms if at.solvent_radius is None)
        if atm_wo_radii:
            atm_wo_radii_text = '\n'.join(atm_wo_radii)
            raise TypeError("There are atoms that could not be assigned any radii value. It may be due to the radii "
                            "you are using or an internal error. Please identify which of the two options it is and "
                            f"report it if necessary.\n The atoms are listed below:\n{atm_wo_radii_text}")

    def _check_radii_names(self):
        """
        Check some problems in the list of radii names defines by de user
        :return:
        """
        # check duplicated
        for r in self.radiinames_list:
            if self.radiinames_list.count(r) > 1:
                raise TypeError(f'Radii {r} is duplicated')

        # check for compatibility
        if 'charmm_radii' in self.radiinames_list and len(self.radiinames_list) > 1:
            raise TypeError('charmm_radii must be defined alone')

    def _get_mapper(self):
        """
        Read the radii specified in a dictionary
        :return:
        """

        # check if radii name exist
        sp = Path(radii_path.__file__).parent.joinpath('radii')

        radii_file_list = []
        custom_path = True
        default_path = False
        for rn in self.radiinames_list:
            if self.searchpath:
                rp = Path(self.searchpath).joinpath(f'{rn}.json')
                if rp.exists():
                    radii_file_list.append(rp)
                    continue
                else:
                    custom_path = False
            rp = sp.joinpath(f'{rn}.json')
            if rp.exists():
                radii_file_list.append(rp)
                default_path = True

        if not custom_path and not default_path:
            raise TypeError('Not exist!')

        mapper = {}
        for radiifile in radii_file_list:
            radiiname = radiifile.stem
            with radiifile.open() as openfile:
                data = json.load(openfile)
                mapper[radiiname] = self._mapper_function(data)

        rbr = []
        atbr = []
        rs_list = []
        for name, value in mapper.items():
            if 'atom_number' in value['type'] and 'residue' in value['type']:
                atbr.append(name)
                rbr.append(name)
            elif 'atom_number' in value['type']:
                atbr.append(name)
            else:
                rbr.append(name)
            rs_list.append(value['radius_set'])

        self.radius_set_text = ' + '.join(rs_list)
        if len(rbr) > 1:
            raise TypeError(f"You defined {'+'.join(rbr)} radii names, which are residue-based. Only one is allowed...")
        if len(atbr) > 1:
            raise TypeError(
                f"You defined {'+'.join(atbr)} radii names, which are residue-based. Only one is allowed...")
        return mapper

    def _mapper_function(self, data):
        radii_mapper = {'data': {}, 'type': data['type']}
        keyword = 'atom_number' if 'atom_number' in data['radii'][0] else 'res_name'
        for a in data['radii']:
            if not radii_mapper['data'].get(a[keyword]):
                radii_mapper['data'][a[keyword]] = [a]
            else:
                radii_mapper['data'][a[keyword]].append(a)

        radii_mapper['radius_set'] = data.get('RADIUS_SET')

        return radii_mapper

    def assign_radii(self, parm):
        # check top and radii compatibility
        non_bmr = None
        if parm.impropers or parm.urey_bradleys:
            # check if charmm topology contains only biomolecular residues
            non_bmr = sum(
                1
                for res in parm.residues
                if not (AminoAcidResidue.has(res.name)
                        or DNAResidue.has(res.name)
                        or RNAResidue.has(res.name)
                        or res.name in ['TIP3P', 'TIP3', 'TP3', 'TIPS3P', 'TIP3o', 'TIP4P', 'TIP4PEW', 'T4E', 'TIP4PD',
                                        'TIP5P',
                                        'SPC', 'SPC/E', 'SPCE', 'WAT', 'OPC', 'HOH'])
            )

        if isinstance(parm, parmed.amber.AmberParm) and 'charmm_radii' in self.radiinames_list:
            raise

        self.reset_radii(parm)

        for r in parm.residues:
            if r.name in ['SOD', 'CLA', 'TIP3']:
                continue
            for radiiname, data in self.mapper.items():
                radiitype, radiidata = data['type'], data['data']
                if radiitype == 'residue' and r.name in radiidata:
                    # check if is charmm terminals to get radii form amber terminal
                    if len(r.name) == 3:
                        if sum(atm.name.startswith('HT') for atm in r.atoms) == 3:
                            temp_resname = f'N{r.name}'
                        elif sum(atm.name.startswith('H') for atm in r.atoms) == 3:
                            temp_resname = f'N{r.name}'
                        elif sum(atm.name.startswith('OT') for atm in r.atoms) == 2:
                            temp_resname = f'C{r.name}'
                        elif sum(atm.name.startswith('OXT') for atm in r.atoms) == 1:
                            temp_resname = f'C{r.name}'
                        else:
                            temp_resname = r.name
                    else:
                        temp_resname = r.name

                    for atm in r.atoms:
                        for atmr in radiidata[temp_resname]:
                            if atmr['atom_name'].startswith('^'):
                                match = re.match(atmr['atom_name'], atm.name)
                            else:
                                match = re.fullmatch(atmr['atom_name'], atm.name)
                            if match:
                                atm.solvent_radius = atmr['radii']
                                break
                    break
                elif radiitype == "atom_number+bonded_atom_number":
                    for atm in r.atoms:
                        if atm.atomic_number in radiidata:
                            for atmr in radiidata[atm.atomic_number]:
                                if atmr['atom_type'] and atmr['mass']:
                                    if re.match(atmr['atom_type'], atm.type) and self.eval_mass(atm.mass, atmr['mass']):
                                        atm.solvent_radius = atmr['radii']
                                elif atmr['bonded_atom_number']:
                                    bondeds = list(atm.bond_partners)
                                    if bondeds[0].atomic_number == atmr['bonded_atom_number']:
                                        atm.solvent_radius = atmr['radii']
                                else:
                                    atm.solvent_radius = atmr['radii']
                        else:
                            atm.solvent_radius = radiidata['*']['radii']
                elif radiitype == "atom_number+bonded_atom_number+residue":
                    cterm = any(atm.name in ['OXT', 'OT2'] for atm in r.atoms)
                    for atm in r.atoms:
                        if atm.atomic_number in radiidata:
                            for atmr in radiidata[atm.atomic_number]:
                                if atmr['atom_type'] and atmr['mass']:
                                    if re.match(atmr['atom_type'], atm.type) and self.eval_mass(atm.mass, atmr['mass']):
                                        atm.solvent_radius = atmr['radii']
                                elif atmr['bonded_atom_number']:
                                    bondeds = list(atm.bond_partners)
                                    if bondeds[0].atomic_number == atmr['bonded_atom_number']:
                                        atm.solvent_radius = atmr['radii']
                                elif atmr.get('res_name'):
                                    if atmr['atom_name'].startswith('^'):
                                        match = re.match(atmr['atom_name'], atm.name)
                                    elif atmr['atom_name'] == 'OTER' and cterm:
                                        match = re.fullmatch('O', atm.name)
                                    else:
                                        match = re.fullmatch(atmr['atom_name'], atm.name)
                                    if match:
                                        if atmr['res_name'] == '*':
                                            atm.solvent_radius = atmr['radii']
                                        elif atm.residue.name == atmr['res_name']:
                                            atm.solvent_radius = atmr['radii']
                                else:
                                    atm.solvent_radius = atmr['radii']
                        else:
                            atm.solvent_radius = radiidata['*']['radii']
                elif radiitype == "atom_number":
                    for atm in r.atoms:
                        if atm.atomic_number in radiidata:
                            for atmr in radiidata[atm.atomic_number]:
                                if atmr['atom_type'] and atmr['mass']:
                                    if re.match(atmr['atom_type'], atm.type) and self.eval_mass(atm.mass, atmr['mass']):
                                        atm.solvent_radius = atmr['radii']
                                else:
                                    atm.solvent_radius = atmr['radii']
                        else:
                            atm.solvent_radius = radiidata['*']
                elif radiitype == 'charmm_radii':
                    # FIXME: include scaling factor ???
                    for atm in r.atoms:
                        for atmr in radiidata['*']:
                            match = None
                            if atmr.get('atom_name'):
                                if atmr['atom_name'].startswith('^'):
                                    match = re.match(atmr['atom_name'], atm.name)
                                else:
                                    match = re.fullmatch(atmr['atom_name'], atm.name)
                            elif atmr.get('atom_type'):
                                if atmr['atom_type'].startswith('^'):
                                    match = re.match(atmr['atom_type'], atm.type)
                                else:
                                    match = re.fullmatch(atmr['atom_type'], atm.type)
                            if match:
                                atm.solvent_radius = atmr['radii']
                        if radiidata.get(r.name):
                            for atmr in radiidata[r.name]:
                                match = None
                                if atmr.get('atom_name'):
                                    if atmr['atom_name'].startswith('^'):
                                        if atmr.get('bonded'):
                                            bondeds = list(atm.bond_partners)
                                            match = re.match(atmr['atom_name'], atm.name) and atmr['bonded'] in bondeds
                                        else:
                                            match = re.match(atmr['atom_name'], atm.name)
                                    elif atmr.get('bonded'):
                                        bondeds = list(atm.bond_partners)
                                        match = re.fullmatch(atmr['atom_name'], atm.name) and atmr['bonded'] in bondeds
                                    else:
                                        match = re.fullmatch(atmr['atom_name'], atm.name)
                                elif atmr.get('atom_type'):
                                    if atmr['atom_type'].startswith('^'):
                                        match = re.match(atmr['atom_type'], atm.type)
                                    else:
                                        match = re.fullmatch(atmr['atom_type'], atm.type)
                                if match:
                                    atm.solvent_radius = atmr['radii']

        self.check_radii(parm)
        _screen1(parm)

    @staticmethod
    def eval_mass(mass, massexp):
        s = ''
        m = ''
        for l in massexp:
            if l in ['>', '>', '=']:
                s += l
            else:
                m += l
        if s == ">=":
            return mass >= float(m)
        elif s == ">":
            return mass > float(m)
        elif s == "<=":
            return mass <= float(m)
        elif s == "<":
            return mass < float(m)
        else:
            return mass == float(m)


#
#
#
# t1 = parmed.load_file(
#     '/home/mario/Documents/gmx_MMPBSA_tests/2.0.x/amber_example/step3_input.parm7')
#
# ChRad(t1, 'charmm_radii')
#
t2 = LoadRadii('tyl06')
# t = parmed.load_file('/home/mario/Documents/gmx_MMPBSA_tests/2.0.x/amber_example/step3_input.parm7')
t = parmed.load_file('/home/mario/Documents/gmx_MMPBSA_tests/2.0.x/Protein_protein/topol.top')
t2.assign_radii(t)

for at in t.residues[0].atoms:
    print(at, at.solvent_radius)

#
# # # t2.assign_radii()
# # for x, y in zip(t1.atoms, t2.top.atoms):
# #     if x.residue.name == 'TIP3':
# #         continue
# #     if x.solvent_radius != y.solvent_radius:
# #         print(x.residue.name, x.name, x.type, x.solvent_radius, y.name, y.type, y.solvent_radius, x.solvent_radius == \
# #               y.solvent_radius)
#
# print(t2.radius_set_text)
