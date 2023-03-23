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
from xBFreE.exceptions import GMXMMPBSA_ERROR
import logging
from math import sqrt
import re
from string import ascii_letters


def _get_restype(resname):
    if resname == 'LYN':
        return 'LYS'
    elif resname == 'ASH':
        return 'ASP'
    elif resname == 'GLH':
        return 'GLU'
    elif resname in ['HIP', 'HIE', 'HID']:
        return 'HIS'
    elif resname in ['CYX', 'CYM']:
        return 'CYS'
    else:
        return resname


def eq_strs(struct1, struct2, noh=False, molid='complex'):
    if len(struct1.atoms) != len(struct2.atoms):
        if not noh:
            return 'atoms', len(struct1.atoms), len(struct2.atoms)
        na1 = sum(not x.startswith('H') for x in struct1.atoms)
        na2 = sum(not x.startswith('H') for x in struct2.atoms)
        if na1 != na2:
            GMXMMPBSA_ERROR(f"The number of atoms in the topology ({len(struct1.atoms)}) and the {molid} structure "
                            f"({len(struct2.atoms)}) are different. Please check these files and verify that they are "
                            f"correct. Otherwise report the error...")

    elif len(struct1.residues) != len(struct2.residues):
        GMXMMPBSA_ERROR(f"The number of residues in the topology ({len(struct1.residues)}) and the {molid} structure "
                        f"({len(struct2.residues)}) are different. Please check these files and verify that they are "
                        f"correct. Otherwise report the error...")
    else:
        return


def check_str(structure, ref=False, skip=False):
    if isinstance(structure, str):
        refstr = parmed.read_PDB(structure)
    else:
        refstr = structure

    previous = 0
    ind = 1
    res_dict = {}
    duplicates = []
    for res in refstr.residues:
        if 'LP' in res.name:
            GMXMMPBSA_ERROR('The LP pseudo-atom is not supported. Please remove them following this instructions: '
                            'https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/Protein_ligand_LPH_atoms_CHARMMff/')
        if res.chain == '':
            if ref:
                GMXMMPBSA_ERROR('The reference structure used is inconsistent. The following residue does not have a '
                                f'chain ID: {res.number}:{res.name}')
            elif not previous:
                res_dict[ind] = [[res.number, res.name, res.insertion_code]]
            elif res.number - previous in [0, 1]:
                res_dict[ind].append([res.number, res.name, res.insertion_code])
            else:
                ind += 1
                res_dict[ind] = [[res.number, res.name, res.insertion_code]]
            previous = res.number
        elif res.chain not in res_dict:
            res_dict[res.chain] = [[res.number, res.name, res.insertion_code]]
        else:
            res_dict[res.chain].append([res.number, res.name, res.insertion_code])

    for chain, resl in res_dict.items():
        res_id_list = [[x, x2] for x, x1, x2 in resl]
        duplicates.extend(
            f'{chain}:{resl[c][0]}:{resl[c][1]}:{resl[c][2]}'
            for c, x in enumerate(res_id_list)
            if res_id_list.count(x) > 1
        )

    if ref:
        if duplicates:
            GMXMMPBSA_ERROR(f'The reference structure used is inconsistent. The following residues are duplicates:\n'
                            f' {", ".join(duplicates)}')
    elif skip:
        if duplicates:
            return refstr
    elif duplicates:
        logging.warning(f'The complex structure used is inconsistent. The following residues are duplicates:\n'
                        f' {", ".join(duplicates)}')
    return refstr


def selector(selection: str):
    string_list = re.split(r"\s|;\s*", selection)
    dist = None
    # exclude = None
    res_selections = []
    if selection == 'all':
        pass
    elif selection.startswith('within'):
        try:
            dist = float(string_list[1])
        except:
            GMXMMPBSA_ERROR(f'Invalid dist, we expected a float value but we get "{string_list[1]}"')
    else:
        # try to process residue selection
        for s in string_list:
            n = re.split(r":\s*|/\s*", s)
            if len(n) != 2 or n[0] not in ascii_letters:
                GMXMMPBSA_ERROR(f'We expected something like this: A/2-10,35,41 B/104 but we get {s} instead')
            chain = n[0]
            resl = n[1].split(',')
            for r in resl:
                rr = r.split('-')
                if len(rr) == 1:
                    ri = [chain, int(rr[0]), ''] if rr[0][-1] not in ascii_letters else [chain, int(rr[0][:-1]), rr[0][-1]]
                    if ri in res_selections:
                        logging.warning('Found duplicated residue in selection: CHAIN:{} RES_NUM:{} ICODE: '
                                        '{}'.format(*ri))
                        continue
                    res_selections.append(ri)
                else:
                    try:
                        start = int(rr[0])
                        end = int(rr[1]) + 1
                    except:
                        GMXMMPBSA_ERROR(f'When residues range is defined, start and end most be integer but we get'
                                        f' {rr[0]} and {rr[1]}')
                    for cr in range(start, end):
                        if [chain, cr, ''] in res_selections:
                            logging.warning('Found duplicated residue in selection: CHAIN:{} RES_NUM:{} ICODE: '
                                            '{}'.format(chain, cr, ''))
                            continue
                        res_selections.append([chain, cr, ''])
    return dist, res_selections

class Residue(object):
    def __init__(self, index, number, chain, mol_id, id_index, name, icode=''):
        self.index = int(index)
        self.number = number
        self.chain = chain
        self.mol_id = mol_id
        self.id_index = id_index
        self.name = name
        self.icode = icode
        self.mutant_label = None
        self.string = f"{mol_id}:{chain}:{name}:{number}:{icode}" if icode else f"{mol_id}:{chain}:{name}:{number}"
        self.mutant_string = None

    def __repr__(self):
        text = f"{type(self).__name__}(index: {self.index}, {self.mol_id}:{self.chain}:{self.name}:{self.number}"
        if self.icode:
            text += f":{self.icode}"
        text += ')'
        return text

    def __str__(self):
        return f"{self.index}"

    def __add__(self, other):
        if isinstance(other, Residue):
            return int(self.index + other.index)
        return int(self.index + other)

    def __sub__(self, other):
        if isinstance(other, Residue):
            return int(self.index - other.index)
        return int(self.index - other)

    def __int__(self):
        return self.index

    def is_mutant(self):
        return bool(self.mutant_label)

    def is_receptor(self):
        return self.mol_id == 'R'

    def is_ligand(self):
        return self.mol_id == 'L'

    def issame(self, other):
        pass

    def set_mut(self, mut):
        self.mutant_label = f'{self.chain}/{self.number}{f":{self.icode}" if self.icode else ""} - {self.name}x{mut}'
        self.mutant_string = (f"{self.mol_id}:{self.chain}:{mut}:{self.number}:{self.icode}" if self.icode
                              else f"{self.mol_id}:{self.chain}:{mut}:{self.number}")


def list2range(input_list):
    """
    Convert a list in list of ranges
    :return: list of ranges, string format of the list of ranges
    """

    def _add(temp):
        if len(temp) == 1:
            ranges_str.append(f"{temp[0]}")
            ranges.append([temp[0], temp[0]])
        else:
            ranges_str.append(f"{str(temp[0])}-{str(temp[-1])}")
            ranges.append([temp[0], temp[-1]])

    ranges = []
    ranges_str = []
    if not input_list:
        return ''
    temp = []
    previous = None

    ilist = sorted(input_list, key=lambda x: x.index if isinstance(x, Residue) else x)

    for x in ilist:
        if not previous:
            temp.append(x)
        elif x == previous + 1:
            temp.append(x)
        else:
            _add(temp)
            temp = [x]
        if x == ilist[-1]:
            _add(temp)
        previous = x
    return {'num': ranges, 'string': ranges_str}


def res2map(indexes, com_file):
    """
    :param com_str:
    :return:
    """
    res_list = []
    rec_list = []
    lig_list = []
    com_len = len(indexes['COM']['COM'])
    if isinstance(com_file, parmed.Structure):
        com_str = com_file
    else:
        com_str = parmed.load_file(com_file)

    resindex = 1
    rec_index = 1
    lig_index = 1
    proc_res = None
    for i in range(com_len):
        res = [com_str.atoms[i].residue.chain, com_str.atoms[i].residue.number, com_str.atoms[i].residue.name,
               com_str.atoms[i].residue.insertion_code]
        # We check who owns the residue corresponding to this atom
        if indexes['COM']['COM'][i] in indexes['COM']['REC']:
            # save residue number in the rec list
            if res != proc_res and resindex not in res_list:
                rec_list.append(resindex)
                res_list.append(Residue(resindex, com_str.atoms[i].residue.number,
                                        com_str.atoms[i].residue.chain, 'R', rec_index,
                                        com_str.atoms[i].residue.name,
                                        com_str.atoms[i].residue.insertion_code))
                resindex += 1
                rec_index += 1
                proc_res = res
        # save residue number in the lig list
        elif res != proc_res and resindex not in res_list:
            lig_list.append(resindex)
            res_list.append(Residue(resindex, com_str.atoms[i].residue.number,
                                    com_str.atoms[i].residue.chain, 'L', lig_index,
                                    com_str.atoms[i].residue.name,
                                    com_str.atoms[i].residue.insertion_code))
            resindex += 1
            lig_index += 1
            proc_res = res

    masks = {'REC': list2range(rec_list), 'LIG': list2range(lig_list)}

    temp = []
    for m, value in masks.items():
        for e in value['num']:
            v = e[0] if isinstance(e, list) else e
            temp.append([v, m])
    temp.sort(key=lambda x: x[0])
    order_list = [c[1] for c in temp]

    return masks, res_list, order_list


def mask2list(com_str, rec_mask, lig_mask):
    rm_list = rec_mask.strip(":").split(',')
    lm_list = lig_mask.strip(':').split(',')
    res_list = []

    for r in rm_list:
        if '-' in r:
            s, e = r.split('-')
            res_list.extend([i, 'R'] for i in range(int(s) - 1, int(e)))
        else:
            res_list.append([int(r) - 1, 'R'])
    for l in lm_list:
        if '-' in l:
            s, e = l.split('-')
            res_list.extend([i, 'L'] for i in range(int(s) - 1, int(e)))
        else:
            res_list.append([int(l) - 1, 'L'])
    res_list = sorted(res_list, key=lambda x: x[0])
    comstr = parmed.load_file(com_str)
    resl = []
    rec_index = 1
    lig_index = 1
    for res, rl in zip(comstr.residues, res_list):
        if rl[1] == 'R':
            resl.append(Residue(rl[0] + 1, res.number, res.chain, rl[1], rec_index, res.name, res.insertion_code))
            rec_index += 1
        else:
            resl.append(Residue(rl[0] + 1, res.number, res.chain, rl[1], lig_index, res.name, res.insertion_code))
            lig_index += 1
    return resl


def get_index_groups(ndx, group):
    groups = []
    with open(ndx) as ndx_file:
        groups.extend(line.split()[1] for line in ndx_file if line.startswith('['))

    if isinstance(group, int):
        if group > len(groups):
            GMXMMPBSA_ERROR('Define a valid index group')
        return group, groups[group]
    else:
        if group not in groups:
            GMXMMPBSA_ERROR('Define a valid index group')
        return groups.index(group), group


def get_indexes(com_ndx, rec_ndx=None, lig_ndx=None):
    ndx_files = {'COM': com_ndx, 'REC': rec_ndx, 'LIG': lig_ndx}
    ndx = {'COM': {'header': [], 'index': []}, 'REC': {'header': [], 'index': []}, 'LIG': {'header': [], 'index': []}}
    for n, f in ndx_files.items():
        if f is None:
            continue
        with open(f) as indexf:
            indexes = []
            for line in indexf:
                if line.startswith('['):
                    header = line.strip('\n[] ')
                    ndx[n]['header'].append(header)
                    if indexes:
                        ndx[n]['index'].append(indexes)
                        indexes = []
                else:
                    indexes.extend(map(int, line.split()))
            ndx[n]['index'].append(indexes)

    comind = ndx['COM']['header'].index('GMXMMPBSA_REC_GMXMMPBSA_LIG')
    crecind = ndx['COM']['header'].index('GMXMMPBSA_REC')
    cligind = ndx['COM']['header'].index('GMXMMPBSA_LIG')
    com_indexes = {'COM': ndx['COM']['index'][comind], 'REC': ndx['COM']['index'][crecind],
                   'LIG': ndx['COM']['index'][cligind]}
    if rec_ndx:
        recind = ndx['REC']['header'].index('GMXMMPBSA_REC')
        rec_indexes = ndx['REC']['index'][recind]
    else:
        rec_indexes = {}
    if lig_ndx:
        ligind = ndx['LIG']['header'].index('GMXMMPBSA_LIG')
        lig_indexes = ndx['LIG']['index'][ligind]
    else:
        lig_indexes = {}
    return {'COM': com_indexes, 'REC': rec_indexes, 'LIG': lig_indexes}


def get_indexes_from_str(com_ind, rec_ind=None, lig_ind=None):
    com_indexes = {'COM': [c for c, x in enumerate(com_ind['com'], start=1) if not x],
                   'REC': [c for c, x in enumerate(com_ind['rec'], start=1) if not x],
                   'LIG': [c for c, x in enumerate(com_ind['lig'], start=1) if not x]}
    rec_indexes = [c for c, x in enumerate(rec_ind, start=1) if not x] if rec_ind else {}
    lig_indexes = [c for c, x in enumerate(lig_ind, start=1) if not x] if lig_ind else {}
    return {'COM': com_indexes, 'REC': rec_indexes, 'LIG': lig_indexes}

def get_dist(coor1, coor2):
    return sqrt((coor2[0] - coor1[0]) ** 2 + (coor2[1] - coor1[1]) ** 2 + (coor2[2] - coor1[2]) ** 2)