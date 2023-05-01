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
from pathlib import Path
import json


def mdout2json(ca, prog):
    mdout_file = Path(ca[ca.index('-o') + 1])
    output_file = mdout_file.parent.joinpath(f'{mdout_file.stem}.json')
    topology = ca[ca.index('-p') + 1]
    t = parmed.load_file(topology)
    res_list = {residue.idx + 1: [atm.idx + 1 for atm in residue.atoms] for residue in t.residues}
    pw = {x: {y: {'TDC': 0.0, 'BDC': 0.0, 'SDC': 0.0} for y in res_list} for x in res_list}

    file_assignments = []
    inputfile = []

    decomp = False
    results_section = []

    with mdout_file.open() as mmfile:
        current_section = None
        while line := mmfile.readline():
            if 'File Assignments:' in line:
                current_section = file_assignments
                line = mmfile.readline()
            elif line.startswith(' Here is the input file:'):
                current_section = inputfile
                line = mmfile.readline()
            if line.startswith('----------------------------------------------------------------------------'):
                line = mmfile.readline()
                if '.  RESULTS' in line:
                    current_section = results_section
                    mmfile.readline()
                    line = mmfile.readline()
                else:
                    current_section = None
            if current_section is not None:
                if line.startswith('DGij'):
                    pw = get_gbnsr6_out(line, pw, t)
                    decomp = True
                else:
                    current_section.append(line)
        energy = _get_energy_gbnsr6(results_section, prog)
        results = {'energy': energy}
        if decomp:
            results['decomp'] = pw
    with open(output_file, "w") as outfile:
        json.dump({'file_assignments': file_assignments, 'inputfile': inputfile, 'results_section': results}, outfile)
    mdout_file.unlink(missing_ok=True)

def get_gbnsr6_out(dgij, pw, t):
    bb = ['CA', 'C', 'O', 'N', 'H', 'OXT', 'H1', 'H2', 'H3']
    kw, at1, at2, energy = dgij.strip('\n').split()
    res_idx = t.atoms[int(at1) - 1].residue.idx + 1
    res2_idx = t.atoms[int(at2) - 1].residue.idx + 1
    if t.atoms[int(at1) - 1].name in bb:
        pw[res_idx][res2_idx]['BDC'] += float(energy)
    else:
        pw[res_idx][res2_idx]['SDC'] += float(energy)
    pw[res_idx][res2_idx]['TDC'] += float(energy)

    if res_idx != res2_idx:
        if t.atoms[int(at2) - 1].name in bb:
            pw[res2_idx][res_idx]['BDC'] += float(energy)
        else:
            pw[res2_idx][res_idx]['SDC'] += float(energy)
        pw[res2_idx][res_idx]['TDC'] += float(energy)
    return pw

def _get_energy_gbnsr6(results_section, prog):
    energy = {}

    store = False
    c = 0
    while True:
        line = results_section[c]
        if "FINAL RESULTS" in line:
            store = True
        if store and line.startswith(' 1-4 NB'):
            words = line.split()
            energy['1-4 EEL'] = float(words[7])
            c += 1
            line = results_section[c]
            words = line.split()
            energy['EEL'] = float(words[2])
            if prog in ['gbnsr6']:
                energy['EGB'] = float(words[5])
            else:
                energy['EPB'] = float(words[5])
            c += 1
            line = results_section[c]
            words = line.split()
            # FIXME: check if is CAVITY (pbsa) o not (gbnsr6), them the value is 1 or 2
            if 'CAVITY' in words[0]:
                energy['ENPOLAR'] = float(words[1])
                energy['EDISPER'] = float(words[4])
            else:
                energy[words[0].strip()] = float(words[2])
        c += 1
        if c == len(results_section):
            break
    return energy