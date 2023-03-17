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

import os
import re
import shutil
from pathlib import Path
import logging
import platform
import sys
import parmed
import re
from string import ascii_letters
import logging
from xBFreE.exceptions import GMXMMPBSA_ERROR


def create_input_args(args: list):
    if not args or 'all' in args:
        return 'general', 'gb', 'gbnsr6', 'pb', 'ala', 'nmode', 'decomp', 'rism'
    elif 'gb' not in args and 'pb' not in args and 'rism' not in args and 'nmode' not in args and 'gbnsr6' not in args:
        GMXMMPBSA_ERROR('You did not specify any type of calculation!')
    elif 'gb' not in args and 'pb' not in args and 'decomp' in args: # FIXME: gbnsr6?
        logging.warning('&decomp calculation is only compatible with &gb and &pb calculations. Will be ignored!')
        args.remove('decomp')
        return ['general'] + args
    else:
        return ['general'] + args


def remove(flag, fnpre='_GMXMMPBSA_'):
    """ Removes temporary files. Allows for different levels of cleanliness """
    # Collect all of the temporary files (those starting with _GMXMMPBSA_)
    allfiles = os.listdir(os.getcwd())

    other_files = ['COM.prmtop', 'REC.prmtop', 'LIG.prmtop', 'MUT_COM.prmtop', 'MUT_REC.prmtop', 'MUT_LIG.prmtop',
                   'leap.log']
    if flag == -1:
        result_files = ['FINAL_RESULTS_MMPBSA.dat', 'FINAL_DECOMP_MMPBSA.dat']
        for fil in allfiles:
            if (
                    fil.startswith(fnpre) or fil.startswith(f"#{fnpre}") or
                    bool(re.match('#?(COM|REC|LIG|MUT_COM|MUT_REC|MUT_LIG)_traj_(\d)\.xtc', fil)) or
                    fil == 'COMPACT_MMXSA_RESULTS.mmxsa' or
                    fil in other_files or
                    fil in result_files):
                if Path(fil).is_dir():
                    shutil.rmtree(fil)
                else:
                    os.remove(fil)

    elif flag == 0:  # remove all temporary files
        for fil in allfiles:

            if fil.startswith(fnpre) or bool(re.match('#?(COM|REC|LIG|MUT_COM|MUT_REC|MUT_LIG)_traj_(\d)\.xtc',
                                                      fil)) or fil in other_files:
                os.remove(fil)


def find_progs(INPUT, mpi_size=0):
    """ Find the necessary programs based in the user INPUT """
    # List all of the used programs with the conditions that they are needed
    logging.info('Checking external programs...')
    used_progs = {'cpptraj': True,
                  'tleap': True,
                  'parmchk2': True,
                  'sander': True,
                  'sander.APBS': INPUT['pb']['sander_apbs'] == 1,
                  'mmpbsa_py_nabnmode': INPUT['nmode']['nmoderun'],
                  # 'rism3d.snglpnt': INPUT['rism']['rismrun']
                  'elsize': INPUT['gb']['alpb'],
                  'gbnsr6': INPUT['gbnsr6']['gbnsr6run']
                  }
    gro_exe = {
        'gmx5': [
            # look for any available gromacs executable
            'gmx', 'gmx_mpi', 'gmx_d', 'gmx_mpi_d'],
        'gmx4': [
            # look for gromacs 4.x
            'make_ndx', 'trjconv', 'editconf']}

    # The returned dictionary:
    my_progs = {}

    for prog, needed in used_progs.items():
        my_progs[prog] = shutil.which(prog, path=os.environ['PATH'])
        if needed:
            if not my_progs[prog]:
                GMXMMPBSA_ERROR(f'Could not find necessary program [{prog}]')
            logging.info(f'{prog} found! Using {str(my_progs[prog])}')

    search_parth = INPUT['general']['gmx_path'] or os.environ['PATH']
    g5 = False
    for gv, g_exes in gro_exe.items():
        if gv == 'gmx5':
            for prog in g_exes:
                if exe := shutil.which(prog, path=search_parth):
                    logging.info('Using GROMACS version > 5.x.x!')
                    my_progs['make_ndx'] = [exe, 'make_ndx']
                    my_progs['editconf'] = [exe, 'editconf']
                    my_progs['trjconv'] = [exe, 'trjconv']
                    g5 = True
                    if prog in ['gmx_mpi', 'gmx_mpi_d'] and mpi_size > 1:
                        GMXMMPBSA_ERROR('gmx_mpi and gmx_mpi_d are not supported when running gmx_MMPBSA in parallel '
                                        'due to incompatibility between the mpi libraries used to compile GROMACS and '
                                        'mpi4py respectively. You can still use gmx_mpi or gmx_mpi_d to run gmx_MMPBSA '
                                        'serial. For parallel calculations use gmx instead')
                    logging.info(f'{prog} found! Using {exe}')
                    break
            if g5:
                break
        else:
            logging.info('Using GROMACS version 4.x.x!')
            for prog in g_exes:
                if exe := shutil.which(prog, path=search_parth):
                    my_progs[prog] = [exe]
                    logging.info(f'{prog} found! Using {str(my_progs[prog])}')

    if 'make_ndx' not in my_progs or 'editconf' not in my_progs or 'trjconv' not in my_progs:
        GMXMMPBSA_ERROR('Could not find necessary program [ GROMACS ]')
    logging.info('Checking external programs...Done.\n')
    return my_progs


def get_sys_info():
    """
    Print relevant system info for debugging proposes in the gmx_MMPBSA.log file
    """
    logging.debug(f"WDIR          : {Path('.').absolute().as_posix()}")
    logging.debug(f"AMBERHOME     : {os.environ['AMBERHOME'] if 'AMBERHOME' in os.environ else ''}")
    logging.debug(f"PYTHON EXE    : {shutil.which('python')}")
    logging.debug("PYTHON VERSION: " + ''.join(sys.version.split('\n')))
    logging.debug(f"MPI           : {shutil.which('mpirun')}")
    logging.debug(f"ParmEd        : {parmed.__version__}")
    logging.debug(f"OS PLATFORM   : {platform.platform()}")
    logging.debug(f"OS SYSTEM     : {platform.system()}")
    logging.debug(f"OS VERSION    : {platform.version()}")
    logging.debug(f"OS PROCESSOR  : {platform.processor()}\n")


def get_warnings():
    info = {'warning': 0, 'error': 0}
    with open('gmx_MMPBSA.log') as logfile:
        for line in logfile:
            if line.startswith('[ERROR  ]'):
                info['error'] += 1
            elif line.startswith('[WARNING]'):
                info['warning'] += 1
    return info


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
                    ci = rr[0].split(':')
                    ri = [chain, int(ci[0]), ''] if len(ci) == 1 else [chain, int(ci[0]), ci[1]]
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


def log_subprocess_output(process):
    while output := process.stdout.readline().decode():
        if output.startswith(' ->  frame'):
            continue
        logging.debug(output.strip('\n'))


def _get_dup_args(args):
    flags_values_list = []
    cv = []
    current_flag = None
    for o in args:
        if o.startswith('-'):
            if current_flag:
                flags_values_list.append([current_flag, cv])
            current_flag = o
            cv = []
        else:
            cv.append(o)

    opt_duplicates = []
    args_duplicates = []
    flags = [a[0] for a in flags_values_list]

    for x in flags:
        if flags.count(x) > 1 and x not in opt_duplicates:
            opt_duplicates.append(x)

    if opt_duplicates:
        GMXMMPBSA_ERROR('Several options are duplicated in the command-line...\n'
                        f"Duplicated options:\n\t{', '.join(opt_duplicates)}")

    flags_values_dict = dict(flags_values_list)
    unique_args = []
    inverted_args_dict = {}
    for k, v in flags_values_dict.items():
        if k in ['-cg', '-rg', '-lg']:
            continue
        for a in v:
            if a not in unique_args:
                unique_args.append(a)
                inverted_args_dict[a] = k
            else:
                args_duplicates.append([k, a])

    text_out = '\n'.join([f"\t{inverted_args_dict[a]} {' '.join(flags_values_dict[inverted_args_dict[a]])} <---> "
                          f"{k} {' '.join(flags_values_dict[k])}"
                          for k, a in args_duplicates])
    if args_duplicates:
        GMXMMPBSA_ERROR('Several args are duplicated in the command-line...\n'
                        f"Duplicated args: \n{text_out}")


class Unbuffered(object):
    """ Takes a stream handle and calls flush() on it after writing """

    def __init__(self, handle):
        self._handle = handle

    def write(self, data):
        self._handle.write(data)
        self._handle.flush()

    def __getattr__(self, attr):
        return getattr(self._handle, attr)



