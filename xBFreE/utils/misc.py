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
import os
import platform
import re
import shutil
import subprocess
import sys
from argparse import ArgumentParser
from pathlib import Path

import parmed
from xBFreE.exceptions import xBFreEErrorLogging


class xBFreE_ArgParser(ArgumentParser):
    """
    ArgumentParser subclass to redirect exit output to the xBFreE logging
    """

    def exit(self, status=0, message=None):
        if message:
            xBFreEErrorLogging(message)
        sys.exit(status)


def create_input_args(args: list, method):
    if method == 'mmpbsa':
        if not args or 'all' in args:
            return 'general', 'gb', 'gbnsr6', 'pb', 'ala', 'nmode', 'decomp', 'rism'
        elif 'gb' not in args and 'pb' not in args and 'rism' not in args and 'nmode' not in args and 'gbnsr6' not in args:
            xBFreEErrorLogging('You did not specify any type of calculation!')
        elif 'gb' not in args and 'pb' not in args and 'decomp' in args:  # FIXME: gbnsr6?
            logging.warning('&decomp calculation is only compatible with &gb and &pb calculations. Will be ignored!')
            args.remove('decomp')
            return ['general'] + args
        else:
            return ['general'] + args
    else:
        print('Not implemented')


def remove(method=None, keeptemp=0):
    """ Removes temporary files. Allows for different levels of cleanliness """

    result_files = {'mmpbsa': ['FINAL_RESULTS_MMPBSA.dat', 'FINAL_DECOMP_MMPBSA.dat', 'COMPACT_RESULTS_MMPBSA.xbfree']
                    # 'FINAL_RESULTS_LIE.dat', 'FINAL_DECOMP_LIE.dat',
                    # 'FINAL_RESULTS_MMPBSA.dat', 'FINAL_DECOMP_MMPBSA.dat',
                    # 'FINAL_RESULTS_MMPBSA.dat', 'FINAL_DECOMP_MMPBSA.dat',
                    }
    rfolder = Path('xBFreE_RESULTS')

    if not method:
        for mf in result_files.values():
            for f in mf:
                if Path(f).exists():
                    os.remove(f)
        if rfolder.exists():
            shutil.rmtree(rfolder)
        # remove log files
        for f in Path().glob('xBFreE*.log'):
            os.remove(f)
    else:
        if keeptemp in [-1, 0]:
            shutil.rmtree(rfolder.joinpath(method))

        if keeptemp == -1:
            for f in result_files[method]:
                if Path(f).exists():
                    os.remove(f)


def find_progs(INPUT, md_prog, mpi_size=0):
    """ Find the necessary programs based in the user INPUT """
    # List all of the used programs with the conditions that they are needed
    logging.info('Checking external programs...')
    used_progs = {'cpptraj': True,
                  'sander': True,
                  'gmx': md_prog == 'gmx',
                  # 'namd': md_prog == 'namd',
                  'sander.APBS': INPUT['pb']['sander_apbs'] == 1,
                  'mmpbsa_py_nabnmode': INPUT['nmode']['nmoderun'],
                  # 'rism3d.snglpnt': INPUT['rism']['rismrun']
                  'elsize': INPUT['gb']['alpb'],
                  'gbnsr6': INPUT['gbnsr6']['gbnsr6run'],
                  }
    # look for any available gromacs executable
    # FIXME: We should use gmx_mpi? It seems to have problem with mpi4py
    gro_exe = ['gmx', 'gmx_mpi', 'gmx_d', 'gmx_mpi_d']

    # The returned dictionary:
    my_progs = {}
    user_path = ([':'.join(INPUT['general']['exe_path'])] if ':' in INPUT['general']['exe_path']
                 else INPUT['general']['exe_path'])
    search_paths = user_path + [os.environ['PATH']]

    for prog, needed in used_progs.items():
        if needed:
            for path in search_paths:
                my_progs[prog] = shutil.which(prog, path=path)
                if my_progs[prog]:
                    break
            if not my_progs[prog]:
                xBFreEErrorLogging(f'Could not find necessary program [{prog}]')
            logging.info(f'{prog} found! Using {str(my_progs[prog])}')

    logging.info('Checking external programs...Done.\n')
    return my_progs


def get_sys_info():
    """
    Print relevant system info for debugging proposes in the xBFreE.log file
    """
    logging.debug(f"WDIR          : {Path('../mmpbsa/utils').absolute().as_posix()}")
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
    with open('xBFreE.log') as logfile:
        for line in logfile:
            if line.startswith('[ERROR  ]'):
                info['error'] += 1
            elif line.startswith('[WARNING]'):
                info['warning'] += 1
    return info


def log_subprocess_output(process):
    while output := process.stdout.readline().decode():
        if output.startswith(' ->  frame'):
            continue
        logging.debug(output.strip('\n'))


def check4dup_args(args):
    flag_index = []
    flags = []

    for i, o in enumerate(args):
        if o.startswith('-'):
            flag_index.append(i)
            flags.append(o)

    opt_duplicates = []
    flags_values = {}
    for i, f in enumerate(flags):
        if flags.count(f) > 1 and f not in opt_duplicates:
            opt_duplicates.append(f)
        if i == len(flags) - 1:
            flags_values[f] = [args[x] for x in range(flag_index[i] + 1, len(args))]
        elif flag_index[i] - flag_index[i + 1]:
            flags_values[f] = [args[x] for x in range(flag_index[i] + 1, flag_index[i + 1])]
        else:
            flags_values[f] = []

    if opt_duplicates:
        xBFreEErrorLogging('Several options are duplicated in the command-line...\n'
                           f"Duplicated options:\n\t{', '.join(opt_duplicates)}")

    args_duplicates = []
    unique_args = []
    inverted_args_dict = {}
    for k, v in flags_values.items():
        # skip this options since they can share the same group number/name
        if k in ['-cg', '-rg', '-lg']:
            continue
        for a in v:
            if a not in unique_args:
                unique_args.append(a)
                inverted_args_dict[a] = k
            else:
                args_duplicates.append([k, a])

    text_out = '\n'.join([f"\t{inverted_args_dict[a]} {' '.join(flags_values[inverted_args_dict[a]])} <---> "
                          f"{k} {' '.join(flags_values[k])}"
                          for k, a in args_duplicates])
    if args_duplicates:
        xBFreEErrorLogging('Several args are duplicated in the command-line...\n'
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
