"""
This module contains classes that back up the state of the MM/PBSA calculation
so that future post-processing of the output files can be done without
re-supplying all of the information again.
"""

# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA                  #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

import re
import warnings
from types import SimpleNamespace

class InfoFile(object):
    """
    This is the class responsible for backing up all the calculation
    information and details so the output files can be re-written or the data
    can be readied for post-processing with external scripts.

    This object is passed a copy of the main app MMPBSA_App class
    """

    # Make a list about which INPUT variables are editable
    EDITABLE_INFO_VARS = {'general': ['exp_ki', 'sys_name', 'ie_segment', 'interaction_entropy', 'c2_entropy',
                                      'verbose'],
                          'decomp': ['csv_format', 'dec_verbose']}

    EDIT_WITH_CARE_VARS = {}

    DEPRECATED_VARS = {}

    def __init__(self, app, make_copy=False):
        """ Instantiate me """
        self.app = app
        self.make_copy = make_copy

    def write_info(self, name=None):
        """ Writes the info file to the info name """
        if name is None:
            name = 'info'
        outfile = open(name, 'w')
        # The data we have to write: INPUT, FILES, and the following attributes:
        # numframes, numframes_nmode, mpi_size (also recognize size),
        # input_file_text, mut_str, mutant_index

        # Start with INPUT (and the editable vars). Allow this to recognize INFO
        # files from the last version of gmx_MMPBSA
        outfile.write('# You can alter the variables below\n')
        for nml, nml_value in InfoFile.EDITABLE_INFO_VARS.items():
            for var in nml_value:
                outfile.write(f"INPUT['{nml}']['{var}'] = {self.write_var(self.app.INPUT[nml][var])}\n")

        outfile.write('#\n# You can alter the variables below with care\n')
        for nml, nml_value in InfoFile.EDIT_WITH_CARE_VARS.items():
            for var in nml_value:
                outfile.write(f"INPUT['{nml}']['{var}'] = {self.write_var(self.app.INPUT[nml][var])}\n")

        if InfoFile.DEPRECATED_VARS:
            outfile.write('#\n# These variables are deprecated and will be remove in next versions\n')
            for nml, nml_value in InfoFile.DEPRECATED_VARS.items():
                for var in nml_value:
                    outfile.write(f"INPUT['{nml}']['{var}'] = {self.write_var(self.app.INPUT[nml][var])}\n")

        outfile.write('#\n# MODIFY NOTHING BELOW HERE, OR GET WHAT YOU DESERVE\n')
        for nml, nml_value in self.app.INPUT.items():
            for var in nml_value:
                if (nml in InfoFile.EDITABLE_INFO_VARS and var in InfoFile.EDITABLE_INFO_VARS[nml] or
                        nml in InfoFile.EDIT_WITH_CARE_VARS and var in InfoFile.EDIT_WITH_CARE_VARS[nml] or
                nml in InfoFile.DEPRECATED_VARS and var in InfoFile.DEPRECATED_VARS[nml]):
                    continue
                outfile.write(f"INPUT['{nml}']['{var}'] = {self.write_var(self.app.INPUT[nml][var])}\n")

        # Now print out the FILES
        for var in dir(self.app.FILES):
            # Skip over __method__ functions and output files
            if var.startswith('_') or var in ('rewrite_output', 'energyout', 'overwrite'):
                continue
            outfile.write("FILES.%s = %s\n" % (var, self.write_var(getattr(self.app.FILES, var))))

        # Now print out the attributes we need to save
        outfile.write('size = %d\n' % self.app.mpi_size)
        outfile.write('numframes = %d\n' % self.app.numframes)
        outfile.write('numframes_nmode = %d\n' % self.app.numframes_nmode)
        outfile.write("mutant_index = %s\n" % self.app.mutant_index)
        outfile.write("mut_str = '%s'\n" % (self.app.resl[self.app.mutant_index].mutant_label
                                            if self.app.mutant_index is not None else ""))
        outfile.write('using_chamber = %s\n' % self.app.using_chamber)
        outfile.write(self.app.input_file_text)

    def read_info(self, name=None):
        """ Reads a _MMPBSA_info file and populates the app """
        # Give self.app a FILES attribute if it does not have one yet before
        # trying to modify it
        if not hasattr(self.app, 'FILES'):
            self.app.FILES = SimpleNamespace()
        if name is None:
            name = 'info'

        # if rewrite-output
        if self.make_copy:
            from pathlib import Path
            from datetime import datetime
            import shutil
            oldfile = Path(name)
            if oldfile.exists():
                newfile = Path(f"{name}_{datetime.now().strftime('%m%d%Y%H%M%S')}")
                shutil.copy(oldfile, newfile)

        inputre_new = re.compile(r'''INPUT\['(\S+)']\['(\S+)'] = (.*)''')
        filesre = re.compile(r'''FILES\.(\S+) = (.*)''')
        otherre = re.compile(r'(\S+) = (.*)')
        infile = open(name, 'r')
        input_text = ''
        for line in infile:
            if line.startswith('#'):
                continue
            # Now load each info file line
            if line.startswith('|'):
                input_text += line
                continue
            if rematch := inputre_new.match(line):
                nl, var, val = rematch.groups()
                val = _determine_type(val.strip())
                # print(nl, var, val)
                nl_dict = self.app.INPUT.get(nl)
                if not nl_dict:
                    self.app.INPUT[nl] = {}
                self.app.INPUT[nl][var] = val
                continue
            if rematch := filesre.match(line):
                var, val = rematch.groups()
                val = _determine_type(val.strip())
                setattr(self.app.FILES, var, val)
                continue
            if rematch := otherre.match(line):
                var, val = rematch.groups()
                val = _determine_type(val.strip())
                if var == 'size':
                    self.app.mpi_size = val
                else:
                    setattr(self.app, var, val)
                continue
        # Determine stability here:
        self.app.stability = self.app.FILES.stability
        # Set app.pre as prefix

        if self.app.FILES.receptor_trajs or self.app.FILES.ligand_trajs:
            self.app.traj_protocol = 'MT'
        else:
            self.app.traj_protocol = 'ST'
        # Load the input file text
        self.app.input_file_text = input_text

    def write_var(self, var):
        """
        Wrapper to return a string in which str vars are enclosed in quotes and
        numeric types (int and float) are not
        """
        return "'%s'" % var if isinstance(var, str) else f"{var}"


def _determine_type(thing):
    """ Determines what type this thing is """
    # If it is a string with quotes, strip off the quotes
    if thing.startswith('"') and thing.endswith('"'):
        return thing[1:-1]
    if thing.startswith("'") and thing.endswith("'"):
        return thing[1:-1]

    # Check for bool or None
    if thing == 'False':
        return False
    elif thing == 'None':
        return None

    elif thing == 'True':
        return True
    # Check for int, then check for float
    try:
        return int(thing)
    except ValueError:
        pass

    try:
        return float(thing)
    except ValueError:
        pass

    # Check for list
    if thing.startswith('[') and thing.endswith(']'):
        return eval(thing)

    # No idea what else it could be! Return string, but warn
    warnings.warn('Encountered unknown type in info file.\n' +
                  str(thing))
    return thing
