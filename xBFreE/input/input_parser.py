"""
This is a module that contains functions responsible for parsing the
input file for xBFreE. It must be included with gmx_MMPBSA to
ensure proper functioning.
"""

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

from xBFreE.exceptions import InputError, InternalError
from xBFreE import __version__
import re


class Variable(object):
    """
    Base variable class. It has a name and a single value
    """

    def __init__(self, varname, dat_type=int, default=None, description='', int_dat_type=str):
        """ Initializes the variable type. Sets the default value as well as
          specifying how many characters are required by the parser to trigger
          recognition
        """
        # Catch illegalities
        if dat_type not in (int, str, float, list, tuple):
            raise InputError('Variable has unknown data type %s' % dat_type.__name__)
        if int_dat_type not in (int, str, float):
            raise InputError('Variable has unknown internal data type %s' % int_dat_type)

        self.name = varname
        self.datatype = dat_type
        self.int_datatype = int_dat_type
        if default is None:
            self.value = None
        elif self.datatype is str:
            self.value = default.replace("'", '').replace('"', '')
        elif self.datatype in [list, tuple]:
            if isinstance(default, str):
                self.value = [self.int_datatype(x.strip()) for x in re.split(';\s*|,\s*', default.replace('"',''))]
            else:
                self.value = default
        else:
            self.value = self.datatype(default)
        self.description = description

    def __str__(self):
        """ Prints statistics for the variable """
        string = 'Variable name:  %s\n' % self.name
        string += 'Variable type:  %s\n' % self.datatype
        string += 'Variable value: %s\n' % self.value
        string += 'Description:    %s\n' % self.description
        return string

    def help_str(self):
        """ returns the string [<name> = <value>.... # description] """
        if self.datatype is str:
            valstring = f'{self.name:20s} = "{self.value:s}"'
        elif self.datatype in [list, tuple]:
            if self.int_datatype == str:
                v = ','.join(self.value)
                valstring = f'{self.name:20s} = "{v}"'
            else:
                v = ','.join(map(str, self.value))
                valstring = f'{self.name:20s} = {v}'
        else:
            valstring = f'{self.name:20s} = {self.value}'
        length = 70
        valstring += ' ' + ' ' * (length - len(valstring) - 2) + ' '
        return valstring + '# %s' % self.description

    def __eq__(self, teststring):
        """ Determines if a variable string matches this variable """
        return self.name == teststring

    def __ne__(self, teststring):
        """ Not equal """
        return not self.__eq__(teststring)

    def SetValue(self, value):
        """ Sets the value of the variable """
        if self.datatype is str:
            self.value = value.replace('"', '').replace("'", '')
        elif self.datatype in [list, tuple]:
            data = value.replace('"', '').replace("'", '')
            self.value = [self.int_datatype(x.strip()) for x in re.split(';\s*|,\s*', data)]
        else:
            self.value = self.datatype(value)


class Namelist(object):
    """ Sets up a namelist. This holds many different Variables, and these
       variables will only be recognized when parsing this particular namelist.
       Set up to mimic the behavior of a Fortran namelist (so that the input is
       similar to the rest of Amber). Some of the known deficiencies:

         o the body of the namelist cannot start on the same line as the start
           or end of the namelist

         o the end of the namelist must be &end or / and must be on its own line

         o It will not (yet) recognize array lengths -- those are simply parsed
           as strings
   """

    def __init__(self, trigger, full_name, to_match=3):
        """ Sets up the list of variables, which is just the trigger for now. The
          trigger is a logical variable that gets set to true if this namelist
          is parsed. Open() trips the trigger if it exists. It can be passed in
          as anything that evaluates to False if it doesn't exist.
      """
        self.trigger = trigger
        self.variables = {}
        if self.trigger is not None:
            self.variables = {self.trigger: False}
        self.open = False
        self.full_name = full_name
        self.to_match = to_match

    def __eq__(self, nml):
        """ A namelist is equal if the name matches properly """
        return nml == self.full_name[:len(nml)] and len(nml) >= min(self.to_match, len(self.full_name))

    def __ne__(self, nml):
        """ Not equal """
        return not self.__eq__(nml)

    def addVariable(self, varname, datatype, default=None, description=None, int_dat_type=str):
        """ Adds a variable to this namelist. It checks to make sure that it's
          going to create a conflict with an existing variable.
        """
        if varname in self.variables:
            raise InternalError(f'Duplicated variable {varname} in Namelist')
        self.variables[varname] = Variable(varname, datatype, default, description, int_dat_type)

    def Open(self):
        """ Signifies that the namelist is open """
        if self.open:
            raise InputError('Namelist already open. Cannot open before closing')

        if self.trigger: self.variables[self.trigger] = True
        self.open = True

    def __str__(self):
        """
        Prints out the full contents of this namelist in the Fortran namelist
        format
        """
        retstr = '&%s\n' % self.full_name
        for variable in self.variables:
            if variable is self.trigger: continue
            retstr += '  %s\n' % self.variables[variable].help_str()
        return f'{retstr}/'


class InputFile(object):
    """ Defines the Input File and parses it. You have to add stuff to the parser
       via addNamelist. Use it as follows:

       input = InputFile()

       input.addNamelist('gb', 'gb', [['saltcon', float, 0], ...],
                         trigger='gbrun')
       input.addNamelist('ala', 'alanine_scanning', [['mutant', int, 0]],
                         trigger='alarun')

       INPUT = input.Parse('mmpbsa.in')
   """

    def __init__(self):
        """ Initializes the input file, sets up empty arrays/dictionaries """
        self.ordered_namelist_keys = []
        self.namelists = {}
        self.text = ''  # text of the given input file

    def __str__(self):
        """ Prints out the input file """
        if not self.text:
            return

        ret_text = self.text.replace('\n', '\n|')  # Add | to start of each line

        return ('|Input file:\n|---------------------------------------' +
                '-----------------------\n|' + ret_text +
                '-----------------------------------------------------' +
                '---------\n')

    def print_contents(self, destination, calc_list=None):
        """ Prints the contents of the input file """
        # Open a file to write to if need be
        # section description
        sd = {'general': '# General namelist variables',
              'gb': '# (AMBER) Generalized-Born namelist variables',
              'gbnsr6': '# GBNSR6 namelist variables',
              'pb': '# (AMBER) Possion-Boltzmann namelist variables',
              'rism': '# 3D-RISM namelist variables',
              'decomp': '# Decomposition namelist variables',
              'ala': '# Alanine scanning namelist variables',
              'nmode': '# Normal Modes Entropy namelist variables'}

        dest = destination if hasattr(destination, 'write') else open(destination, 'w')
        if calc_list:
            dest.write(f'Input file generated by gmx_MMPBSA ({__version__})\n'
                       f'Be careful with the variables you modify, some can have severe consequences on the results '
                       f'you obtain.\n\n')
        for namelist in self.ordered_namelist_keys:
            if calc_list and namelist in calc_list or not calc_list:
                dest.write(f'{sd[namelist]}\n')
                dest.write('%s\n\n' % self.namelists[namelist])

        # Close the file if we opened it.
        if dest is not destination:
            dest.close()

    def addNamelist(self, name, full_name, variable_list, trigger=None):
        """ Adds a namelist to the input file that will be parsed. Variable list
          must be an array of arrays. Each array element must be an array that
          has the information [varname, datatype, default, chars to match]. The
          'trigger' is the logical variable that gets set to true if this
          namelist is specified in the input file.
        """

        if name in self.ordered_namelist_keys:
            raise InputError('Namelist %s defined multiple times' % name)

        self.ordered_namelist_keys.append(name)
        self.namelists[name] = Namelist(trigger, full_name)

        for var in variable_list:

            if not isinstance(var, (list, tuple)) or len(var) not in [4, 5]:
                raise InputError('variables in variable_list must be lists of ' +
                                 'length 4 or 5. [varname, datatype, default, description, internal_datatype ('
                                 'Optional)]')
            if len(var) == 4:
                self.namelists[name].addVariable(var[0], var[1], var[2], var[3])
            else:
                self.namelists[name].addVariable(var[0], var[1], var[2], var[3], var[4])

    def _full_namelist_name(self, nml):
        """ Determines what the full namelist name is. We try to make as many
          allowances as possible. We will match the first 3 characters and
          replace all _'s with
      """
        nml = nml.replace(' ', '_')  # replaces spaces with _'s
        for key in self.ordered_namelist_keys:
            if self.namelists[key] == nml: return key

        raise InputError('Unrecognized namelist %s' % nml)

    def Parse(self, filename):
        """
        This subroutine parses the input file. Only data in namelists are
          parsed, and all namelists must be set prior to calling this routine.

          It will create a dictionary of Input variables for all variables in
          all namelists. They all flood the same namespace. If there are any
          conflicts between variables in namelists, an error will be raised.
          Make sure all input variables are unique!
        """
        from pathlib import Path

        # Make sure our file exists

        if filename is None:
            raise InputError("No input file was provided!")
        if not Path(filename).exists():
            raise InputError("Can't find input file (%s)" % filename)

        # Load the whole thing into memory. This should be plenty short enough.
        lines = open(filename, 'r').readlines()

        # Save the text of the input file so we can echo it back later
        self.text = ''.join(lines)

        # We will loop through the input file three times:
        #
        # 1st: Load all of the data into an array (namelist_fields)
        # 2nd: Combine multi-element values (to allow commas in input variables)
        # 3rd: Loop through the fields to change the values of the variables.

        declared_namelists = []  # keep track of the namelists we've found so far
        namelist_fields = []  # entries in a given namelist
        innml = False  # are we in a namelist now? Don't enter multiple

        # split up the input file into separate fields by comma

        for line in lines:
            # Skip title lines (we are flexible here) and comments
            if not innml and not line.strip().startswith('&'):
                continue
            if line.strip().startswith('#') or line.strip().startswith('!'):
                continue

            # Catch some errors
            if innml and line.strip().startswith('&') and line.strip() != '&end':
                raise InputError('Invalid input. Terminate each namelist prior to starting another one.')

            # End of a namelist
            elif innml and line.strip() in ['/', '&end']:
                innml = False

            # Now if we finally find a namelist
            elif not innml and line.strip().startswith('&'):
                innml = True
                namelist = line.strip()[1:].lower()
                namelist = self._full_namelist_name(namelist)

                if namelist in declared_namelists:
                    raise InputError('Namelist %s specified multiple times' % namelist)

                self.namelists[namelist].Open()
                declared_namelists.append(namelist)
                namelist_fields.append([])

            # We are in a namelist here, now fill in the fields
            elif innml:
                line = line[:line.strip().index('#')] if '#' in line else line.strip('\n')
                items = line.strip().split(',')
                # Screen any blank fields
                j = 0
                while j < len(items):
                    items[j] = items[j].strip()
                    if len(items[j]) == 0:
                        items.pop(j)
                    else:
                        j += 1
                namelist_fields[-1].extend(items)
        # # Combine any multi-element fields into the last field that has a = in it
        begin_field = -1
        for i, _ in enumerate(namelist_fields):
            for j, _ in enumerate(namelist_fields[i]):
                if '=' in namelist_fields[i][j]:
                    begin_field = j
                elif begin_field == -1:
                    raise f'Invalid input file! Error reading namelist {declared_namelists[i]}'
                else:
                    namelist_fields[i][begin_field] += f',{namelist_fields[i][j]}'
        # Now parse through the items to add them to the master dictionary. Note
        # that thanks to the last step, all data in namelist_fields will be
        # contained within fields that have a '='. All others can be ignored
        for i in range(len(namelist_fields)):
            for j in range(len(namelist_fields[i])):
                if '=' not in namelist_fields[i][j]:
                    continue
                var = namelist_fields[i][j].split('=')
                var[0] = var[0].strip()
                var[1] = var[1].strip()

                # Now we have to loop through all variables in that namelist to
                # see if this is the variable we want.
                found = False
                for key in self.namelists[declared_namelists[i]].variables:
                    if self.namelists[declared_namelists[i]].variables[key] == var[0]:
                        self.namelists[declared_namelists[i]].variables[key].SetValue(var[1])
                        found = True
                        break

                if not found:
                    raise InputError(f'Unknown variable {var[0]} in &{declared_namelists[i]}')
        # Now it's time to fill the INPUT dictionary
        INPUT = {}
        for nml in self.ordered_namelist_keys:
            INPUT[nml] = {}
            for var in self.namelists[nml].variables:
                # Here, the triggers are just bool types, so protect from accessing
                # an attribute that doesn't exist! We only allow Variable types and
                # bool types
                var_object = self.namelists[nml].variables[var]
                try:
                    INPUT[nml][var] = self.namelists[nml].variables[var].value
                except AttributeError:
                    if isinstance(var_object, bool):
                        INPUT[nml][var] = var_object
                    else:
                        raise InputError('Disallowed namelist variable type')
        return INPUT


