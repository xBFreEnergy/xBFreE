"""
This is a module that contains the exceptions thrown by xBFreE
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
import logging
import sys


class xBFreE_Error(Exception):
    """ Base xBFreE error class """

    def __init__(self, msg='xBFreE error'):
        self.msg = msg

    def __str__(self):
        return self.msg


class xBFreE_Warning(Warning):
    """ Base xBFreE warning class """

    def __init__(self, msg='xBFreE warning'):
        self.msg = msg

    def __str__(self):
        return self.msg


class TimerException(xBFreE_Error):
    """ Error in timer module """
    pass


class CommandlineError(xBFreE_Error):
    """ Error parsing the command-line """
    pass


class CalcError(xBFreE_Error):
    """ Error when running calculations """
    pass


class SetupError(xBFreE_Error):
    """ Error in standard setup; i.e. environment variables """
    pass


class PrmtopError(xBFreE_Error):
    """ Error in one of the prmtops """
    pass


class SelectionError(xBFreE_Error):
    """ Error in which a residue selection is illegal """
    pass


class TrajError(xBFreE_Error):
    """ Error in trajectory processing """
    pass


class InputError(xBFreE_Error):
    """ Error in the Input File """
    pass


class IllegalOption(xBFreE_Error):
    """ Error captured when looking for incompatibilities """
    pass


class InterruptError(xBFreE_Error):
    """ When the process is interrupted """
    pass


class OutputError(xBFreE_Error):
    """ Error in parsing the output """
    pass


class InconsistentError1(xBFreE_Error):
    """ Error when internal potential terms are inconsistent. Specifically
        BOND, ANGLE, and DIHED terms
    """
    pass


class InconsistentError2(xBFreE_Error):
    """ Error when internal potential terms are inconsistent. Specifically
        1-4 VDW, and 1-4 EEL terms
    """
    pass


class ConvergenceError(xBFreE_Error):
    """ Error when nmode frame doesn't minimize within tolerance """
    pass


class MutateError(xBFreE_Error):
    """ Error mutating a trajectory in alanine scanning """
    pass


class MutantResError(MutateError):
    """ Error if we're mutating to a bad residue """
    pass


class CreateInputError(xBFreE_Error):
    """ Error creating MDIN file """
    pass


class LengthError(xBFreE_Error):
    """
    Error that occurs when trying to subtract different length EnergyVectors
    """
    pass


class DecompError(xBFreE_Error):
    """ Error parsing decomp results """
    pass


class InternalError(xBFreE_Error):
    """ Error from buggy coding """
    pass


class NoFileExists(xBFreE_Error):
    """ Error if we don't have a file we need """
    pass


class InputWarning(xBFreE_Warning):
    """ If we have a non-fatal warning """
    pass

class StabilityWarning(xBFreE_Warning):
    """
    When define stability calculation and protein or ligand
    """
    pass

class xBFreEErrorLogging:
    def __init__(self, msg='xBFreE error', exc=xBFreE_Error, exception=True):
        logging.error(f"{exc.__name__} {msg}.\n           Check the xBFreE.log file to report the problem.")
        if exception:
            raise exc(msg + '\nCheck the xBFreE.log file to report the problem.')
        else:
            print(msg+ '\nCheck the xBFreE.log file to report the problem.')
            sys.exit(1)

