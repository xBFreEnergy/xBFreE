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

from xBFreE.input.input_parser import InputFile

# Define the MM/PBSA input file here
input_file = InputFile()

# Add namelists with a list of variables. The variables are added to the
# namelists in lists. The entries are:
# [<variable name> <variable type> <default value> <description>]

'''
Variable structure:
|       Index       | Description                                                                                     |
|:-----------------:|:------------------------------------------------------------------------------------------------|
|   variable name   | The name of the variable. Must be easy to type, remember no space and with maximum one "_" if   |
|                   |  is needed.                                                                                     |
|   variable type   | Define the variable type though python types (int, float, str and list).  When a list is        |
|                   | defined the elements type must be defined as well. For example, list[str] define a list of      |
|                   | string elements, list[int] define a list which contain integers, and so on.                     |
|  v. description   | A brief description of the variable. This description is used to print out the inputfile and    |
|                   | terminal output for the --input-file-help                                                       |
For example:
['sys_name', str, '', 'System name'],
'''

input_file.addNamelist('general', 'general', 'General namelist variables',
                       [
                           # Basic options
                           ['sys_name', str, '', 'System name'],
                           ['startframe', int, 1, 'First frame to analyze'],
                           ['endframe', int, 9999999, 'Last frame to analyze'],
                           ['interval', int, 1, 'Number of frames between adjacent frames analyzed'],

                           # Parameters options
                           ['PBRadii', str, 'mbondi', 'Define PBRadii to build amber topology from GROMACS files'],
                           ['temperature', float, 298.15, 'Temperature'],
                           ['explicit_water', int, 0, 'Define the number of explicit water molecules'],

                           # Entropy options
                           ['qh_entropy', int, 0, 'Do quasi-harmonic calculation'],
                           ['interaction_entropy', int, 0, 'Do Interaction Entropy calculation'],
                           ['ie_segment', int, 25, 'Trajectory segment to calculate interaction entropy'],
                           ['c2_entropy', int, 0, 'Do C2 Entropy calculation'],

                           # Miscellaneous options
                           ['assign_chainID', int, 0, 'Assign chains ID'],
                           ['exp_ki', list[float], [0.0], 'Experimental Ki in nM'],
                           ['full_traj', int, 0, 'Print a full traj. AND the thread trajectories'],
                           ['gmx_path', str, '', 'Force to use this path to get GROMACS executable'],
                           ['exe_path', list[str], [], 'Add a custom path to search executables (i.e, gmx, apbs)'],
                           ['radii_path', str, '', 'Add a custom path to search radii file'],

                           ['keep_files', int, 1, 'How many files to keep after successful completion'],

                           ['netcdf', int, 0, 'Use NetCDF intermediate trajectories'],
                           ['solvated_trajectory', int, 1, 'Define if it is necessary to cleanup the trajectories'],
                           ['verbose', int, 1, 'How many energy terms to print in the final output']
                       ], trigger=None)

input_file.addNamelist('gb', 'gb', '(AMBER) Generalized-Born namelist variables',
                       [
                           ['igb', int, 5, 'GB model to use'],
                           ['intdiel', float, 1.0, 'Internal dielectric constant for sander'],
                           ['extdiel', float, 78.5, 'External dielectric constant for sander'],

                           ['saltcon', float, 0, 'Salt concentration (M)'],
                           ['surften', float, 0.0072, 'Surface tension'],
                           ['surfoff', float, 0.0, 'Surface tension offset'],
                           ['molsurf', int, 0, 'Use Connelly surface (\'molsurf\' program)'],
                           ['msoffset', float, 0.0, 'Offset for molsurf calculation'],
                           ['probe', float, 1.4, 'Solvent probe radius for surface area calc'],

                            # Options for QM
                           ['ifqnt', int, 0, 'Use QM on part of the system'],
                           ['qm_theory', str, '', 'Semi-empirical QM theory to use'],
                           ['qm_residues', str, '', 'Residues to treat with QM'],

                           # TODO: deprecated since 1.5.0. Automatic charge assignment
                           ['qmcharge_com', int, 0, 'Charge of QM region in complex'],
                           ['qmcharge_lig', int, 0, 'Charge of QM region in ligand'],
                           ['qmcharge_rec', int, 0, 'Charge of QM region in receptor'],
                           ['qmcut', float, 9999, 'Cutoff in the QM region'],
                           ['scfconv', float, 1.0e-8, 'Convergence criteria for the SCF calculation, in kcal/mol'],
                           ['peptide_corr', int, 0, 'Apply MM correction to peptide linkages'],
                           ['writepdb', int, 1, 'Write a PDB file of the selected QM region'],
                           ['verbosity', int, 0, 'Controls the verbosity of QM/MM related output'],

                           # Options for alpb
                           ['alpb', int, 0, 'Use Analytical Linearized Poisson-Boltzmann (ALPB)'],
                           ['arad_method', int, 1, 'Selected method to estimate the effective electrostatic size']
                       ], trigger='gbrun')

input_file.addNamelist('gbnsr6', 'gbnsr6', 'GBNSR6 namelist variables',
                       [
                           ['b', float, 0.028, 'Specifies the value of uniform offset to the (inverse) effective '
                                               'radii'],
                           ['alpb', int, 1, 'Specifies if ALBP correction is to be used.'],
                           ['epsin', float, 1.0, 'Sets the dielectric constant of the solute region'],
                           ['epsout', float, 78.5, 'Sets the implicit solvent dielectric constant for the solvent'],
                            # FIXME: convert to M
                           ['istrng', float, 0.0, 'Sets the ionic strength in M for the GB equation'],
                           ['rs', float, 0.52, 'Dielectric boundary shift compared to the '
                                               'molecular surface (only relevant for the -chagb option)'],
                           ['dprob', float, 1.4, 'Sets the radius of the solvent probe'],
                           ['space', float, 0.5, 'Sets the grid spacing that determines the resolution of the solute '
                                                 'molecular surface'],
                           ['arcres', float, 0.2, 'Sets the arc resolution used for numerical integration over '
                                                  'molecular surface'],
                           ['radiopt', int, 0, 'Specifies the set of intrinsic atomic radii to be used with the chagb'
                                               'option.'],
                           ['chagb', int, 0, 'Define if CHAGB is used'],
                           ['roh', int, 1, 'Sets the value of RzOH for CHA GB model'],
                           ['tau', float, 1.47, 'Sets the value of τ in the CHAGB model'],
                           ['cavity_surften', float, 0.005, 'Surface tension parameter for nonpolar '
                                                            'solvation calculation'],
                       ], trigger='gbnsr6run')

input_file.addNamelist('pb', 'pb', '(AMBER) Possion-Boltzmann namelist variables',
                       [
                            # Basic input options
                           ['ipb', int, 2, 'Dielectric model for PB'],
                           ['inp', int, 1, 'Nonpolar solvation method'],
                           ['sander_apbs', int, 0, 'Use sander.APBS?'],

                           # Options to define the physical constants
                           ['indi', float, 1, 'Internal dielectric constant'],
                           ['exdi', float, 80, 'External dielectric constant'],
                           ['emem', float, 4.0, 'Membrane dielectric constant'],
                           ['smoothopt', int, 1, 'Set up dielectric values for finite-difference grid edges that are '
                                                 'located across the solute/solvent dielectric boundary'],
                           ['istrng', float, 0.0, 'Ionic strength (M)'],
                           ['radiopt', int, 1, 'Use optimized radii?'],
                           ['prbrad', float, 1.4, 'Probe radius'],
                           ['iprob', float, 2.0, 'Mobile ion probe radius (Angstroms) for ion accessible surface used '
                                                 'to define the Stern layer'],
                           ['sasopt', int, 0, 'Molecular surface in PB implict model'],
                           ['arcres', float, 0.25, 'The resolution (Å) to compute solvent accessible arcs'],

                           # Options for Implicit Membranes
                           ['memopt', int, 0, 'Use PB optimization for membrane'],
                           ['mprob', float, 2.70, 'Membrane probe radius in Å'],
                           ['mthick', float, 40.0, 'Membrane thickness'],
                           ['mctrdz', float, 0.0, 'Distance to offset membrane in Z direction'],
                           ['poretype', int, 1, 'Use exclusion region for channel proteins'],

                           # Options to select numerical procedures
                           ['npbopt', int, 0, 'Use NonLinear PB solver?'],
                           ['solvopt', int, 1, 'Select iterative solver'],
                           ['accept', float, 0.001, 'Sets the iteration convergence criterion (relative to the initial '
                                                    'residue)'],
                           ['linit', int, 1000, 'Number of SCF iterations'],
                           ['fillratio', float, 4, 'Ratio between the longest dimension of the rectangular '
                                                   'finite-difference grid and that of the solute'],
                           ['scale', float, 2.0, '1/scale = grid spacing for the finite difference solver (default = '
                                                 '1/2 Å)'],
                           ['nbuffer', float, 0, 'Sets how far away (in grid units) the boundary of the finite '
                                                 'difference grid is away from the solute surface'],
                           ['nfocus', int, 2, 'Electrostatic focusing calculation'],
                           ['fscale', int, 8, 'Set the ratio between the coarse and fine grid spacings in an '
                                              'electrostatic focussing calculation'],
                           ['npbgrid', int, 1, 'Sets how often the finite-difference grid is regenerated'],

                            # Options to compute energy and forces
                           ['bcopt', int, 5, 'Boundary condition option'],
                           ['eneopt', int, 2, 'Compute electrostatic energy and forces'],
                           ['frcopt', int, 0, 'Output for computing electrostatic forces'],
                           ['scalec', int, 0, 'Option to compute reaction field energy and forces'],
                           ['cutfd', float, 5.0, 'Cutoff for finite-difference interactions'],
                           ['cutnb', float, 0.0, 'Cutoff for nonbonded interations'],
                           ['nsnba', int, 1, 'Sets how often atom-based pairlist is generated'],

                            # Options to select a non-polar solvation treatment
                           ['decompopt', int, 2, 'Option to select different decomposition schemes when INP = 2'],
                           ['use_rmin', int, 1, 'The option to set up van der Waals radii'],
                           ['sprob', float, 0.557, 'Solvent probe radius for SASA used to compute the dispersion term'],
                           ['vprob', float, 1.300, 'Solvent probe radius for molecular volume (the volume enclosed by '
                                                   'SASA)'],
                           ['rhow_effect', float, 1.129, 'Effective water density used in the non-polar dispersion '
                                                         'term calculation'],
                           ['use_sav', int, 1, 'Use molecular volume (the volume enclosed by SASA) for cavity term '
                                               'calculation'],
                           ['cavity_surften', float, 0.0378, 'Surface tension'],
                           ['cavity_offset', float, -0.5692, 'Offset for nonpolar solvation calc'],
                           ['maxsph', int, 400, 'Approximate number of dots to represent the maximum atomic solvent '
                                                'accessible surface'],
                           ['maxarcdot', int, 1500, 'Number of dots used to store arc dots per atom'],

                           # Options for output
                           ['npbverb', int, 0, 'Option to turn on verbose mode']
                       ], trigger='pbrun')

input_file.addNamelist('rism', 'rism', '3D-RISM namelist variables',
                       [
                           ['xvv', str, 'spc', 'XVV file for 3D-RISM'],
                           ['closure', list[str], ['kh'], 'Closure equation to use'],
                           ['gfcorrection', int, 0, 'Compute the Gaussian fluctuation excess chemical potential '
                                                    'functional'],
                           ['pcpluscorrection', int, 0, 'Compute the PC+/3D-RISM excess chemical potential functional'],
                           ['noasympcorr', int, 1, 'Turn off long range asymptotic corrections for thermodynamic '
                                                   'output only'],
                           ['buffer', float, 14, 'Distance between solute and edge of grid'],
                           ['solvcut', float, -1, 'Cutoff of the box'],
                           ['grdspc', list[float], [0.5, 0.5, 0.5], 'Grid spacing'],
                           ['ng', list[int], [-1, -1, -1], 'Number of grid points'],
                           ['solvbox', list[int], [-1, -1, -1], 'Box limits'],
                           ['tolerance', list[float], [1.0e-5], 'Convergence tolerance'],
                           ['ljTolerance', float, -1.0, 'Determines the Lennard-Jones cutoff distance based on the '
                                                        'desired accuracy of the calculation'],
                           ['asympKSpaceTolerance', float, -1.0, 'Determines the reciprocal space long range '
                                                                 'asymptotics cutoff distance based on the desired '
                                                                 'accuracy of the calculation'],
                           ['treeDCF', int, 1, 'Use the treecode approximation to calculate the direct '
                                               'correlation function (DCF) long-range asymptotic correction'],
                           ['treeTCF', int, 1, 'Use the treecode approximation to calculate the total '
                                               'correlation function (TCF) long-range asymptotic correction'],
                           ['treeCoulomb', int, 0, 'Use direct sum or the treecode approximation to calculate the '
                                                   'Coulomb potential energy'],
                           ['treeDCFMAC', float, 0.1, 'Treecode multipole acceptance criterion for the DCF long-range '
                                                      'asymptotic correction'],
                           ['treeTCFMAC', float, 0.1, 'Treecode multipole acceptance criterion for the TCF long-range '
                                                      'asymptotic correction'],
                           ['treeCoulombMAC', float, 0.1, 'Treecode multipole acceptance criterion for the Coulomb '
                                                          'potential energy'],
                           ['treeDCFOrder', int, 2, 'Treecode Taylor series order for the DCF long-range asymptotic '
                                                    'correction'],
                           ['treeTCFOrder', int, 2, 'Treecode Taylor series order for the TCF long-range asymptotic '
                                                    'correction'],
                           ['treeCoulombOrder', int, 2, 'Treecode Taylor series order for the Coulomb potential '
                                                        'energy'],
                           ['treeDCFN0', int, 500, 'Maximum number of grid points contained within the treecode leaf '
                                                   'clusters for the DCF'],
                           ['treeTCFN0', int, 500, 'Maximum number of grid points contained within the treecode leaf '
                                                   'clusters for the  TCF'],
                           ['treeCoulombN0', int, 500, 'Maximum number of grid points contained within the treecode '
                                                       'leaf clusters for the Coulomb potential energy'],
                           ['mdiis_del', float, 0.7, 'MDIIS step size'],
                           ['mdiis_nvec', int, 5, 'Number of previous iterations MDIIS uses to predict a new solution'],
                           ['mdiis_restart', float, 10.0, 'Use lowest residual solution in memory if '
                                                          'current residual is mdiis_restart times larger than '
                                                          'the smallest residual in memory'],
                           ['maxstep', int, 10000, 'Maximum number of iterative steps per solution'],
                           ['npropagate', int, 5, 'Number of previous solutions to use in predicting a new solution'],
                           ['polardecomp', int, 0, 'Break solv. energy into polar and nonpolar terms'],
                           # TODO: work with entropicDecomp? need more tests...
                           ['entropicdecomp', int, 0, 'Decomposes solvation free energy into energy and entropy '
                                                      'components'],
                           # ['centering', int, 1, 'Select how solute is centered in the solvent box'],
                           ['rism_verbose', int, 0, 'Control how much 3D-RISM info to print']
                       ], trigger='rismrun')

input_file.addNamelist('ala', 'alanine_scanning', 'Decomposition namelist variables',
                       [
                           ['mutant_res', str, '', 'Which residue will be mutated'],
                           ['mutant', str, 'ALA', 'Defines if Alanine or Glycine scanning will be performed'],
                           ['mutant_only', int, 0, 'Only compute mutant energies'],
                           ['cas_intdiel', int, 0, 'Change the intdiel value based on which aa is mutated'],
                           ['intdiel_nonpolar', int, 1, 'intdiel for nonpolar residues'],
                           ['intdiel_polar', int, 3, 'intdiel for polar residues'],
                           ['intdiel_positive', int, 5, 'intdiel for positive charged residues'],
                           ['intdiel_negative', int, 5, 'intdiel for negative charged residues']
                       ], trigger='alarun')

input_file.addNamelist('decomp', 'decomposition', 'Alanine scanning namelist variables',
                       [
                           ['idecomp', int, 0, 'Which type of decomposition analysis to do'],
                           ['dec_verbose', int, 0, 'Control energy terms are printed to the output'],
                           ['print_res', str, 'within 6', 'Which residues to print decomposition data for'],
                           ['csv_format', int, 1, 'Write decomposition data in CSV format']
                       ], trigger='decomprun')

input_file.addNamelist('nmode', 'nmode', 'Normal Modes Entropy namelist variables',
                       [
                            # Basic Options
                           ['nmstartframe', int, 1, 'First frame to analyze for normal modes'],
                           ['nmendframe', int, 1000000, 'Last frame to analyze for normal modes'],
                           ['nminterval', int, 1, 'Interval to take snapshots for normal mode analysis'],
                            # Parameters options
                           ['nmode_igb', int, 1, 'GB model to use for normal mode calculation'],
                           ['nmode_istrng', float, 0, 'Ionic strength for GB model (M)'],
                           ['dielc', float, 1, 'Dielectric constant'],
                            # Minimization options
                           ['drms', float, 0.001, 'Minimization gradient cutoff'],
                           ['maxcyc', int, 10000, 'Maximum number of minimization cycles'],
                       ], trigger='nmoderun')
