"""
 This is a module that contains the class of the main gmx_MMPBSA
 Application.
"""

# Import system modules

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

import os
import sys
import logging
from pathlib import Path
from xBFreE.utils import misc
from xBFreE.mmpbsa.output.amber import (QHout, NMODEout, QMMMout, GBout, PBout, PolarRISM_std_Out, RISM_std_Out,
                                     PolarRISM_gf_Out, RISM_gf_Out, PolarRISM_pcplus_Out, RISM_pcplus_Out,
                                     BindingStatistics, IEout, C2out, DeltaDeltaStatistics, DeltaIEC2Statistic,
                                     DeltaDeltaQH, GBNSR6out, MMout)
from xBFreE.mmpbsa.calculation import (CalculationList, EnergyCalculation, PBEnergyCalculation,
                                   NmodeCalc, QuasiHarmCalc, CopyCalc, PrintCalc, LcpoCalc, MolsurfCalc,
                                   InteractionEntropyCalc, C2EntropyCalc, MergeOut, ListEnergyCalculation)
from xBFreE.mmpbsa.createinput import create_inputs
from xBFreE.exceptions import (xBFreE_Error, InternalError, InputError, xBFreEErrorLogging)
from xBFreE.mmpbsa.infofile import InfoFile
from xBFreE.fake_mpi import MPI as FakeMPI
from xBFreE.input.mmpbsa import input_file as _input_file
from xBFreE.mmpbsa.make_trajs import make_trajectories, make_mutant_trajectories
from xBFreE.mmpbsa.output.output_file import (write_outputs, write_decomp_output, data2pkl)
from xBFreE.mmpbsa.parm_setup import MMPBSA_System
from xBFreE.utils.timer import Timer

# Global variables for the excepthook replacement at the bottom. Override these
# in the MMPBSA_App constructor and input file reading
_unbuf_stdout = misc.Unbuffered(sys.stdout)  # unbuffered stdout
_unbuf_stderr = misc.Unbuffered(sys.stderr)  # unbuffered stderr
_stdout = sys.stdout
_stderr = sys.stderr
_mpi_size = 1
_rank = 0
_MPI = FakeMPI()


# Main class

class MMPBSA_App(object):
    """ Main MM/PBSA application for driving the entire calculation """
    # The command line parser and input file objects are class attributes here
    input_file = _input_file

    def __init__(self, MPI, stdout=None, stderr=None, size=None):
        """
        Sets up the main gmx_MMPBSA driver class. All we set up here is the output
        and error streams (unbuffered by default) and the prefix for the
        intermediate files. Also set up empty INPUT dict
        """
        self.md_prog = None
        global _rank, _stdout, _stderr, _mpi_size, _MPI
        _MPI = self.MPI = MPI
        self.mmpbsa_folder = Path('xBFreE_RESULTS', 'mmpbsa')
        self.INPUT = {}
        if stdout is None:
            _stdout = self.stdout = _unbuf_stdout
        else:
            _stdout = self.stdout = stdout

        if stderr is None:
            _stderr = self.stderr = _unbuf_stderr
        else:
            _stderr = self.stderr = stderr

        # MPI-related variables. Squash output for non-master threads
        _rank = self.mpi_rank = self.MPI.COMM_WORLD.Get_rank()
        self.master = self.mpi_rank == 0
        _mpi_size = self.mpi_size = self.MPI.COMM_WORLD.Get_size()
        if not self.master:
            self.stdout = open(os.devnull, 'w')

        # create the mmpbsa folder
        if self.master:
            if self.mmpbsa_folder.exists():
                self.remove()
            self.mmpbsa_folder.mkdir()

        # Set up timers
        timers = [Timer() for _ in range(self.mpi_size)]
        self.timer = timers[self.mpi_rank]

        # Support possible threading for those that don't use MPI. However, if
        # mpi_size is > 1, just use the MPI mechanism instead
        if size is not None and self.mpi_size == 1:
            self.mpi_size = size

    def file_setup(self):
        """ Sets up the trajectories and input files """
        # If we are rewriting the output file only, bail out here
        if self.FILES.rewrite_output:
            return
        # This work belongs to the 'setup' timer
        self.timer.start_timer('setup')
        if not hasattr(self, 'normal_system'):
            xBFreEErrorLogging('MMPBSA_App not set up and parms not checked!', InternalError)
        # Set up some local refs for convenience
        FILES, INPUT, master = self.FILES, self.INPUT, self.master

        # # Now we're getting ready, remove existing intermediate files
        # elif master and not FILES.rewrite_output:
        #     self.remove(0)

        # Create input files based on INPUT dict
        if master:
            create_inputs(INPUT, self.normal_system)
        self.timer.stop_timer('setup')

        # Now create our trajectory files

        self.timer.add_timer('cpptraj', 'Creating trajectories with cpptraj:')
        self.timer.start_timer('cpptraj')

        if master:
            logging.info('Preparing trajectories for simulation...\n')
            (self.numframes, rec_frames,
             lig_frames, self.numframes_nmode) = make_trajectories(INPUT, FILES, self.external_progs['cpptraj'],
                                                                   self.mpi_size, self.devices_rank)
            if self.traj_protocol == 'MT' and not self.numframes == rec_frames == lig_frames:
                xBFreEErrorLogging('The complex, receptor, and ligand trajectories must be the same length. Since v1.5.0 '
                                'we have simplified a few things to make the code easier to maintain. Please check the '
                                'documentation')

        self.MPI.COMM_WORLD.Barrier()

        self.timer.stop_timer('cpptraj')

        self.timer.add_timer('muttraj', 'Mutating trajectories:')
        self.timer.start_timer('muttraj')

        if INPUT['ala']['alarun'] and self.master:
            logging.info('Mutating trajectories...')
        _, mutant_residue = make_mutant_trajectories(INPUT, FILES, self.mpi_rank, self.external_progs['cpptraj'],
                                                     self.normal_system, self.mutant_system)

        self.MPI.COMM_WORLD.Barrier()

        if master:
            logging.info('%d frames were processed by cpptraj for use in calculation.' % self.numframes)
            if INPUT['nmode']['nmoderun']:
                logging.info('%d frames were processed by cpptraj for nmode calculations.' % self.numframes_nmode)

        self.timer.stop_timer('muttraj')

        # Add all of the calculation timers
        self.timer.add_timer('calc', 'Total calculation time:')
        if INPUT['gb']['gbrun']:
            self.timer.add_timer('gb', 'Total GB calculation time:')
        if INPUT['gbnsr6']['gbnsr6run']:
            self.timer.add_timer('gbnsr6', 'Total GBNSR6 calculation time:')
        if INPUT['pb']['pbrun']:
            self.timer.add_timer('pb', 'Total PB calculation time:')
        if INPUT['rism']['rismrun']:
            self.timer.add_timer('rism', 'Total 3D-RISM calculation time:')
        if INPUT['nmode']['nmoderun']:
            self.timer.add_timer('nmode', 'Total normal mode calculation time:')
        if INPUT['general']['qh_entropy']:
            self.timer.add_timer('qh', 'Total quasi-harmonic calculation time:')

        self.sync_mpi()

    def run_mmpbsa(self, rank=None):
        """
        Runs the MM/PBSA analysis. This assumes FILES and INPUT are already set.
        """

        if not hasattr(self, 'external_progs'):
            xBFreEErrorLogging('external_progs not declared in run_mmpbsa!', InternalError)

        if rank is None:
            rank = self.mpi_rank
        master = rank == 0

        # Load the list of calculations we need to do, then run them.

        if master:
            self.timer.start_timer('calc')

        self.load_calc_list()

        self.stdout.write('\n')

        self.calc_list.run(rank, self.stdout)

        self.sync_mpi()

        if master:
            self.timer.stop_timer('calc')

    def load_calc_list(self):
        """
        Sets up all of the calculations to be run. When adding a new
        calculation type, add a class to calculation.py, import it at the top of
        the file here, then append it to the calc list appropriately
        """
        nframes = self.numframes if self.master else 0
        nmframes = self.numframes_nmode if self.master else 0
        self.calc_list = CalculationList(self.timer, nframes, nmframes, self.mpi_size)
        if self.master:
            gpu_size = sum([1 for d in self.devices_rank.values() if d])
            gpu_size_text = f' and {gpu_size} GPUs' if gpu_size else ''
            text = f"Starting calculations in {self.mpi_size} CPUs{gpu_size_text}..."
            logging.info(text)
            if (self.INPUT['pb']['pbrun'] or self.INPUT['rism']['rismrun'] or
                self.INPUT['nmode']['nmoderun']) and self.mpi_size > 1:
                logging.warning('PB/RISM/NMODE will be calculated with multiple threads, make sure you have enough RAM.')
        if not self.INPUT['ala']['mutant_only']:
            self.calc_list.append(PrintCalc('Running calculations on normal system...'), timer_key=None)
            self._load_calc_list(self.normal_system)
        if self.INPUT['ala']['alarun']:
            self.calc_list.append(PrintCalc('Running calculations on mutant system...'), timer_key=None)
            self._load_calc_list(self.mutant_system, True)

    def _load_calc_list(self, parm_system, mutant=False):
        """
        Internal routine to handle building calculation list. Called separately
        for mutant and normal systems
        """
        # Set up a dictionary of external programs to use based one external progs

        progs = {'gb': self.external_progs.get('sander'),
                 'gbnsr6': self.external_progs.get('gbnsr6'),
                 'sa': self.external_progs.get('cpptraj'),
                 'pb': self.external_progs.get('sander'),
                 'pbcuda': self.external_progs.get('pbcuda'),
                 'rism': self.external_progs.get('sander'),
                 'qh': self.external_progs.get('cpptraj'),
                 'nmode': self.external_progs.get('mmpbsa_py_nabnmode')
                 }
        if self.INPUT['pb']['sander_apbs']:
            progs['pb'] = self.external_progs['sander.APBS']

        # NetCDF or ASCII intermediate trajectories?
        trj_sfx = 'nc' if self.INPUT['general']['netcdf'] else 'mdcrd'

        # Determine if we just copy the receptor files. This only happens if we
        # are doing mutant calculations, we're not only doing the mutant, and the
        # receptor/mutant receptor topologies are equal. Same for the ligand
        copy_receptor = (mutant and not self.INPUT['ala']['mutant_only'] and
                         self.FILES.receptor_prmtop == self.FILES.mutant_receptor_prmtop)
        copy_ligand = (mutant and not self.INPUT['ala']['mutant_only'] and
                       self.FILES.ligand_prmtop == self.FILES.mutant_ligand_prmtop)

        prefix = 'mutant_' if mutant else ''


        mm_com_calculated = False
        mm_rec_calculated = False
        mm_lig_calculated = False

        # First load the GB calculations
        if self.INPUT['gb']['gbrun']:
            incrd = '%sdummy%%s.inpcrd' % prefix

            # See whether we are doing molsurf or LCPO. Reduce # of arguments
            # needed to 3, filling in the others here
            if self.INPUT['gb']['molsurf']:
                SAClass = lambda a1, a2, a3: MolsurfCalc(progs['sa'], a1, a2, a3,
                                                         self.INPUT['gb']['probe'], self.INPUT['gb']['msoffset'])
            else:
                SAClass = lambda a1, a2, a3: LcpoCalc(progs['sa'], a1, a2, a3,
                                                      self.INPUT['gb']['probe'], self.INPUT['gb']['msoffset'])

            # Mdin depends on decomp or not
            if self.INPUT['decomp']['decomprun']:
                mdin_template = 'gb_decomp_%s.mdin'
            elif self.INPUT['gb']['ifqnt']:
                mdin_template = 'gb_qmmm_%s.mdin'
            elif self.INPUT['gb']['alpb']:
                mdin_template = 'gb_%s.mdin'
            else:
                mdin_template = 'gb.mdin'

            # Now do complex-specific stuff
            try:
                mdin = mdin_template % 'com'
            except TypeError:
                mdin = mdin_template

            self.calc_list.append(PrintCalc(f"Beginning GB calculations with {progs['gb']}"), timer_key='gb')

            c = EnergyCalculation(progs['gb'], parm_system.complex_prmtop,
                                  incrd % 'complex',
                                  '%scomplex.%s.%%d' % (prefix, trj_sfx),
                                  mdin, '%scomplex_gb.mdout.%%d' % (prefix),
                                  'restrt.%d')
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='gb',
                                  output_basename='%scomplex_gb.mdout.%%d' % (prefix))
            c = SAClass(parm_system.complex_prmtop,
                        '%scomplex.%s.%%d' % (prefix, trj_sfx),
                        '%scomplex_gb_surf.dat.%%d' % prefix)
            self.calc_list.append(c, '', timer_key='gb')

            if not self.stability:
                try:
                    mdin = mdin_template % 'rec'
                except TypeError:
                    mdin = mdin_template

                # Either copy the existing receptor if the mutation is in the ligand
                # or perform a receptor calculation
                if copy_receptor:
                    c = CopyCalc('receptor_gb.mdout.%%d', '%sreceptor_gb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='gb')
                    c = CopyCalc('receptor_gb_surf.dat.%%d', '%sreceptor_gb_surf.dat.%%d' % prefix)
                    self.calc_list.append(c, '', timer_key='gb')
                else:
                    c = EnergyCalculation(progs['gb'], parm_system.receptor_prmtop,
                                          incrd % 'receptor',
                                          '%sreceptor.%s.%%d' % (prefix, trj_sfx),
                                          mdin, '%sreceptor_gb.mdout.%%d' % (prefix),
                                          'restrt.%d')
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='gb', output_basename='%sreceptor_gb.mdout.%%d' % (prefix))
                c = SAClass(parm_system.receptor_prmtop,
                            '%sreceptor.%s.%%d' % (prefix, trj_sfx),
                            '%sreceptor_gb_surf.dat.%%d' % prefix)
                self.calc_list.append(c, '', timer_key='gb')

                try:
                    mdin = mdin_template % 'lig'
                except TypeError:
                    mdin = mdin_template

                # Either copy the existing ligand if the mutation is in the receptor
                # or perform a ligand calculation
                if copy_ligand:
                    c = CopyCalc('ligand_gb.mdout.%%d', '%sligand_gb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='gb')
                    c = CopyCalc('ligand_gb_surf.dat.%%d', '%sligand_gb_surf.dat.%%d' % prefix)
                    self.calc_list.append(c, '', timer_key='gb')
                else:
                    c = EnergyCalculation(progs['gb'], parm_system.ligand_prmtop,
                                          incrd % 'ligand',
                                          '%sligand.%s.%%d' % (prefix, trj_sfx),
                                          mdin, '%sligand_gb.mdout.%%d' % (prefix),
                                          'restrt.%d')
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='gb', output_basename='%sligand_gb.mdout.%%d' % (prefix))
                c = SAClass(parm_system.ligand_prmtop,
                            '%sligand.%s.%%d' % (prefix, trj_sfx),
                            '%sligand_gb_surf.dat.%%d' % prefix)
                self.calc_list.append(c, '', timer_key='gb')
        # end if self.INPUT['gb']['gbrun']

        # GBNSR6 calculation
        if self.INPUT['gbnsr6']['gbnsr6run']:
            incrd = '%sdummy%%s.inpcrd' % prefix
            mdin = 'gbnsr6.mdin'

            # Mdin depends on decomp or not
            if self.INPUT['decomp']['decomprun']:
                mdin_template = 'mm_gbnsr6_decomp_%s.mdin'
            else:
                mdin_template = 'mm_gbnsr6.mdin'

            # Now do complex-specific stuff
            try:
                mm_mdin = mdin_template % 'com'
            except TypeError:
                mm_mdin = mdin_template


            self.calc_list.append(PrintCalc(f"Beginning GBNSR6 calculations with {progs['gbnsr6']}"),
                                  timer_key='gbnsr6')
            if not mm_com_calculated:
                self.calc_list.append(PrintCalc("  calculating complex contribution..."),
                                      timer_key='gbnsr6')
                c = EnergyCalculation(progs['gb'], parm_system.complex_prmtop,
                                      incrd % 'complex',
                                      '%scomplex.%s.%%d' % (prefix, trj_sfx),
                                      mm_mdin,
                                      f'{prefix}complex_mm.mdout.%d',
                                      'restrt.%d')
                self.calc_list.append(c, '    calculating MM...', timer_key='gbnsr6',
                                      output_basename=f'{prefix}complex_mm.mdout.%d')
            # use pre directly to have only one folder per rank
            files = sorted(list(Path(f"inpcrd_{self.mpi_rank}").glob(f"{prefix}complex*.inpcrd")),
                           key=lambda x: int(x.stem.split('.')[1]))
            mdouts = [file.parent.joinpath(f"{file.name.split('.')[0]}_gbnsr6{file.suffixes[0]}.mdout").as_posix()
                      for file in files]
            inpcrds = [file.as_posix() for file in files]

            c = ListEnergyCalculation(progs['gbnsr6'], parm_system.complex_prmtop, mdin, inpcrds, mdouts)
            self.calc_list.append(c, '    calculating GB...', timer_key='gbnsr6',
                                      output_basename=f"inpcrd_%d/{prefix}complex_gbnsr6.mdout")

            c = MergeOut(self.FILES.complex_prmtop, f"{prefix}complex_gbnsr6.mdout.%d",
                         f'{prefix}complex_mm.mdout.%d', mdouts, self.INPUT['decomp']['idecomp'],
                         self.INPUT['decomp']['dec_verbose'])
            self.calc_list.append(c, '', timer_key='gbnsr6')

            if not self.stability:
                try:
                    mm_mdin = mdin_template % 'rec'
                except TypeError:
                    mm_mdin = mdin_template
                # Either copy the existing receptor if the mutation is in the ligand
                # or perform a receptor calculation
                if copy_receptor:
                    c = CopyCalc('receptor_mm.mdout.%d',
                                 f'{prefix}receptor_mm.mdout.%d')
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='gbnsr6')
                    c = CopyCalc('receptor_gbnsr6.mdout.%d',
                                 f'{prefix}receptor_gbnsr6.mdout.%d')
                    self.calc_list.append(c, '', timer_key='gbnsr6')
                else:
                    if not mm_rec_calculated:
                        self.calc_list.append(PrintCalc("  calculating receptor contribution..."),
                                              timer_key='gbnsr6')
                        c = EnergyCalculation(progs['gb'], parm_system.receptor_prmtop,
                                              incrd % 'receptor',
                                              f'{prefix}receptor.{trj_sfx}.%d',
                                              mm_mdin,
                                              f'{prefix}receptor_mm.mdout.%d',
                                              'restrt.%d')

                        self.calc_list.append(c, '    calculating MM...', timer_key='gbnsr6',
                                              output_basename=f'{prefix}receptor_mm.mdout.%d')
                    files = sorted(list(Path(f"inpcrd_{self.mpi_rank}").glob(f"{prefix}receptor*.inpcrd")),
                                   key=lambda x: int(x.stem.split('.')[1]))
                    mdouts = [
                        file.parent.joinpath(f"{file.name.split('.')[0]}_gbnsr6{file.suffixes[0]}.mdout").as_posix()
                        for file in files]
                    inpcrds = [file.as_posix() for file in files]

                    c = ListEnergyCalculation(progs['gbnsr6'], parm_system.receptor_prmtop, mdin, inpcrds, mdouts)
                    self.calc_list.append(c, '    calculating GB...', timer_key='gbnsr6',
                                          output_basename=f"inpcrd_%d/{prefix}receptor_gbnsr6.mdout")

                    c = MergeOut(self.FILES.receptor_prmtop, f"{prefix}receptor_gbnsr6.mdout.%d",
                                 f'{prefix}receptor_mm.mdout.%d', mdouts, self.INPUT['decomp']['idecomp'],
                                 self.INPUT['decomp']['dec_verbose'])
                    self.calc_list.append(c, '', timer_key='gbnsr6')

                try:
                    mm_mdin = mdin_template % 'lig'
                except TypeError:
                    mm_mdin = mdin_template
                # Either copy the existing ligand if the mutation is in the receptor
                # or perform a ligand calculation
                if copy_ligand:
                    c = CopyCalc('ligand_mm.mdout.%d',
                                 f'{prefix}ligand_mm.mdout.%d')
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='gbnsr6')
                    c = CopyCalc('ligand_gbnsr6.mdout.%d',
                                 f'{prefix}ligand_gbnsr6.mdout.%d')
                    self.calc_list.append(c, '', timer_key='gbnsr6')
                else:
                    if not mm_lig_calculated:
                        self.calc_list.append(PrintCalc("  calculating ligand contribution..."),
                                              timer_key='gbnsr6')
                        c = EnergyCalculation(progs['gb'], parm_system.ligand_prmtop,
                                              incrd % 'ligand',
                                              f'{prefix}ligand.{trj_sfx}.%d',
                                              mm_mdin,
                                              f'{prefix}ligand_mm.mdout.%d',
                                              'restrt.%d')

                        self.calc_list.append(c, '    calculating MM...', timer_key='gbnsr6',
                                              output_basename=f'{prefix}ligand_mm.mdout.%d')
                    files = sorted(list(Path(f"inpcrd_{self.mpi_rank}").glob(f"{prefix}ligand*.inpcrd")),
                                   key=lambda x: int(x.stem.split('.')[1]))
                    mdouts = [
                        file.parent.joinpath(f"{file.name.split('.')[0]}_gbnsr6{file.suffixes[0]}.mdout").as_posix()
                        for file in files]
                    inpcrds = [file.as_posix() for file in files]

                    c = ListEnergyCalculation(progs['gbnsr6'], parm_system.ligand_prmtop, mdin, inpcrds, mdouts)
                    self.calc_list.append(c, '    calculating GB...', timer_key='gbnsr6',
                                          output_basename=f"inpcrd_%d/{prefix}ligand_gbnsr6.mdout")
                    c = MergeOut(self.FILES.ligand_prmtop, f"{prefix}ligand_gbnsr6.mdout.%d",
                                 f'{prefix}ligand_mm.mdout.%d', mdouts, self.INPUT['decomp']['idecomp'],
                                 self.INPUT['decomp']['dec_verbose'])
                    self.calc_list.append(c, '', timer_key='gbnsr6')

        # Next load the PB calculations
        if self.INPUT['pb']['pbrun']:
            incrd = '%sdummy%%s.inpcrd' % prefix

            # Mdin depends on decomp or not
            if self.INPUT['decomp']['decomprun']:
                mdin_template = 'pb_decomp_%s.mdin'
            else:
                mdin_template = 'pb.mdin'

            # Now do complex-specific stuff
            try:
                mdin = mdin_template % 'com'
            except TypeError:
                mdin = mdin_template

            self.calc_list.append(PrintCalc(f"Beginning PB calculations with {progs['pb']}"), timer_key='pb')

            c = PBEnergyCalculation(progs['pb'], parm_system.complex_prmtop,
                                    incrd % 'complex',
                                    '%scomplex.%s.%%d' % (prefix, trj_sfx),
                                    mdin, '%scomplex_pb.mdout.%%d' % prefix,
                                    'restrt.%d')
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='pb',
                                  output_basename='%scomplex_pb.mdout.%%d' % (prefix))
            if not self.stability:
                try:
                    mdin = mdin_template % 'rec'
                except TypeError:
                    mdin = mdin_template

                # Either copy the existing receptor if the mutation is in the ligand
                # or perform a receptor calculation
                if copy_receptor:
                    c = CopyCalc('receptor_pb.mdout.%%d', '%sreceptor_pb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = PBEnergyCalculation(progs['pb'], parm_system.receptor_prmtop,
                                            incrd % 'receptor',
                                            '%sreceptor.%s.%%d' % (prefix, trj_sfx),
                                            mdin, '%sreceptor_pb.mdout.%%d' % prefix,
                                            'restrt.%d')
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='pb', output_basename='%sreceptor_pb.mdout.%%d' % (prefix))

                try:
                    mdin = mdin_template % 'lig'
                except TypeError:
                    mdin = mdin_template

                # Either copy the existing ligand if the mutation is in the receptor
                # or perform a ligand calculation
                if copy_ligand:
                    c = CopyCalc('ligand_pb.mdout.%%d', '%sligand_pb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = PBEnergyCalculation(progs['pb'], parm_system.ligand_prmtop,
                                            incrd % 'ligand',
                                            '%sligand.%s.%%d' % (prefix, trj_sfx),
                                            mdin, '%sligand_pb.mdout.%%d' % (prefix),
                                            'restrt.%d')
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='pb', output_basename='%sligand_pb.mdout.%%d' % (prefix))
        # end if self.INPUT['pb']['pbrun']

        if self.INPUT['rism']['rismrun']:
            mdin = 'rism.mdin'
            # get xvv file from INPUT
            from xBFreE.mmpbsa.data import xvv_files

            self.FILES.xvvfile = None
            if self.INPUT['rism']['xvv'] in xvv_files:
                self.FILES.xvvfile = xvv_files[self.INPUT['rism']['xvv']]
            else:
                p = Path(self.INPUT['rism']['xvv'])
                if p.exists() and p.suffix == '.xvv':
                    self.FILES.xvvfile = p.absolute().as_posix()

            if not self.FILES.xvvfile:
                xBFreEErrorLogging('Please, define a valid xvv file')

            self.calc_list.append(
                PrintCalc('Beginning 3D-RISM calculations with %s' % progs['rism']), timer_key='rism')

            c = EnergyCalculation(progs['rism'], parm_system.complex_prmtop,
                                  '%sdummycomplex.inpcrd' % prefix,
                                  '%scomplex.%s.%%d' % (prefix, trj_sfx), mdin,
                                  '%scomplex_rism.mdout.%%d' % prefix,
                                  'restrt.%d', self.FILES.xvvfile)
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='rism',
                                  output_basename='%scomplex_rism.mdout.%%d' % (prefix))

            if not self.stability:
                if copy_receptor:
                    c = CopyCalc('receptor_rism.mdout.%%d', '%sreceptor_rism.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = EnergyCalculation(progs['rism'], parm_system.receptor_prmtop,
                                          '%sdummyreceptor.inpcrd' % prefix,
                                          '%sreceptor.%s.%%d' % (prefix, trj_sfx), mdin,
                                          '%sreceptor_rism.mdout.%%d' % prefix,
                                          'restrt.%d', self.FILES.xvvfile)
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='rism', output_basename='%sreceptor_rism.mdout.%%d' % (prefix))

                if copy_ligand:
                    c = CopyCalc('ligand_rism.mdout.%%d', '%sligand_rism.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = EnergyCalculation(progs['rism'], parm_system.ligand_prmtop,
                                          '%sdummyligand.inpcrd' % prefix,
                                          '%sligand.%s.%%d' % (prefix, trj_sfx), mdin,
                                          '%sligand_rism.mdout.%%d' % prefix,
                                          'restrt.%d', self.FILES.xvvfile)
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='rism', output_basename='%sligand_rism.mdout.%%d' % (prefix))

        if self.INPUT['nmode']['nmoderun']:
            self.calc_list.append(
                PrintCalc('Beginning nmode calculations with %s' % progs['nmode']), timer_key='nmode')

            c = NmodeCalc(progs['nmode'], parm_system.complex_prmtop,
                          '%scomplex.pdb' % prefix,
                          '%scomplex_nm.%s.%%d' % (prefix, trj_sfx),
                          '%scomplex_nm.out.%%d' % prefix, self.INPUT)
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='nmode',
                                  output_basename='%scomplex_nm.out.%%d' % (prefix))

            if not self.stability:
                if copy_receptor:
                    c = CopyCalc('receptor_nm.out.%%d', '%sreceptor_nm.out.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='nmode')
                else:
                    c = NmodeCalc(progs['nmode'], parm_system.receptor_prmtop,
                                  '%sreceptor.pdb' % prefix,
                                  '%sreceptor_nm.%s.%%d' % (prefix, trj_sfx),
                                  '%sreceptor_nm.out.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='nmode', output_basename='%sreceptor_nm.out.%%d' % (prefix))

                if copy_ligand:
                    c = CopyCalc('ligand_nm.out.%%d', '%sligand_nm.out.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='nmode')
                else:
                    c = NmodeCalc(progs['nmode'], parm_system.ligand_prmtop,
                                  '%sligand.pdb' % prefix,
                                  '%sligand_nm.%s.%%d' % (prefix, trj_sfx),
                                  '%sligand_nm.out.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='nmode', output_basename='%sligand_nm.out.%%d' % (prefix))

        # Only master does entropy calculations
        if self.INPUT['general']['qh_entropy']:
            self.calc_list.append(
                PrintCalc('Beginning quasi-harmonic calculations with %s' %
                          progs['qh']), timer_key='qh')

            c = QuasiHarmCalc(progs['qh'], parm_system.complex_prmtop,
                              '%scomplex.%s' % (prefix, trj_sfx),
                              '%scpptrajentropy.in' % prefix,
                              '%scpptraj_entropy.out' % prefix,
                              self.INPUT['general']['receptor_mask'],
                              self.INPUT['general']['ligand_mask'])
            self.calc_list.append(c, '', timer_key='qh')

    def make_prmtops(self):
        self.timer.add_timer('setup_top', 'Total GROMACS setup time:')
        self.timer.start_timer('setup_top')
        # if not self.FILES.rewrite_output and self.master:
        #     self.remove(-1)

        # Find external programs IFF we are doing a calc
        external_progs = misc.find_progs(self.INPUT, self.md_prog, self.mpi_size) if self.master else {}
        external_progs = self.MPI.COMM_WORLD.bcast(external_progs, root=0)
        # Make external_progs an instance attribute
        self.external_progs = external_progs
        if self.master:
            # Make amber topologies
            if self.md_prog == 'gmx':
                logging.info(f'Building AMBER topologies from GROMACS files using {self.traj_protocol} approach...')
                from xBFreE.mmpbsa.topology.gmx import BuildTopGromacs
                top_builder = BuildTopGromacs
            elif self.md_prog == 'amber':
                logging.info(f'Building AMBER topologies from AMBER files using {self.traj_protocol} approach...')
                from xBFreE.mmpbsa.topology.amber import BuildTopAmber
                top_builder = BuildTopAmber
            elif self.md_prog == 'namd':
                logging.info(f'Building AMBER topologies from NAMD files using {self.traj_protocol} approach...')
                if self.FILES.complex_top.split('.')[-1] == 'psf':
                    from xBFreE.mmpbsa.topology.namd import BuildTopNAMD_CHARMM
                    top_builder = BuildTopNAMD_CHARMM
                else:
                    from xBFreE.mmpbsa.topology.namd import BuildTopNAMD_Amber
                    top_builder = BuildTopNAMD_Amber
            else:
                logging.info(f'Building AMBER topologies from CHARMM files using {self.traj_protocol} approach...')
                from xBFreE.mmpbsa.topology.charmm import BuildTopCHARMM
                top_builder = BuildTopCHARMM

            maketop = top_builder(self.FILES, self.INPUT, self.external_progs)
            (self.FILES.complex_prmtop, self.FILES.receptor_prmtop, self.FILES.ligand_prmtop,
             self.FILES.mutant_complex_prmtop,
             self.FILES.mutant_receptor_prmtop, self.FILES.mutant_ligand_prmtop) = maketop.buildTopology()
            logging.info('Building AMBER topologies from GROMACS files... Done.\n')
            self.INPUT['general']['receptor_mask'], self.INPUT['general']['ligand_mask'], self.resl = maketop.get_masks()
            self.mutant_index = maketop.com_mut_index
            self.mut_str = self.resl[maketop.com_mut_index].mutant_label if self.mutant_index is not None else ''
            # self.FILES.complex_fixed = 'COM_FIXED.pdb'
        self.FILES = self.MPI.COMM_WORLD.bcast(self.FILES, root=0)
        self.INPUT = self.MPI.COMM_WORLD.bcast(self.INPUT, root=0)
        self.sync_mpi()
        self.timer.stop_timer('setup_top')

    def loadcheck_prmtops(self):
        """ Loads the topology files and checks their consistency """
        # Start setup timer and make sure we've already set up our input
        self.timer.add_timer('setup', 'Total AMBER setup time:')
        self.timer.start_timer('setup')
        if not hasattr(self, 'FILES') or not hasattr(self, 'INPUT'):
            xBFreEErrorLogging('MMPBSA_App not set up! Cannot check parms yet!', InternalError)
        # create local aliases to avoid abundant selfs
        FILES, INPUT = self.FILES, self.INPUT
        if self.master:
            # Now load the parms and check them
            logging.info('Loading and checking parameter files for compatibility...')
        self.normal_system = MMPBSA_System(FILES.complex_prmtop, FILES.receptor_prmtop, FILES.ligand_prmtop)
        self.using_chamber = self.normal_system.complex_prmtop.chamber
        self.mutant_system = None
        if INPUT['ala']['alarun']:
            if (FILES.mutant_receptor_prmtop is None and FILES.mutant_ligand_prmtop is None and not self.stability):
                xBFreEErrorLogging('Alanine scanning requires either a mutated receptor or mutated ligand topology '
                                'file!')
            if FILES.mutant_receptor_prmtop is None:
                FILES.mutant_receptor_prmtop = FILES.receptor_prmtop
            elif FILES.mutant_ligand_prmtop is None:
                FILES.mutant_ligand_prmtop = FILES.ligand_prmtop
            self.mutant_system = MMPBSA_System(FILES.mutant_complex_prmtop, FILES.mutant_receptor_prmtop,
                                               FILES.mutant_ligand_prmtop)
        # If we have a chamber prmtop, force using sander
        if self.using_chamber:
            if INPUT['rism']['rismrun']:
                xBFreEErrorLogging('CHAMBER prmtops cannot be used with 3D-RISM')
            if INPUT['nmode']['nmoderun']:
                xBFreEErrorLogging('CHAMBER prmtops cannot be used with NMODE')

        self.normal_system.Map(INPUT['general']['receptor_mask'], INPUT['general']['ligand_mask'])
        self.normal_system.CheckConsistency()
        if INPUT['ala']['alarun']:
            self.mutant_system.Map(INPUT['general']['receptor_mask'], INPUT['general']['ligand_mask'])
            self.mutant_system.CheckConsistency()
        if (INPUT['general']['ligand_mask'] is None or INPUT['general']['receptor_mask'] is None):
            com_mask, INPUT['general']['receptor_mask'], INPUT['general']['ligand_mask'] = \
                self.normal_system.Mask('all', in_complex=True)
        self.sync_mpi()
        self.timer.stop_timer('setup')

    def write_final_outputs(self):
        """ Writes the final output files for gmx_MMPBSA """
        self.timer.add_timer('output', 'Statistics calculation & output writing:')
        self.timer.start_timer('output')
        if (not hasattr(self, 'input_file_text') or not hasattr(self, 'FILES') or
                not hasattr(self, 'INPUT') or not hasattr(self, 'normal_system')):
            xBFreEErrorLogging('I am not prepared to write the final output file!', InternalError)
        # Only the master does this, so bail out if we are not master
        if not self.master:
            return
        # If we haven't already parsed our output files, do that now
        # FIXME: does this make sense?
        if not hasattr(self, 'calc_types'):
            self.parse_output_files()
        # Do the output files now
        write_outputs(self)
        if self.INPUT['decomp']['decomprun']:
            write_decomp_output(self)
        if self.INPUT['general']['keep_files'] in [0, 2]:
            data2pkl(self)

        info = InfoFile(self)
        info.write_info('info')

        os.chdir(self.FILES.wdir)

        self.timer.stop_timer('output')

    def finalize(self):
        """ We are done. Finish up timers and print out timing info """
        self.timer.done()
        if not self.master:
            self.MPI.Finalize()
            sys.exit(0)
        logging.info('Timing:')
        if not self.FILES.rewrite_output:
            self.timer.print_('setup_top')
        self.timer.print_('setup')

        if not self.FILES.rewrite_output:
            self._finalize_timers()
        self.timer.print_('output')
        self.timer.print_('global', True)

        self.remove(self.INPUT['general']['keep_files'])

        exe_info = misc.get_warnings()
        logging.info(f"\n   Finalizing gmx_MMPBSA: [ERROR  ] = {exe_info['error']}; [WARNING] = {exe_info['warning']}\n"
                     f"   Check the gmx_MMPBSA.log file for more details...\n")


        logging.info(
            '\n Thank you for using gmx_MMPBSA. Please consider supporting gmx_MMPBSA by citing our publication:'
            '\n    Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. '
            '\n    gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. '
            '\n    J Chem Theory Comput., 2021, 17 (10):6281-6291. Epub 2021 Sep 29. PMID: 34586825.'
            '\n    https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645'
            '\n\nAlso consider citing MMPBSA.py:'
            '\n    Miller III, B. R., McGee Jr., T. D., Swails, J. M. Homeyer, N. Gohlke, H. and Roitberg, A. E.'
            '\n    MMPBSA.py: An Efficient Program for End-State Free Energy Calculations.'
            '\n    J. Chem. Theory Comput., 2012, 8 (9) pp 3314-3321\n')
        self.MPI.Finalize()

        end = 0
        if self.FILES.gui and not self.FILES.rewrite_output:
            import subprocess
            from pathlib import Path
            logging.info('Opening gmx_MMPBSA_ana to analyze results...\n')
            ifile = Path(f'{self.FILES.prefix}info')
            if not ifile.exists():
                ifile = Path('COMPACT_RESULTS_MMPBSA.xbfree')

            g = subprocess.Popen(['gmx_MMPBSA_ana', '-f', ifile.as_posix()])
            if g.wait():
                end = 1
        if end:
            logging.error('Unable to start gmx_MMPBSA_ana...')
        logging.info('Finalized...')
        sys.exit(end)

    def _finalize_timers(self):
        self.timer.print_('cpptraj')

        if self.INPUT['ala']['alarun']:
            self.timer.print_('muttraj', True)

        # self.stdout.write('\n')
        self.timer.print_('calc')

        if self.INPUT['gb']['gbrun']:
            self.timer.print_('gb')

        if self.INPUT['gbnsr6']['gbnsr6run']:
            self.timer.print_('gbnsr6')

        if self.INPUT['pb']['pbrun']:
            self.timer.print_('pb')

        if self.INPUT['nmode']['nmoderun']:
            self.timer.print_('nmode')

        if self.INPUT['general']['qh_entropy']:
            self.timer.print_('qh', True)

        # self.stdout.write('\n')

    def get_parser(self, parser):
        if self.master:
            self.FILES = parser
        else:
            self.FILES = object()
        os.chdir(self.FILES.subwdir)
        # Broadcast the FILES
        self.FILES = self.MPI.COMM_WORLD.bcast(self.FILES)
        if self.FILES.receptor_trajs or self.FILES.ligand_trajs:
            self.traj_protocol = 'MT'  # multiple traj protocol
        else:
            self.traj_protocol = 'ST'  # single traj protocol
        # change by explicit argument
        self.stability = self.FILES.stability

        # define md_prog to identify file type and workflow
        if 'gmx' in self.FILES.subparser:
            self.md_prog = 'gmx'
        elif 'amber' in self.FILES.subparser:
            self.md_prog = 'amber'
        elif 'namd' in self.FILES.subparser:
            self.md_prog = 'namd'
        elif 'charmm' in self.FILES.subparser:
            self.md_prog = 'charmm'
        else:
            xBFreEErrorLogging('Unknown MD program!')


    def read_input_file(self, infile=None):
        """ Reads the input file, pull it from FILES if not provided here """
        if infile is None:
            if not hasattr(self, 'FILES'):
                xBFreEErrorLogging('FILES not present and no input file given!', InternalError)
            infile = self.FILES.input_file
        self.INPUT = self.input_file.Parse(infile)
        self.input_file_text = str(self.input_file)
        if self.master:
            for line in self.input_file_text.split('\n'):
                logging.debug(line)

    def process_input(self):
        """
        This handles processing of the INPUT dict if necessary if this is a 'new'
        calculation (i.e., not re-writing output). This does the following prep:
           - invert scale
           - determine trajectory file suffix
           - set decomp-dependent GBSA default
           - adjust verbose for stability calcs
           - 3D-RISM setup
           - Set temperature. Don't put it in namelist, because temp change
             for entropy requires changes to nmode and cpptraj calcs, meaning it
             is not as easily changed here.
        """
        # Invert scale
        self.INPUT['pb']['scale'] = 1 / self.INPUT['pb']['scale']

        # Set up netcdf variables and decide trajectory suffix
        if self.INPUT['general']['netcdf'] == 0:
            self.INPUT['general']['netcdf'] = ''
            self.trjsuffix = 'mdcrd'
        else:
            self.INPUT['general']['netcdf'] = 'netcdf'
            self.trjsuffix = 'nc'

        # Set default GBSA for Decomp
        if self.INPUT['decomp']['decomprun']:
            self.INPUT['gb']['gbsa'] = 2

        # Stability: no terms cancel, so print them all
        # if self.stability:
        #     self.INPUT['verbose'] = 2

        # 3D-RISM stuff (keywords are case-insensitive)
        if self.INPUT['rism']['solvcut'] is None:
            self.INPUT['rism']['solvcut'] = self.INPUT['rism']['buffer']

        self.INPUT['rism']['rismrun_std'] = bool(self.INPUT['rism']['rismrun'])
        self.INPUT['rism']['rismrun_gf'] = self.INPUT['rism']['rismrun'] and self.INPUT['rism']['gfcorrection']
        self.INPUT['rism']['rismrun_pcplus'] = self.INPUT['rism']['rismrun'] and self.INPUT['rism']['pcpluscorrection']

        # Default temperature
        # self.INPUT['temp'] = 298.15

        # check for cuda if any calculation require it
        if self.INPUT['pbcuda']['pbcudarun']:
            from xBFreE.utils.cuda_check import get_cuda_devices
            # create new communicators for each node according to the rank type shared (for size > 1)
            if self.mpi_size > 1:
                nodes_comm = self.MPI.COMM_WORLD.Split_type(self.MPI.COMM_TYPE_SHARED, key=self.mpi_rank)
                nodes_rank = nodes_comm.Get_rank()
                if nodes_rank == 0:
                    dds = [self.mpi_rank, {'nodes_rank': nodes_rank, 'host_name': self.MPI.Get_processor_name(),
                                           'devices': get_cuda_devices()}]
                else:
                    dds = [self.mpi_rank]
                self.devices = self.MPI.COMM_WORLD.gather(dds)
            else:
                self.devices = [[self.mpi_rank, {'nodes_rank': 0, 'host_name': 'host', 'devices': get_cuda_devices()}]]

            self.devices_rank = {}
            if self.master:
                for d in self.devices:
                    if len(d) > 1:
                        self.devices_rank[d[0]] = d[1]
                    else:
                        self.devices_rank[d[0]] = None

            self.devices_rank = self.MPI.COMM_WORLD.bcast(self.devices_rank)
            self.sync_mpi()
            print(f"{self.mpi_rank = }, {self.devices_rank = }, {all(self.devices_rank.values()) = }")
            # if self.mpi_rank == 7:
            # print(f"{self.devices_rank = }")


    def check_for_bad_input(self, INPUT=None):
        """ Checks for bad user input """
        if INPUT is None:
            INPUT = self.INPUT
        if not self.master:
            return
        # Check deprecated variables
        # check force fields

        logging.info(f'Checking {self.FILES.input_file} input file...')


        if INPUT['gb']['igb'] not in [1, 2, 5, 7, 8]:
            xBFreEErrorLogging('Invalid value for IGB (%s)! ' % INPUT['gb']['igb'] + 'IGB must be 1, 2, 5, 7, or 8.', InputError)
        if INPUT['gb']['intdiel'] < 0:
            xBFreEErrorLogging('INDI must be non-negative!', InputError)
        if INPUT['gb']['extdiel'] < 0:
            xBFreEErrorLogging('EXDI must be non-negative!', InputError)
        if INPUT['gb']['saltcon'] < 0:
            xBFreEErrorLogging('SALTCON must be non-negative!', InputError)
        if INPUT['gb']['surften'] < 0:
            xBFreEErrorLogging('SURFTEN must be non-negative!', InputError)
        if INPUT['gb']['alpb'] == 1 and INPUT['gb']['igb'] == 8:
            xBFreEErrorLogging('IGB=8 is incompatible with ALPB=1! IGB must be 1, 2, 5, or 7 if ALPB=1.', InputError)
        if INPUT['gb']['arad_method'] not in [1, 2, 3]:
            xBFreEErrorLogging('ARAD_METHOD must be 1, 2, or 3!', InputError)
        if INPUT['pb']['indi'] < 0:
            xBFreEErrorLogging('INDI must be non-negative!', InputError)
        if INPUT['pb']['exdi'] < 0:
            xBFreEErrorLogging('EXDI must be non-negative!', InputError)
        if INPUT['pb']['memopt'] > 0 and (INPUT['pb']['emem'] < INPUT['pb']['indi'] or
                                          INPUT['pb']['emem'] > INPUT['pb']['exdi']):
            logging.warning(
                "Membrane dielectric constant (emem) should be between indi and exdi's or there may be errors."
            )
        if INPUT['pb']['scale'] < 0:
            xBFreEErrorLogging('SCALE must be non-negative!', InputError)
        if INPUT['pb']['linit'] < 0:
            xBFreEErrorLogging('LINIT must be a positive integer!', InputError)
        if INPUT['pb']['prbrad'] not in [1.4, 1.6]:
            xBFreEErrorLogging('PRBRAD (%s) must be 1.4 and 1.6!' % INPUT['prbrad'], InputError)
        if INPUT['pb']['istrng'] < 0:
            xBFreEErrorLogging('ISTRNG must be non-negative!', InputError)
        if INPUT['pb']['inp'] not in [1, 2]:
            xBFreEErrorLogging('INP/NPOPT (%s) must be 1, or 2!' % INPUT['inp'], InputError)
        if INPUT['pb']['cavity_surften'] < 0:
            xBFreEErrorLogging('CAVITY_SURFTEN must be non-negative!', InputError)
        if INPUT['pb']['fillratio'] <= 0:
            xBFreEErrorLogging('FILL_RATIO must be positive!', InputError)
        if INPUT['pb']['radiopt'] not in [0, 1]:
            xBFreEErrorLogging('RADIOPT (%s) must be 0 or 1!' % INPUT['radiopt'], InputError)
        if INPUT['nmode']['dielc'] <= 0:
            xBFreEErrorLogging('DIELC must be positive!', InputError)
        if INPUT['nmode']['maxcyc'] < 1:
            xBFreEErrorLogging('MAXCYC must be a positive integer!', InputError)
        if INPUT['pb']['sander_apbs'] not in [0, 1]:
            xBFreEErrorLogging('SANDER_APBS must be 0 or 1!', InputError)
        if INPUT['decomp']['idecomp'] not in [0, 1, 2, 3, 4]:
            xBFreEErrorLogging('IDECOMP (%s) must be 1, 2, 3, or 4!' % INPUT['decomp']['idecomp'], InputError)
        if INPUT['decomp']['idecomp'] != 0 and INPUT['pb']['sander_apbs'] == 1:
            xBFreEErrorLogging('IDECOMP cannot be used with sander.APBS!', InputError)
        if INPUT['decomp']['decomprun'] and INPUT['decomp']['idecomp'] == 0:
            xBFreEErrorLogging('IDECOMP cannot be 0 for Decomposition analysis!', InputError)

        if INPUT['ala']['alarun'] and INPUT['general']['netcdf'] != '':
            xBFreEErrorLogging('Alanine scanning is incompatible with NETCDF != 0!', InputError)
        # if INPUT['general']['PBRadii'] not in range(1, 8):
        #     print(INPUT['general']['PBRadii'])
        #     xBFreEErrorLogging('PBRadii must be 1, 2, 3, 4, 5, 6, or 7!', InputError)
        if INPUT['general']['solvated_trajectory'] not in [0, 1]:
            xBFreEErrorLogging('SOLVATED_TRAJECTORY must be 0 or 1!', InputError)
        if INPUT['gb']['ifqnt'] not in [0, 1]:
            xBFreEErrorLogging('QMMM must be 0 or 1!', InputError)
        if INPUT['gb']['ifqnt'] == 0 and (INPUT['gb']['qm_theory'] or INPUT['gb']['qm_residues']):
            logging.warning('qm_theory/qm_residues variable has been defined, however the potential function is '
                            'strictly classical (ifqnt=0). Please, set ifqnt=1 if you want to use Use QM/MM')
        if (
                not INPUT['gb']['molsurf']
                and (INPUT['gb']['msoffset'] != 0 or INPUT['gb']['probe'] != 1.4)
                and self.master
        ):
            logging.warning('offset and probe are molsurf-only options')
        if INPUT['gb']['ifqnt'] == 1:
            if INPUT['gb']['qm_theory'] not in ['PM3', 'AM1', 'MNDO', 'PDDG-PM3', 'PM3PDDG',
                                                'PDDG-MNDO', 'PDDGMNDO', 'PM3-CARB1',
                                                'PM3CARB1', 'DFTB', 'SCC-DFTB', 'RM1', 'PM6',
                                                'PM3-ZnB', 'PM3-MAIS', 'PM6-D', 'PM6-DH+',
                                                'AM1-DH+', 'AM1-D*', 'PM3ZNB', 'MNDO/D',
                                                'MNDOD']:
                xBFreEErrorLogging('Invalid QM_THEORY (%s)! ' % INPUT['gb']['qm_theory'] +
                                'This variable must be set to allowable options.\n' +
                                'PM3, AM1, MNDO, PDDG-PM3, PM3PDDG, PDDG-MNDO, PDDGMNDO, \n'
                                'PM3-CARB1, PM3CARB1, DFTB, SCC-DFTB, RM1, PM6, PM3-ZnB, \n'
                                'PM3-MAIS, PM6-D, PM6-DH+, AM1-DH+, AM1-D*, PM3ZNB, MNDO/D, MNDOD', InputError)
            if INPUT['gb']['qm_residues'] == '':
                xBFreEErrorLogging('QM_RESIDUES must be specified for IFQNT = 1!', InputError)
            if INPUT['decomp']['decomprun']:
                xBFreEErrorLogging('QM/MM and decomposition are incompatible!', InputError)
            if (INPUT['gb']['qmcharge_lig'] + INPUT['gb']['qmcharge_rec'] !=
                    INPUT['gb']['qmcharge_com'] and not self.stability):
                xBFreEErrorLogging('The total charge of the ligand and receptor ' +
                                'does not equal the charge of the complex!', InputError)
            if INPUT['gb']['scfconv'] < 1.0e-12:
                logging.warning('There is a risk of convergence problems when the requested convergence is less than '
                                '1.0e-12 kcal/mol')
            if INPUT['gb']['writepdb']:
                logging.info('Writing qmmm_region.pdb PDB file of the selected QM region...')
            if INPUT['gb']['verbosity'] not in [0, 1, 2, 3, 4, 5]:
                xBFreEErrorLogging('VERBOSITY must be 0, 1, 2, 3, 4 or 5!', InputError)
            if INPUT['general']['verbose'] >= 2:
                logging.warning('VERBOSITY values of 2 or higher will produce a lot of output')

        if INPUT['rism']['rismrun']:
            if INPUT['rism']['rism_verbose'] not in [0, 1, 2]:
                xBFreEErrorLogging('RISM_VERBOSE must be 0, 1, or 2!', InputError)
            if INPUT['rism']['buffer'] < 0 and INPUT['rism']['solvcut'] < 0:
                xBFreEErrorLogging('If BUFFER < 0, SOLVCUT must be > 0!', InputError)
            for tol in INPUT['rism']['tolerance']:
                if tol <= 0:
                    xBFreEErrorLogging('TOLERANCE must be positive!', InputError)
            if INPUT['rism']['tolerance'][-1] > 0.00001:
                logging.warning(
                    f"Default TOLERANCE value is 0.00001! However {INPUT['rism']['tolerance'][-1]} is been used. "
                    f"Check documentation for more details...")
            if INPUT['rism']['buffer'] < 0 and INPUT['rism']['ng'] == '':
                xBFreEErrorLogging('You must specify NG if BUFFER < 0!', InputError)
            if INPUT['rism']['polardecomp'] not in [0, 1]:
                xBFreEErrorLogging('POLARDECOMP must be either 0 or 1!', InputError)
            if INPUT['rism']['entropicdecomp'] not in [0, 1]:
                xBFreEErrorLogging('ENTROPICDECOMP must be either 0 or 1!', InputError)
            for i in zip(['treeDCF', 'treeTCF', 'treeCoulomb'], [INPUT['rism']['treeDCF'], INPUT['rism']['treeTCF'],
                                                                 INPUT['rism']['treeCoulomb']]):
                if i[1] not in [0, 1]:
                    xBFreEErrorLogging(f'{i[0]} must be either 0 or 1!', InputError)
            # if INPUT['thermo'] not in ['std', 'gf', 'both']:
            #     xBFreEErrorLogging('THERMO must be "std", "gf", "both"!', InputError)
            # TODO: include other corrections? pc+?
            # if INPUT['thermo'] not in ['std', 'gf', 'pc+', 'all']:
            #     xBFreEErrorLogging('THERMO must be "std", "gf", "pc+" or "all"!', InputError)
        if (
                not INPUT['gb']['gbrun']
                and not INPUT['pb']['pbrun']
                and not INPUT['rism']['rismrun']
                and not INPUT['nmode']['nmoderun']
                and not INPUT['general']['qh_entropy']
                and not INPUT['gbnsr6']['gbnsr6run']
        ):
            xBFreEErrorLogging('You did not specify any type of calculation!', InputError)

        if INPUT['gb']['gbrun'] and INPUT['general']['PBRadii'] == 7:
            xBFreEErrorLogging('PBRadii = 7 (charmm_radii) is compatible only with &pb!', InputError)

        if INPUT['decomp']['decomprun'] and \
                not (INPUT['gb']['gbrun'] or INPUT['pb']['pbrun'] or INPUT['gbnsr6']['gbnsr6run']):
            xBFreEErrorLogging('DECOMP must be run with either GB or PB!', InputError)

        if '-deo' in sys.argv and not INPUT['decomp']['decomprun']:
            logging.warning("&decomp namelist has not been defined in the input file. Ignoring '-deo' flag... ")

        if (
                not INPUT['gb']['molsurf']
                and (INPUT['gb']['msoffset'] != 0 or INPUT['gb']['probe'] != 1.4)
                and self.master
        ):
            logging.warning('offset and probe are molsurf-only options')
        if INPUT['ala']['cas_intdiel'] not in [0, 1]:
            xBFreEErrorLogging('cas_intdiel must be set to 0 or 1!', InputError)
        # User warning when intdiel > 10
        if self.INPUT['gb']['intdiel'] > 10:
            logging.warning('Intdiel is greater than 10...')
        # check mutant definition
        if self.INPUT['ala']['mutant'].upper() not in ['ALA', 'A', 'GLY', 'G']:
            xBFreEErrorLogging('The mutant most be ALA (or A) or GLY (or G)', InputError)

        # fixed the error when try to open gmx_MMPBSA_ana in the issue
        # https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/33
        if self.INPUT['general']['startframe'] < 1:
            # xBFreEErrorLogging('The startframe variable must be >= 1')
            logging.warning(f"The startframe variable must be >= 1. Changing startframe from"
                            f" {self.INPUT['general']['startframe']} to 1")
            self.INPUT['general']['startframe'] = 1
        if INPUT['nmode']['nmoderun']:
            if self.INPUT['nmode']['nmstartframe'] < 1:
                logging.warning(f"The nmstartframe variable must be >= 1. Changing nmstartframe from"
                                f" {self.INPUT['nmode']['nmstartframe']} to 1")
                self.INPUT['nmode']['nmstartframe'] = 1
            if INPUT['nmode']['drms'] > 0.001:
                logging.warning(
                    f"Default DRMS value is 0.001! However {INPUT['nmode']['drms']} is being used. Check "
                    f'documentation for more details...')
            if INPUT['nmode']['maxcyc'] < 10000:
                logging.warning(f"Default MAXCYC value is 10000! However {INPUT['nmode']['maxcyc']} is being used. "
                                f"Check documentation for more details...")

        # Assigning variables
        # set the pbtemp = temperature
        self.INPUT['pb']['pbtemp'] = self.INPUT['general']['temperature']
        # self.INPUT['gbnsr6']['istrng'] = self.INPUT['gbnsr6']['istrng'] * 1000

        logging.info(f'Checking {self.FILES.input_file} input file... Done.\n')

    def remove(self, removetemp=0):
        """ Removes temporary files """
        if not self.master:
            return
        misc.remove('mmpbsa', removetemp)

    def sync_mpi(self):
        """ Throws up a barrier """
        self.MPI.COMM_WORLD.Barrier()

    def parse_output_files(self, from_calc=True):
        """
        This parses the output files and loads them into dicts for easy access
        """
        # Only the master does this
        from types import SimpleNamespace
        if not self.master:
            return
        logging.info('Parsing results to output files...\n')
        self.calc_types = SimpleNamespace(normal={}, mutant={}, mut_norm={}, decomp_normal={}, decomp_mutant={})
        INPUT, FILES = self.INPUT, self.FILES
        # Quasi-harmonic analysis is a special-case, so handle that separately
        if INPUT['general']['qh_entropy']:
            if not INPUT['ala']['mutant_only']:
                self.calc_types.normal['qh'] = QHout('cpptraj_entropy.out', INPUT['general']['temperature'])
            if INPUT['alarun']:
                self.calc_types.mutant['qh'] = QHout('mutant_cpptraj_entropy.out', INPUT['general']['temperature'])
            if INPUT['ala']['alarun'] and not INPUT['ala']['mutant_only']:
                self.calc_types.mut_norm['qh'] = DeltaDeltaQH(self.calc_types.mutant['qh'],
                                                              self.calc_types.normal['qh'])

        # Determine if our GB is QM/MM or not
        GBClass = QMMMout if INPUT['gb']['ifqnt'] else GBout
        # Determine which kind of RISM output class we are based on std/gf and
        # polardecomp
        if INPUT['rism']['polardecomp']:
            RISM_GF = PolarRISM_gf_Out
            RISM_Std = PolarRISM_std_Out
            RISM_PCplus = PolarRISM_pcplus_Out
        else:
            RISM_GF = RISM_gf_Out
            RISM_Std = RISM_std_Out
            RISM_PCplus = RISM_pcplus_Out
        # Now we make a list of the other calculation types, their INPUT triggers,
        # their key in the calc_types dict, the base name of their output files
        # without the prefix (with %s-substitution for complex, receptor, or
        # ligand), and the class for their output
        nmls = ('nmode', 'gb', 'pb', 'rism', 'rism', 'rism', 'gbnsr6')
        triggers = ('nmoderun', 'gbrun', 'pbrun', 'rismrun_std', 'rismrun_gf', 'rismrun_pcplus', 'gbnsr6run')
        outclass = (NMODEout, GBClass, PBout, RISM_Std, RISM_GF, RISM_PCplus, GBNSR6out)
        outkey = ('nmode', 'gb', 'pb', 'rism std', 'rism gf', 'rism pcplus', 'gbnsr6')
        basename = ('%s_nm.out', '%s_gb.mdout', '%s_pb.mdout', '%s_rism.mdout', '%s_rism.mdout', '%s_rism.mdout',
                    '%s_gbnsr6.mdout')

        for i, key in enumerate(outkey):
            if not INPUT.get(nmls[i]) or not INPUT[nmls[i]].get(triggers[i]) or not INPUT[nmls[i]][triggers[i]]:
                continue
            numframes = self.numframes_nmode if key == 'nmode' else self.numframes
            # Non-mutant
            if not INPUT['ala']['mutant_only']:
                self.calc_types.normal[key] = {'complex': outclass[i]('complex', self.INPUT, self.using_chamber)}
                self.calc_types.normal[key]['complex'].parse_from_file(basename[i] % 'complex',
                                                                       self.mpi_size, numframes)
                # check if the nmode output is valid
                if self.calc_types.normal[key]['complex'].no_nmode_convergence:
                    self.INPUT['nmode']['nmoderun'] = False
                    del self.calc_types.normal[key]
                    continue

                if not self.stability:
                    self.calc_types.normal[key]['receptor'] = outclass[i]('receptor', self.INPUT, self.using_chamber)
                    self.calc_types.normal[key]['receptor'].parse_from_file(basename[i] % 'receptor',
                                                                            self.mpi_size, numframes)
                    self.calc_types.normal[key]['ligand'] = outclass[i]('ligand', self.INPUT, self.using_chamber)
                    self.calc_types.normal[key]['ligand'].parse_from_file(basename[i] % 'ligand',
                                                                          self.mpi_size, numframes)
                    self.calc_types.normal[key]['delta'] = BindingStatistics(self.calc_types.normal[key]['complex'],
                                                                             self.calc_types.normal[key]['receptor'],
                                                                             self.calc_types.normal[key]['ligand'],
                                                                             self.using_chamber, self.traj_protocol)
            # Time for mutant
            if INPUT['ala']['alarun']:
                self.calc_types.mutant[key] = {'complex': outclass[i]('Mutant-Complex', self.INPUT, self.using_chamber)}
                self.calc_types.mutant[key]['complex'].parse_from_file('mutant_' + basename[i] % 'complex',
                                                                       self.mpi_size, numframes)
                if not self.stability:
                    self.calc_types.mutant[key]['receptor'] = outclass[i]('Mutant-Receptor', self.INPUT,
                                                                          self.using_chamber)
                    self.calc_types.mutant[key]['receptor'].parse_from_file('mutant_' + basename[i] %
                                                                            'receptor', self.mpi_size, numframes)
                    self.calc_types.mutant[key]['ligand'] = outclass[i]('Mutant-Ligand', self.INPUT,
                                                                        self.using_chamber)
                    self.calc_types.mutant[key]['ligand'].parse_from_file('mutant_' + basename[i] % 'ligand',
                                                                          self.mpi_size, numframes)
                    self.calc_types.mutant[key]['delta'] = BindingStatistics(self.calc_types.mutant[key]['complex'],
                                                                             self.calc_types.mutant[key]['receptor'],
                                                                             self.calc_types.mutant[key]['ligand'],
                                                                             self.using_chamber, self.traj_protocol)

            if INPUT['ala']['alarun'] and not INPUT['ala']['mutant_only']:
                self.calc_types.mut_norm[key] = {'complex': DeltaDeltaStatistics(
                    self.calc_types.mutant[key]['complex'], self.calc_types.normal[key]['complex'])}
                if not self.stability:
                    if self.FILES.receptor_prmtop != self.FILES.mutant_receptor_prmtop:
                        self.calc_types.mut_norm[key]['receptor'] = DeltaDeltaStatistics(
                            self.calc_types.mutant[key]['receptor'], self.calc_types.normal[key]['receptor'])
                    else:
                        self.calc_types.mut_norm[key]['ligand'] = DeltaDeltaStatistics(
                            self.calc_types.mutant[key]['ligand'], self.calc_types.normal[key]['ligand'])
                    self.calc_types.mut_norm[key]['delta'] = DeltaDeltaStatistics(
                        self.calc_types.mutant[key]['delta'], self.calc_types.normal[key]['delta'])

            self.get_iec2entropy(from_calc)

        if not hasattr(self, 'resl'):
            from xBFreE.utils.molecule import mask2list
            self.resl = mask2list(FILES.complex_fixed, INPUT['general']['receptor_mask'], INPUT['general']['ligand_mask'])
            if INPUT['ala']['alarun']:
                self.resl[self.mutant_index].set_mut(INPUT['ala']['mutant'])

        if INPUT['decomp']['decomprun']:
            self._get_decomp()

    def get_iec2entropy(self, from_calc):
        allowed_met = ['gb', 'pb', 'rism std', 'rism gf', 'rism pcplus', 'gbnsr6']
        if self.INPUT['general']['interaction_entropy']:
            self.calc_types.normal['ie'] = {}
        if self.INPUT['general']['c2_entropy']:
            self.calc_types.normal['c2'] = {}

        for key in allowed_met:
            if self.INPUT['general']['interaction_entropy']:
                if key in self.calc_types.normal:
                    if from_calc:
                        edata = self.calc_types.normal[key]['delta']['GGAS']
                        logging.info('Beginning Interaction Entropy calculations...')
                        ie = InteractionEntropyCalc(edata, self.INPUT, key)
                        ie.save_output(f'normal_{key}_IE.dat')

                    self.calc_types.normal['ie'][key] = IEout(self.INPUT, key)
                    self.calc_types.normal['ie'][key].parse_from_file(f'normal_{key}_IE.dat',
                                                                        self.numframes)
                if key in self.calc_types.mutant:
                    if from_calc:
                        edata = self.calc_types.mutant[key]['delta']['GGAS']
                        logging.info('Beginning Mutant Interaction Entropy calculations...')
                        mie = InteractionEntropyCalc(edata, self.INPUT, key)
                        mie.save_output(f'mutant_{key}_IE.dat')

                    self.calc_types.mutant['ie'][key] = IEout(self.INPUT, key)
                    self.calc_types.mutant['ie'][key].parse_from_file(f'mutant_{key}_IE.dat',
                                                                 self.numframes)

                if self.INPUT['ala']['alarun'] and not self.INPUT['ala']['mutant_only'] and key in self.calc_types.normal:
                    self.calc_types.mut_norm['ie'][key] = DeltaIEC2Statistic(
                        self.calc_types.mutant['ie'][key], self.calc_types.normal['ie'][key])

            if self.INPUT['general']['c2_entropy']:
                if not self.INPUT['ala']['mutant_only'] and key in self.calc_types.normal:
                    if from_calc:
                        edata = self.calc_types.normal[key]['delta']['GGAS']
                        logging.info('Beginning C2 Entropy calculations...')
                        c2 = C2EntropyCalc(edata, self.INPUT, key)
                        c2.save_output('normal_{key}_c2_entropy.dat')

                    self.calc_types.normal['c2'][key] = C2out(key)
                    self.calc_types.normal['c2'][key].parse_from_file(f'normal_{key}_c2_entropy.dat')
                if key in self.calc_types.mutant:
                    if from_calc:
                        edata = self.calc_types.mutant[key]['delta']['GGAS']
                        logging.info('Beginning Mutant C2 Entropy calculations...')
                        c2 = C2EntropyCalc(edata, self.INPUT, key)
                        c2.save_output(f'mutant_{key}_c2_entropy.dat')

                    self.calc_types.mutant['c2'][key] = C2out(key)
                    self.calc_types.mutant['c2'][key].parse_from_file(f'mutant_{key}_c2_entropy.dat')

                if self.INPUT['ala']['alarun'] and not self.INPUT['ala']['mutant_only'] and key in self.calc_types.normal:
                    self.calc_types.mut_norm['c2'][key] = DeltaIEC2Statistic(
                        self.calc_types.mutant['c2'][key], self.calc_types.normal['c2'][key])

    def _get_mm_data(self):

        # # FIXME: add pbsa.cuda, APBS and PBDelphi
        # if not self.INPUT['gbnsr6']['gbnsr6run']:
        #     return

        mm_data = {}

        if not self.INPUT['ala']['mutant_only']:
            mm_data['normal'] = {'complex': MMout('complex', self.INPUT, self.using_chamber)}
            mm_data['normal']['complex'].parse_from_file("complex_mm.mdout", self.mpi_size, self.numframes)
            if not self.stability:
                mm_data['normal']['receptor'] = MMout('receptor', self.INPUT, self.using_chamber)
                mm_data['normal']['receptor'].parse_from_file("receptor_mm.mdout", self.mpi_size,
                                                              self.numframes)
                mm_data['normal']['ligand'] = MMout('ligand', self.INPUT, self.using_chamber)
                mm_data['normal']['ligand'].parse_from_file("ligand_mm.mdout", self.mpi_size, self.numframes)

        # Time for mutant
        if self.INPUT['ala']['alarun']:
            mm_data['mutant'] = {'complex': MMout('mutant-complex', self.INPUT, self.using_chamber)}
            mm_data['mutant']['complex'].parse_from_file("mutant_complex_mm.mdout", self.mpi_size,
                                                         self.numframes)
            if not self.stability:
                mm_data['mutant']['receptor'] = MMout('mutant-receptor', self.INPUT, self.using_chamber)
                mm_data['mutant']['receptor'].parse_from_file("mutant_receptor_mm.mdout", self.mpi_size,
                                                              self.numframes)
                mm_data['mutant']['ligand'] = MMout('mutant-ligand', self.INPUT, self.using_chamber)
                mm_data['mutant']['ligand'].parse_from_file("mutant_ligand_mm.mdout", self.mpi_size,
                                                            self.numframes)
        return mm_data

    def _res2print(self):
        """
        Get residues list from print_res variable
        Returns: residues to print list
        """
        print_res = []
        for x in self.INPUT['decomp']['print_res'].split(','):
            r = list(map(int, x.split('-')))
            s = r[0]
            e = r[-1] + 1
            print_res.extend(range(s, e))
        return print_res

    def _get_decomp(self):
        from xBFreE.mmpbsa.output.amber import (DecompOut, PairDecompOut, DecompBinding, PairDecompBinding)
        nmls = ('gb', 'pb', 'gbnsr6')
        outkey = ('gb', 'pb', 'gbnsr6')
        triggers = ('gbrun', 'pbrun', 'gbnsr6run')
        basename = ('%s_gb.mdout', '%s_pb.mdout', '%s_gbnsr6.mdout')
        INPUT, FILES = self.INPUT, self.FILES
        headers = {'gb': 'Generalized Born', 'pb': 'Poisson Boltzmann', 'gbnsr6': 'Generalized Born (R6)'}
        if INPUT['decomp']['idecomp'] in [1, 2]:
            DecompBindingClass = DecompBinding
            DecompClass = DecompOut
        # Pairwise
        else:
            DecompBindingClass = PairDecompBinding
            DecompClass = PairDecompOut

        # get residues list from print_res variable
        print_res = self._res2print()
        com_list = {}
        rec_list = {}
        lig_list = {}
        for x in self.resl:
            if x.index in print_res:
                com_list[x.index] = x
                if x.is_receptor():
                    rec_list[x.id_index] = x
                else:
                    lig_list[x.id_index] = x

        for i, key in enumerate(outkey):
            if not INPUT.get(nmls[i]) or not INPUT[nmls[i]].get(triggers[i]) or not INPUT[nmls[i]][triggers[i]]:
                continue
            # FIXME
            if key == 'gb':
                surften = INPUT['gb']['surften']
            elif key == 'pb':
                surften = INPUT['pb']['cavity_surften']
            else: # gbnsr6
                surften = INPUT['gbnsr6']['cavity_surften']

            if not self.INPUT['ala']['mutant_only']:
                self.calc_types.decomp_normal[key] = {'complex': DecompClass('complex')}
                self.calc_types.decomp_normal[key]['complex'].parse_from_file(basename[i] % 'complex',
                                                                              com_list, INPUT, surften,
                                                                              self.mpi_size, self.numframes)
                if not self.stability:
                    self.calc_types.decomp_normal[key]['receptor'] = DecompClass('receptor')
                    self.calc_types.decomp_normal[key]['receptor'].parse_from_file(basename[i] % 'receptor',
                                                                                   rec_list, INPUT, surften,
                                                                                   self.mpi_size, self.numframes)
                    self.calc_types.decomp_normal[key]['ligand'] = DecompClass('ligand')
                    self.calc_types.decomp_normal[key]['ligand'].parse_from_file(basename[i] % 'ligand',
                                                                                 lig_list, INPUT, surften,
                                                                                 self.mpi_size, self.numframes)
                    self.calc_types.decomp_normal[key]['delta'] = DecompBindingClass(
                        self.calc_types.decomp_normal[key]['complex'], self.calc_types.decomp_normal[key]['receptor'],
                        self.calc_types.decomp_normal[key]['ligand'], INPUT,
                        f'Energy Decomposition Analysis (All units kcal/mol): {headers[key]} model')

            if INPUT['ala']['alarun']:
                # Do mutant
                self.calc_types.decomp_mutant[key] = {'complex': DecompClass('Mutant-Complex')}
                self.calc_types.decomp_mutant[key]['complex'].parse_from_file(
                    ('mutant_' + basename[i] % 'complex'),
                    com_list,
                    INPUT,
                    surften,
                    self.mpi_size,
                    self.numframes,
                    True
                )

                if not self.stability:
                    self.calc_types.decomp_mutant[key]['receptor'] = DecompClass('Mutant-Receptor')
                    self.calc_types.decomp_mutant[key]['receptor'].parse_from_file(
                        ('mutant_' + basename[i] % 'receptor'),
                        rec_list,
                        INPUT,
                        surften,
                        self.mpi_size,
                        self.numframes,
                        True
                    )

                    self.calc_types.decomp_mutant[key]['ligand'] = DecompClass('Mutant-Ligand')
                    self.calc_types.decomp_mutant[key]['ligand'].parse_from_file(
                        ('mutant_' + basename[i] % 'ligand'),
                        lig_list,
                        INPUT,
                        surften,
                        self.mpi_size,
                        self.numframes,
                        True
                    )

                    self.calc_types.decomp_mutant[key]['delta'] = DecompBindingClass(
                        self.calc_types.decomp_mutant[key]['complex'], self.calc_types.decomp_mutant[key]['receptor'],
                        self.calc_types.decomp_mutant[key]['ligand'], INPUT,
                        f'Energy Decomposition Analysis (All units kcal/mol): {headers[key]} model ({self.mut_str})')


# Local methods

def excepthook(exception_type, exception_value, tb):
    """
    Replaces sys.excepthook so fatal exceptions kill all MPI threads and we can
    control the printing of tracebacks. Those are helpful for debugging purposes,
    but may be unsightly to users. debug_printlevel set above controls this
    behavior
    """
    import traceback
    global _stderr, _mpi_size, _rank
    if not isinstance(exception_type, xBFreE_Error):
        traceback.print_tb(tb)
    _stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
    if _mpi_size > 1:
        _stderr.write('Error occurred on rank %d.' % _rank + os.linesep)
    _stderr.write('Exiting. All files have been retained.' + os.linesep)
    _MPI.COMM_WORLD.Abort(1)


def interrupt_handler(signal, frame):
    """ Handles interrupt signals for a clean exit """
    global _MPI, _stderr
    _stderr.write('\n%s interrupted! Program terminated. All files are kept.\n' %
                  os.path.split(sys.argv[0])[1])
    _MPI.COMM_WORLD.Abort(1)

