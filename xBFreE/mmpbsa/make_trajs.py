"""
This module contains classes and such that are responsible for generating
trajectories for MM/PBSA based on input, starting trajectories. It will
strip trajectories and output the ones that will be used for the calculations.
Needed for proper operation of xBFreE

Methods:
         make_trajectories(INPUT, FILES, size): Makes all the non-mutant
            trajectories needed for the calculation, as well as all dummy
            files (i.e. restarts and PDBs)

         make_mutant_trajectories(INPUT, FILES, rank): Mutates the trajectories

Classes:
         Trajectory: Class for manipulating Amber trajectories through cpptraj
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

from warnings import warn
from xBFreE.exceptions import (TrajError, xBFreE_Error, InternalError, MutantResError)
from pathlib import Path

strip_mask = ':WAT,Cl*,CIO,Cs+,IB,K*,Li+,MG*,Na+,Rb+,CS,RB,NA,F,CL'


def make_trajectories(INPUT, FILES, cpptraj, size, devices=None):
    """
    This function creates the necessary trajectory files, and creates thread-specific trajectories for parallel
    calculations
    """
    stability = FILES.stability

    # File suffix is dependent on file type
    if INPUT['general']['netcdf']:
        trj_suffix = 'nc'
    else:
        trj_suffix = 'mdcrd'

    traj = Trajectory(FILES.complex_prmtop, FILES.complex_trajs, cpptraj)
    traj.Setup(INPUT['general']['startframe'], INPUT['general']['endframe'], INPUT['general']['interval'])
    # RMS fit
    traj.rms('!(%s)' % strip_mask)

    com_frames = int(traj.processed_frames)
    rec_frames = 0
    lig_frames = 0
    num_frames_nmode = 0

    # Sanity check
    if traj.processed_frames < size:
        raise xBFreE_Error('Must have at least as many frames as processors!')

    # We now know how many frames we have in total, so make a list that lists the
    # number of frames found for each rank, and assign extra frames incrementally
    frames_per_rank = traj.processed_frames // size
    extras = traj.processed_frames - frames_per_rank * size
    frame_count = [frames_per_rank for i in range(size)]
    for i in range(size):
        if i < extras: frame_count[i] += 1

    # Dump our complex trajectories
    if INPUT['general']['full_traj'] or INPUT['general']['qh_entropy']:
        traj.Outtraj('complex.%s' % trj_suffix, filetype=INPUT['general']['netcdf'])
    traj.Outtraj('complex.pdb', frames='1', filetype='pdb')
    traj.Outtraj('dummycomplex.inpcrd', frames='1', filetype='restart')

    # Now dump thread-specific trajectories
    last_frame = 1
    for i in range(size):
        frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
        traj.Outtraj('complex.%s.%d' % (trj_suffix, i), frames=frame_string, filetype=INPUT['general']['netcdf'])
        # TODO: include pbsa.cuda. For APBS and pbsa.cuda we need to generate pqr instead
        if INPUT['gbnsr6']['gbnsr6run']:
            temp_dir = Path(f"inpcrd_{i}")
            temp_dir.mkdir()
            traj.Outtraj(f"inpcrd_{i}/complex.inpcrd", frames=frame_string, filetype='restart', options=['keepext'])
        last_frame += frame_count[i]

    if devices:
        # FIXME: we make this here because there is not support for md in gpu processes
        devices_size = sum([1 for d in devices.values() if d is not None])
        frames_per_rank_gpu = traj.processed_frames // devices_size
        extras_gpu = traj.processed_frames - frames_per_rank_gpu * devices_size
        frame_count_gpu = {r: frames_per_rank_gpu for r, d in devices.items() if d}
        for i in range(devices_size):
            if i < extras_gpu: frame_count_gpu[i] += 1

        last_frame_gpu = 1
        for r, d in devices.items():
            if d:
                frame_string = '%d-%d' % (last_frame_gpu, last_frame_gpu + frame_count_gpu[r] - 1)
                # TODO: only for pbsa.cuda?
                temp_dir = Path(f"inpcrd_gpu_{r}")
                temp_dir.mkdir()
                traj.Outtraj(f"inpcrd_gpu_{r}/complex.inpcrd", frames=frame_string, filetype='restart',
                             options=['keepext'])
                last_frame_gpu += frame_count_gpu[r]

    # Now create the receptor/ligand trajectories if we're taking them from
    # the complex trajectory

    if not stability and not FILES.receptor_trajs:
        traj.Strip(INPUT['general']['ligand_mask'])
        if INPUT['general']['full_traj'] or INPUT['general']['qh_entropy']:
            traj.Outtraj('receptor.%s' % trj_suffix, filetype=INPUT['general']['netcdf'])
        traj.Outtraj('receptor.pdb', frames='1', filetype='pdb')
        traj.Outtraj('dummyreceptor.inpcrd', frames='1', filetype='restart')
        last_frame = 1
        for i in range(size):
            frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
            traj.Outtraj('receptor.%s.%d' % (trj_suffix, i), frames=frame_string,
                         filetype=INPUT['general']['netcdf'])
            # FIXME: include pbsa.cuda. For APBS and pbsa.cuda we need to generate pqr instead
            if INPUT['gbnsr6']['gbnsr6run']:
                traj.Outtraj(f"inpcrd_{i}/receptor.inpcrd", frames=frame_string, filetype='restart',
                             options=['keepext'])
            last_frame += frame_count[i]

        if devices:
            last_frame_gpu = 1
            for r, d in devices.items():
                if d:
                    frame_string = '%d-%d' % (last_frame_gpu, last_frame_gpu + frame_count_gpu[r] - 1)
                    # TODO: only for pbsa.cuda?
                    traj.Outtraj(f"inpcrd_gpu_{r}/receptor.inpcrd", frames=frame_string, filetype='restart',
                                 options=['keepext'])
                    last_frame_gpu += frame_count_gpu[r]

        traj.Unstrip(restrip_solvent=True)
        traj.rms('!(%s)' % strip_mask)

    if not stability and not FILES.ligand_trajs:
        traj.Strip(INPUT['general']['receptor_mask'])
        if INPUT['general']['full_traj'] or INPUT['general']['qh_entropy']:
            traj.Outtraj('ligand.%s' % trj_suffix, filetype=INPUT['general']['netcdf'])
        traj.Outtraj('ligand.pdb', frames='1', filetype='pdb')
        traj.Outtraj('dummyligand.inpcrd', frames='1', filetype='restart')
        last_frame = 1
        for i in range(size):
            frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
            traj.Outtraj('ligand.%s.%d' % (trj_suffix, i), frames=frame_string,
                         filetype=INPUT['general']['netcdf'])
            # FIXME: include pbsa.cuda. For APBS and pbsa.cuda we need to generate pqr instead
            if INPUT['gbnsr6']['gbnsr6run']:
                traj.Outtraj(f"inpcrd_{i}/ligand.inpcrd", frames=frame_string, filetype='restart',
                             options=['keepext'])
            last_frame += frame_count[i]

        if devices:
            last_frame_gpu = 1
            for r, d in devices.items():
                if d:
                    frame_string = '%d-%d' % (last_frame_gpu, last_frame_gpu + frame_count_gpu[r] - 1)
                    # TODO: only for pbsa.cuda?
                    traj.Outtraj(f"inpcrd_gpu_{r}/ligand.inpcrd", frames=frame_string, filetype='restart',
                                 options=['keepext'])
                    last_frame_gpu += frame_count_gpu[r]

        traj.Unstrip(restrip_solvent=True)
        traj.rms('!(%s)' % strip_mask)

    # Run cpptraj to get the trajectory
    traj.Run('normal_traj_cpptraj.out')


    # FIXME: since the com, rec and lig must have the same number of frames in MT, we should reformat this code
    # Go back and do the receptor and ligand if we used a multiple
    # trajectory approach
    if not stability and FILES.receptor_trajs:
        rectraj = Trajectory(FILES.receptor_prmtop, FILES.receptor_trajs, cpptraj)
        rectraj.Setup(INPUT['general']['startframe'], INPUT['general']['endframe'], INPUT['general']['interval'])
        rec_frames = int(rectraj.processed_frames)
        rectraj.rms('!(%s)' % strip_mask)
        if INPUT['general']['full_traj'] or INPUT['general']['qh_entropy']:
            rectraj.Outtraj('receptor.%s' % trj_suffix, filetype=INPUT['general']['netcdf'])
        rectraj.Outtraj('receptor.pdb', frames='1', filetype='pdb')
        rectraj.Outtraj('dummyreceptor.inpcrd', frames='1',
                        filetype='restart')

        # Now do the same split-up of workload as we did for complex, but don't
        # assume the same number of frames as we had for the complex
        if rectraj.processed_frames < size:
            raise xBFreE_Error('Too many procs for receptor snapshots')
        frames_per_rank = rectraj.processed_frames // size
        extras = rectraj.processed_frames - frames_per_rank * size
        frame_count = [frames_per_rank for i in range(size)]
        for i in range(size):
            if i < extras: frame_count[i] += 1
        last_frame = 1
        for i in range(size):
            frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
            rectraj.Outtraj('receptor.%s.%d' % (trj_suffix, i), frames=frame_string,
                            filetype=INPUT['general']['netcdf'])
            # FIXME: include pbsa.cuda. For APBS and PBDelphi we need to generate pqr instead
            if INPUT['gbnsr6']['gbnsr6run']:
                traj.Outtraj(f"inpcrd_{i}/receptor.inpcrd", frames=frame_string, filetype='restart',
                             options=['keepext'])
            last_frame += frame_count[i]

        rectraj.Run('receptor_traj_cpptraj.out')

    # end if not stability and FILES.receptor_trajs

    if not stability and FILES.ligand_trajs:
        ligtraj = Trajectory(FILES.ligand_prmtop, FILES.ligand_trajs, cpptraj)
        ligtraj.Setup(INPUT['general']['startframe'], INPUT['general']['endframe'], INPUT['general']['interval'])
        lig_frames = int(ligtraj.processed_frames)
        ligtraj.rms('!(%s)' % strip_mask)
        ligtraj.Outtraj('ligand.%s' % trj_suffix, filetype=INPUT['general']['netcdf'])
        ligtraj.Outtraj('ligand.pdb', frames='1', filetype='pdb')
        ligtraj.Outtraj('dummyligand.inpcrd', frames='1', filetype='restart')

        # Now do the same split-up of workload as we did for complex, but don't
        # assume the same number of frames as we had for the complex
        if ligtraj.processed_frames < size:
            raise xBFreE_Error('Too many procs for ligand snapshots')
        frames_per_rank = ligtraj.processed_frames // size
        extras = ligtraj.processed_frames - frames_per_rank * size
        frame_count = [frames_per_rank for i in range(size)]
        for i in range(size):
            if i < extras: frame_count[i] += 1
        last_frame = 1
        for i in range(size):
            frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
            ligtraj.Outtraj('ligand.%s.%d' % (trj_suffix, i), frames=frame_string,
                            filetype=INPUT['general']['netcdf'])
            # FIXME: include pbsa.cuda. For APBS and apbs.cuda we need to generate pqr instead
            if INPUT['gbnsr6']['gbnsr6run']:
                traj.Outtraj(f"inpcrd_{i}/ligand.inpcrd", frames=frame_string, filetype='restart',
                             options=['keepext'])
            last_frame += frame_count[i]

        ligtraj.Run('ligand_traj_cpptraj.out')

    # end if not stability and FILES.ligand_trajs

    # Now make the nmode trajectories
    if INPUT['nmode']['nmoderun']:
        nmtraj = Trajectory(FILES.complex_prmtop, ['complex.%s.%d' %
                                                   (trj_suffix, i) for i in range(size)], cpptraj)
        nmtraj.Setup(INPUT['nmode']['nmstartframe'], INPUT['nmode']['nmendframe'], INPUT['nmode']['nminterval'])

        num_frames_nmode = int(nmtraj.processed_frames)

        # Now split up the complex trajectory by thread
        if nmtraj.processed_frames < size:
            raise xBFreE_Error('More processors than complex nmode frames!')

        frames_per_rank = nmtraj.processed_frames // size
        extras = nmtraj.processed_frames - frames_per_rank * size
        frame_count = [frames_per_rank for i in range(size)]
        for i in range(size):
            if i < extras: frame_count[i] += 1
        last_frame = 1
        for i in range(size):
            frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
            nmtraj.Outtraj('complex_nm.%s.%d' % (trj_suffix, i),
                           frames=frame_string, filetype=INPUT['general']['netcdf'])
            last_frame += frame_count[i]

        nmtraj.Run('com_nm_traj_cpptraj.out')

        if not stability:
            nmtraj = Trajectory(FILES.receptor_prmtop, [f'receptor.{trj_suffix}.{i:d}' for i in range(size)],
                                cpptraj)
            nmtraj.Setup(INPUT['nmode']['nmstartframe'], INPUT['nmode']['nmendframe'], INPUT['nmode']['nminterval'])
            # Now split up the complex trajectory by thread
            if nmtraj.processed_frames < size:
                raise xBFreE_Error('More processors than receptor nmode frames!')

            frames_per_rank = nmtraj.processed_frames // size
            extras = nmtraj.processed_frames - frames_per_rank * size
            frame_count = [frames_per_rank for i in range(size)]
            for i in range(size):
                if i < extras: frame_count[i] += 1
            last_frame = 1
            for i in range(size):
                frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
                nmtraj.Outtraj('receptor_nm.%s.%d' % (trj_suffix, i),
                               frames=frame_string, filetype=INPUT['general']['netcdf'])
                last_frame += frame_count[i]

            nmtraj.Run('rec_nm_traj_cpptraj.out')

            nmtraj = Trajectory(FILES.ligand_prmtop, [f'ligand.{trj_suffix}.{i:d}' for i in range(size)], cpptraj)
            nmtraj.Setup(INPUT['nmode']['nmstartframe'], INPUT['nmode']['nmendframe'], INPUT['nmode']['nminterval'])
            # Now split up the complex trajectory by thread
            if nmtraj.processed_frames < size:
                raise xBFreE_Error('More processors than ligand nmode frames!')

            frames_per_rank = nmtraj.processed_frames // size
            extras = nmtraj.processed_frames - frames_per_rank * size
            frame_count = [frames_per_rank for i in range(size)]
            for i in range(size):
                if i < extras: frame_count[i] += 1
            last_frame = 1
            for i in range(size):
                frame_string = '%d-%d' % (last_frame, last_frame + frame_count[i] - 1)
                nmtraj.Outtraj(f'ligand_nm.{trj_suffix}.{i:d}', frames=frame_string,
                               filetype=INPUT['general']['netcdf'])
                last_frame += frame_count[i]

            nmtraj.Run('lig_nm_traj_cpptraj.out')

        # end if not stability

    # end if INPUT['nmode']['nmoderun']
    return com_frames, rec_frames, lig_frames, num_frames_nmode


def make_mutant_trajectories(INPUT, FILES, rank, cpptraj, norm_sys, mut_sys):
    """ Mutates given trajectories and outputs dummy files for mutants """
    from xBFreE.mmpbsa.alamdcrd import MutantMdcrd, GlyMutantMdcrd
    import shutil
    if not INPUT['ala']['alarun']: return None, None

    stability = FILES.stability

    if INPUT['general']['netcdf']:
        trj_suffix = 'nc'
        raise TypeError('Alanine/Glycine scanning requires ASCII trajectories (netcdf=0)')
    else:
        trj_suffix = 'mdcrd'

    if (
        not stability
        and (
            FILES.ligand_prmtop != FILES.mutant_ligand_prmtop
            or FILES.receptor_prmtop == FILES.mutant_receptor_prmtop
        )
        and (
            FILES.ligand_prmtop == FILES.mutant_ligand_prmtop
            or FILES.receptor_prmtop != FILES.mutant_receptor_prmtop
        )
    ):
        raise xBFreE_Error('Alanine/Glycine scanning requires either a mutated ligand or receptor topology file with '
                           'only 1 mutant residue, but not both')

    master = rank == 0
    # FIXME: create folders for mutants

    # Have each rank mutate our rank's normal complex trajectory
    try:
        com_mut = MutantMdcrd('complex.%s.%d' % (trj_suffix, rank), norm_sys.complex_prmtop,
                              mut_sys.complex_prmtop)
    except MutantResError:
        com_mut = GlyMutantMdcrd(f'complex.{trj_suffix}.{rank:d}', norm_sys.complex_prmtop,
                                 mut_sys.complex_prmtop)
    com_mut.MutateTraj('mutant_complex.%s.%d' % (trj_suffix, rank))

    # Have each rank mutate our rank's normal receptor or ligand trajectory
    # and copy the normal one to the mutant if the mutated residue is *not*
    # present in there
    if not stability:
        if FILES.receptor_prmtop != FILES.mutant_receptor_prmtop:
            try:
                rec_mut = MutantMdcrd('receptor.%s.%d' % (trj_suffix, rank),
                                      norm_sys.receptor_prmtop, mut_sys.receptor_prmtop)
            except MutantResError:
                rec_mut = GlyMutantMdcrd('receptor.%s.%d' % (trj_suffix, rank),
                                         norm_sys.receptor_prmtop, mut_sys.receptor_prmtop)
            rec_mut.MutateTraj('mutant_receptor.%s.%d' % (trj_suffix, rank))
            shutil.copyfile('ligand.%s.%d' % (trj_suffix, rank),
                            'mutant_ligand.%s.%d' % (trj_suffix, rank))

        elif FILES.ligand_prmtop != FILES.mutant_ligand_prmtop:
            try:
                lig_mut = MutantMdcrd('ligand.%s.%d' % (trj_suffix, rank),
                                      norm_sys.ligand_prmtop, mut_sys.ligand_prmtop)
            except MutantResError:
                lig_mut = GlyMutantMdcrd('ligand.%s.%d' % (trj_suffix, rank),
                                         norm_sys.ligand_prmtop, mut_sys.ligand_prmtop)
            lig_mut.MutateTraj('mutant_ligand.%s.%d' % (trj_suffix, rank))
            shutil.copyfile('receptor.%s.%d' % (trj_suffix, rank),
                            'mutant_receptor.%s.%d' % (trj_suffix, rank))

    if INPUT['gbnsr6']['gbnsr6run']:
        temp_dir = Path(f"inpcrd_{rank}")
        if not temp_dir.exists():
            temp_dir.mkdir()
        com_traj = Trajectory(FILES.mutant_complex_prmtop, f"mutant_complex.{trj_suffix}.{rank}", cpptraj)
        com_traj.Setup()
        com_traj.Outtraj(f"inpcrd_{rank}/mutant_complex.inpcrd", filetype='restart', options=['keepext'])
        com_traj.Run('commutant_gbnsr6_traj_cpptraj.out')

        if not stability:
            # receptor
            rec_traj = Trajectory(FILES.mutant_receptor_prmtop, f"mutant_receptor.{trj_suffix}.{rank}", cpptraj)
            rec_traj.Setup()
            rec_traj.Outtraj(f"inpcrd_{rank}/mutant_receptor.inpcrd", filetype='restart', options=['keepext'])
            rec_traj.Run('recmutant_gbnsr6_traj_cpptraj.out')
            # ligand
            lig_traj = Trajectory(FILES.mutant_ligand_prmtop, f"mutant_ligand.{trj_suffix}.{rank}", cpptraj)
            lig_traj.Setup()
            lig_traj.Outtraj(f"inpcrd_{rank}/mutant_ligand.inpcrd", filetype='restart', options=['keepext'])
            lig_traj.Run('ligmutant_gbnsr6_traj_cpptraj.out')

    # Have our master dump out dummy files
    if master:
        com_traj = Trajectory(FILES.mutant_complex_prmtop, 'mutant_complex.%s.0' % trj_suffix, cpptraj)
        com_traj.Setup(1, 1, 1)
        com_traj.Outtraj('mutant_complex.pdb', frames='1', filetype='pdb')
        com_traj.Outtraj('mutant_dummycomplex.inpcrd', frames='1', filetype='restart')
        com_traj.Run('mutant_complex_cpptraj.out')
        if not stability:
            rec_traj = Trajectory(FILES.mutant_receptor_prmtop, 'mutant_receptor.%s.0' % trj_suffix, cpptraj)
            rec_traj.Setup(1, 1, 1)
            rec_traj.Outtraj('mutant_receptor.pdb', frames='1', filetype='pdb')
            rec_traj.Outtraj('mutant_dummyreceptor.inpcrd', frames='1', filetype='restart')
            rec_traj.Run('mutant_receptor_cpptraj.out')

            lig_traj = Trajectory(FILES.mutant_ligand_prmtop, 'mutant_ligand.%s.0' % trj_suffix, cpptraj)
            lig_traj.Setup(1, 1, 1)
            lig_traj.Outtraj('mutant_ligand.pdb', frames='1', filetype='pdb')
            lig_traj.Outtraj('mutant_dummyligand.inpcrd', frames='1', filetype='restart')
            lig_traj.Run('mutant_ligand_cpptraj.out')

    # Mutate our nmode trajectories if need be
    if INPUT['nmode']['nmoderun']:
        com_mut = MutantMdcrd('complex_nm.%s.%d' % (trj_suffix, rank),
                              norm_sys.complex_prmtop, mut_sys.complex_prmtop)
        com_mut.MutateTraj('mutant_complex_nm.%s.%d' % (trj_suffix, rank))
        if not stability and FILES.receptor_prmtop != FILES.mutant_receptor_prmtop:
            rec_mut = MutantMdcrd('receptor_nm.%s.%d' % (trj_suffix, rank),
                                  norm_sys.receptor_prmtop, mut_sys.receptor_prmtop)
            rec_mut.MutateTraj('mutant_receptor_nm.%s.%d' %
                               (trj_suffix, rank))
            shutil.copyfile('ligand_nm.%s.%d' % (trj_suffix, rank),
                            'mutant_ligand_nm.%s.%d' % (trj_suffix, rank))

        if not stability and FILES.ligand_prmtop != FILES.mutant_ligand_prmtop:
            lig_mut = MutantMdcrd('ligand_nm.%s.%d' % (trj_suffix, rank),
                                  norm_sys.ligand_prmtop, mut_sys.ligand_prmtop)
            lig_mut.MutateTraj('mutant_ligand_nm.%s.%d' %
                               (trj_suffix, rank))
            shutil.copyfile('ligand_nm.%s.%d' % (trj_suffix, rank),
                            'mutant_ligand_nm.%s.%d' % (trj_suffix, rank))

    # If we're doing a quasi-harmonic approximation we need the full com traj
    if (INPUT['general']['full_traj'] or INPUT['general']['qh_entropy']) and master:
        com_mut = MutantMdcrd('complex.%s' % trj_suffix, norm_sys.complex_prmtop, mut_sys.complex_prmtop)
        com_mut.MutateTraj('mutant_complex.%s' % trj_suffix)

    return str(com_mut), com_mut.mutres


class Trajectory(object):
    """ Base Trajectory class:
        Methods: __init__(prmtop, traj_files)
                 Setup(startframe, endframe, interval)
                 Query()
                 Image()
                 Strip(mask)
                 StripSolvent(strip_mask)
                 Unstrip(restrip_solvent=True)
                 Outtraj(fname, startframe, endframe, interval, filetype, nobox)
                 Run(output)

        The way this class is typically used is initializing it by giving it a
        list of trajectory files. It will then automatically run Query to make
        sure that all trajectory files are valid, and it counts how many frames
        are present in each trajectory. Then you can strip the solvent (this
        just tags what is the solvent in case you want to re-strip it after every
        Unstrip() command). Then, add masks to strip out. Then, you can run with
        the current actions via the Run() routine. Typical sequence:

        trajsystem = Trajectory(solvated_prmtop, mdcrds)
        trajsystem.Setup(startframe, endframe, interval)
        trajsystem.Image()
        trajsystem.StripSolvent(strip_mask)
        trajsystem.Strip(receptor_mask)
        trajsystem.Outtraj('_MMPBSA_ligand.mdcrd', filetype=INPUT['netcdf'])
        trajsystem.Unstrip(restrip_solvent=True)
        trajsystem.Strip(ligand_mask)
        trajsystem.Outtraj('_MMPBSA_receptor.mdcrd', filetype=INPUT['netcdf'])
        trajsystem.Run('_MMPBSA_create_trajectories.out')
    """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, prmtop, traj_files, cpptraj='cpptraj'):
        """ Sets up a basic trajectory """

        # Make sure we have a list or string, and force it to be a list
        if isinstance(traj_files, list):
            self.traj_files = traj_files
        elif isinstance(traj_files, tuple):
            self.traj_files = list(traj_files)
        else:
            self.traj_files = [str(traj_files)]

        self.prmtop = prmtop

        # Find cpptraj
        self.exe = cpptraj

        self.strip_solvent = False

        self.Query()

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Setup(self, startframe=1, endframe=9999999, interval=1):
        """ Finds the program and clears the queue and other attributes assigned
            during other function calls
        """

        self.strip_solvent = False
        orig_endframe = endframe
        orig_startframe = startframe

        # Delete the stripmask, but make its non-existence non-fatal, since it
        # won't exist unless Setup is called a second time.

        try:
            del self.stripmask
        except AttributeError:
            pass

        self.actions = []

        # Catch stupid choices
        if startframe > endframe:
            raise TrajError('Cannot have startframe (%d) > endframe (%d)' %
                            (startframe, endframe))

        if startframe < 0:
            raise TrajError('Startframe (%d) < 0' % startframe)

        if interval <= 0:
            raise TrajError('Interval (%d) <= 0' % interval)

        if startframe > self.total_frames:
            raise TrajError('start frame (%d) > total frames (%d)' %
                            (startframe, self.total_frames))

        # Bring endframe down if it's > the total number of frames

        endframe = min(endframe, self.total_frames)

        self.analyzed_frames = int((endframe - startframe) / interval) + 1

        # If we have an interval != 1, then it's possible that the endframe that
        # the user set is not the actual one that will be used. To make things
        # simpler, I will adjust the endframe to what is *actually* used

        endframe = startframe + interval * int((endframe - startframe) // interval)

        # Set up start and end arrays (for each trajectory)

        traj_starts = [-1 for i in range(len(self.traj_files))]
        traj_ends = [0 for i in range(len(self.traj_files))]

        # Now determine where each trajectory starts and ends

        for i in range(len(self.traj_files)):

            # skip over any trajectories that lie entirely before startframe

            if startframe > self.traj_sizes[i]:
                traj_starts[i] = -1  # this will tell us to skip this traj file
                startframe -= self.traj_sizes[i]
                continue

            # Now we start at our startframe

            traj_starts[i] = startframe

            # Now we figure out our last frame, adjusting for interval

            last_frame = startframe + interval * int((self.traj_sizes[i] - startframe) / interval)

            traj_ends[i] = min(endframe, last_frame, self.traj_sizes[i])

            # Now determine our new start frame for the next trajectory. First
            # we need to know how many frames are left over at the end of the
            # last frame, and subtract those from the interval to determine
            # where our next trajectory should start from. Also move our endframe
            # down, but only by our last_frame! (since the last couple frames in
            # our trajectory file could be leftover due to interval > 1)

            startframe = interval + last_frame - self.traj_sizes[i]

            endframe -= last_frame

            if endframe < self.traj_sizes[i]:
                break  # we're done now

        # end for i len(traj_files)

        # We now have to trajin each of the analyzed trajectories

        for i in range(len(self.traj_files)):
            if traj_starts[i] < 0: continue  # skip -1's

            self.actions.append('trajin %s %d %d %d' % (self.traj_files[i], traj_starts[i], traj_ends[i], interval))

        self.actions.append('noprogress')  # quash the progress bar

        self.processed_frames = (min(orig_endframe, self.total_frames) - orig_startframe) / interval + 1

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def StripSolvent(self, stripmask):
        """ Strips the solvent """
        self.stripmask = stripmask
        self.strip_solvent = True
        self.Strip(stripmask)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Query(self):
        """ Finds out how many frames are in the given traj files """
        from subprocess import Popen, PIPE
        import re

        framere = re.compile(r'Frames: (\d+)')

        self.traj_sizes = []

        # Determine how many frames are in each trajectory of the list

        for traj in self.traj_files:

            process = Popen([self.exe, '-p', str(self.prmtop), '-y', traj, '-tl'], stdin=PIPE, stdout=PIPE)

            (output, error) = process.communicate(b'')

            if process.wait():  # if it quits with return code != 0
                raise TrajError('%s failed when querying %s' % (self.exe, traj))

            output = output.decode()

            # Now parse the output to find out how many frames are there. We are
            # looking for "Coordinate processing will occur on x frames."

            num_frames = framere.findall(output)
            if len(num_frames) < 1:
                raise TrajError('Could not find number of frames in ' + traj)
            elif len(num_frames) > 1:
                raise RuntimeError('Unexpected output from cpptraj. Has format changed?')
            self.traj_sizes.append(int(num_frames[0]))

        # end for traj in self.traj_files

        self.total_frames = sum(self.traj_sizes)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Image(self):
        """ Images the trajectories """
        from parmed.amber import LoadParm

        solvated_prmtop = LoadParm(self.prmtop)

        ifbox = solvated_prmtop.ptr('ifbox')
        if ifbox == 0:
            warn('Solvated topology %s has IFBOX == 0' % ifbox)

        self.actions.append('autoimage')

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def rms(self, mask):
        """ Does an RMS fit around a specific mask """
        self.actions.append('rmsd %s mass first' % mask)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Strip(self, mask):
        """ Strips a mask from the coordinates. """
        if mask is None or len(str(mask).strip()) == 0:
            raise InternalError('Cannot pass no mask to Strip!')
        self.actions.append('strip %s' % mask)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Unstrip(self, restrip_solvent=True):
        """ Returns to unstripped state """
        self.actions.append('unstrip')
        if self.strip_solvent and restrip_solvent:
            self.actions.append('strip %s' % self.stripmask)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Run(self, output_file=None):
        """ Runs cpptraj to actually create the files """
        from sys import stdout as sys_stdout
        from subprocess import Popen, PIPE

        # Accept output_file as both a file object and string object
        own_handle = False
        if output_file is None:
            stdout = sys_stdout
        else:
            try:
                stdout = open(output_file, 'w')
                own_handle = True
            except TypeError:
                stdout = output_file

        # Now it's time to run the program

        try:
            input_string = ''
            for action in self.actions:
                input_string += action.strip() + '\n'

            process = Popen([self.exe, self.prmtop], stdout=stdout, stdin=PIPE)

            process.communicate(input_string.encode())

            if process.wait():
                raise TrajError('Error running %s' % self.exe)
        finally:
            if own_handle: stdout.close()


    def Outtraj(self, filename, frames=None, filetype='', nobox='nobox', options=None):
        """ This adds an outtraj command to the action stack, and you can specify
            the type of trajectory file to output (such as restart/pdb for input
            files, etc.)
        """
        if not frames: frames = '1-%d' % self.total_frames
        if options:
            self.actions.append(
                f"outtraj {filename} {filetype} onlyframes {frames} {nobox} {' '.join(options)}"
            )
        else:
            self.actions.append(
                f"outtraj {filename} {filetype} onlyframes {frames} {nobox}"
            )
