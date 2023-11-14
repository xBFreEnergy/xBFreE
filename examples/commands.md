# commands:

# Protein-protein

## gmx_MMPBSA

python ../../../run_cmd.py gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -ci index.ndx -cg 3 4 -cp topol.top

## amber_MMPBSA

python ../../../run_cmd.py amber_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-30 :31-39 -cp com.parm7

## namd_MMPBSA (amber)

python ../../../run_cmd.py namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-30 :31-39 -cp com.parm7

## namd_MMPBSA (charmm)

python ../../../run_cmd.py namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.dcd -cg :1-30 :31-39 -cp com.psf


# Protein-ligand

## gmx_MMPBSA

python ../../../run_cmd.py gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -ci index.ndx -cg 1 18 -cp topol.top

## amber_MMPBSA

python ../../../run_cmd.py amber_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-58 :59 -cp com.parm7

## namd_MMPBSA (amber)

python ../../../run_cmd.py namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-58 :59 -cp com.parm7

## namd_MMPBSA (charmm)

python ../../../run_cmd.py namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.dcd -cg :1-58 :59 -cp com.psf
