---
template: main.html
title: Protein-protein
---

# Protein-protein binding free energy calculations

!!! info
    This example can be found in the [examples/Protein_protein][6] directory in the repository folder. If you didn't
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    xBFreE GitHub repository.

## Requirements

Depending on the MD engine and force field used, `xBFreE` requires specific files:

=== "GROMACS"

    | Input File required            |                             Required                              |       Type        | Description                                                                                                      |
    |:-------------------------------|:-----------------------------------------------------------------:|:-----------------:|:-----------------------------------------------------------------------------------------------------------------|
    | Input parameters file          |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |       `in`        | Input file containing all the specifications regarding the type of calculation that is going to be performed     |
    | The MD Structure+mass(db) file |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |    `tpr` `pdb`    | Structure file containing the system coordinates                                                                 |
    | An index file                  |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |       `ndx`       | File containing the receptor and ligand in separated groups                                                      |
    | Receptor and ligand group      |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |    `integers`     | Group numbers in the index files                                                                                 |
    | A trajectory file              |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and with no pbc.                                                             |
    | A topology file                |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |       `top`       | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder                     |
    | A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |       `pdb`       | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |

=== "AMBER"

    | Input File required            |                             Required                              |        Type        | Description                                                                                                                                                                                                                 |
    |:-------------------------------|:-----------------------------------------------------------------:|:------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | Input parameters file          |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |        `in`        | Input file containing all the specifications regarding the type of calculation that is going to be performed                                                                                                                |
    | The MD Structure+mass(db) file |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |       `pdb`        | Structure file containing the system coordinates                                                                                                                                                                            |
    | Receptor and ligand masks      |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |      `range`       | The "mask" selection expressions start with ":". Residues can be selected by numbers (given as ranges separated by a dash) _e.g._ :1-30 :31-39                                                                                                |
    | A trajectory file              |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     | `nc` `pdb` `mdcrd` | Final AMBER MD trajectory, fitted and with no pbc.                                                                                                                                                                          |
    | A topology file                |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |  `parm7` `prmtop`  | AMBER topology file                                                                                                                                                                                                         |

=== "NAMD (amberff)"

    This is based on the standard output of CHARMM-GUI online server
    
    | Input File required            |                          Required                          |       Type       | Description                                                                                                                                    |
    |:-------------------------------|:----------------------------------------------------------:|:----------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------|
    | Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |       `in`       | Input file containing all the specifications regarding the type of calculation that is going to be performed                                   |
    | The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |      `pdb`       | Structure file containing the system coordinates                                                                                               |
    | Receptor and ligand masks      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |     `range`      | The "mask" selection expressions start with ":". Residues can be selected by numbers (given as ranges separated by a dash) _e.g._ :1-30 :31-39 |
    | A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |   `pdb` `dcd`    | Final NAMD MD trajectory, fitted and with no pbc.                                                                                              |
    | A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `parm7` `prmtop` | Topology file in amber format                                                                                                                  |       

=== "NAMD (charmmff)"

    This is based on the standard output of CHARMM-GUI online server
    
    | Input File required            |                          Required                          |    Type     | Description                                                                                                                                    |
    |:-------------------------------|:----------------------------------------------------------:|:-----------:|:-----------------------------------------------------------------------------------------------------------------------------------------------|
    | Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `in`     | Input file containing all the specifications regarding the type of calculation that is going to be performed                                   |
    | The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `pdb`    | Structure file containing the system coordinates                                                                                               |
    | Receptor and ligand masks      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |   `range`   | The "mask" selection expressions start with ":". Residues can be selected by numbers (given as ranges separated by a dash) _e.g._ :1-30 :31-39 |
    | A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `pdb` `dcd` | Final NAMD MD trajectory, fitted and with no pbc.                                                                                              |
    | A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `psf`    | Topology file in psf format                                                                                                                    |       

:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "GROMACS"
        
    === "Serial"

            xbfree gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -ci index.ndx -cg 3 4 -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

    === "With MPI"

            mpirun -np 2 xbfree gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -ci index.ndx -cg 3 4 -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "AMBER"

    === "Serial"

            xbfree amber_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-30 :31-39 -cp com.parm7 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

    === "With MPI"

            mpirun -np 2 xbfree amber_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-30 :31-39 -cp com.parm7 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "NAMD (amberff)"

    === "Serial"

            xbfree namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-30 :31-39 -cp com.parm7 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

    === "With MPI"

            mpirun -np 2 xbfree namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.pdb -cg :1-30 :31-39 -cp com.parm7 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv


=== "NAMD (charmmff)"

    === "Serial"

            xbfree namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.dcd -cg :1-30 :31-39 -cp com.psf -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

    === "With MPI"

            mpirun -np 2 xbfree namd_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.dcd -cg :1-30 :31-39 -cp com.psf -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv


where the `mmpbsa.in` input file, is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for GB calculation"
Sample input file for GB calculation
This input file is meant to show only that xBFreE works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Prot-Prot",
startframe=1,
endframe=10,
/
&gb
igb=2, saltcon=0.150,
/
```

!!! info "Keep in mind"
    See a detailed list of all the options in `xBFreE` input file [here][2] as well as several [examples][3]. 
    These examples are meant only to show that xBFreE works. It is recommended to go over these variables, even 
    the ones that are not included in this input file but are available for the calculation that it's performed and
    see the values they can take (check the [input file section](../../input_file.md)). This will allow you to 
    tackle a number of potential problems or simply use fancier approximations in your calculations.

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand (in this case, the 
ligand is also another protein) amber format topologies and trajectories will be obtained from that of the complex. To 
do so, an MD Structure+mass(db) file, a trajectory file, an index file (for those using GROMACS) and 
both the receptor and ligand group numbers or masks are needed. The `mmpbsa.in` input file will 
contain all the parameters needed for the MM/PB(GB)SA calculation. In this case, 10 frames 
are going to be used when performing the MM/PB(GB)SA calculation with the igb2 (GB-OBC1) model and a salt 
concentration = 0.15M.

A plain text output file with all the statistics (default: `FINAL_RESULTS_MMPBSA.dat`) and a CSV-format 
output file containing all energy terms for every frame in every calculation will be saved. The file name in 
'-eo' flag will be forced to end in [.csv] (`FINAL_RESULTS_MMPBSA.csv` in this case). This file is only written when 
specified on the command-line.

!!! note
    Once the calculation is done, the results can be analyzed in `xBFreE-Analyzer` (if `-nogui` flag was not used in the command-line). 
    Please, check the [xBFreE-Analyzer][4] section for more information

  [1]: ../../mmpbsa/command-line.md#gmx_mmpbsa-command-line  
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Protein_protein
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  
