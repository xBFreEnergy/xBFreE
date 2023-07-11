---
template: main.html
title:
---

## **xBFreE** command-line
**xBFreE** uses subcommands to execute different BFE methods for different MD programs:

<div class="termy">
    ```console
    // All flags available in **xBFreE** are shown below:

    $ xbfree -h
    
    usage: xbfree [-h] [-v] {gmx_MMPBSA,amber_MMPBSA,namd_MMPBSA,charmm_MMPBSA} ...
    
    xBFreEnergy is a tool to compute Binding Free Energy with different methods
    
    positional arguments:
      {gmx_MMPBSA,amber_MMPBSA,namd_MMPBSA,charmm_MMPBSA}
                            Methods to compute Binding Free Energy
        gmx_MMPBSA          PB and other implicit solvent-based calculations for GROMACS
        amber_MMPBSA        PB and other implicit solvent-based calculations for AMBER
        namd_MMPBSA         PB and other implicit solvent-based calculations for NAMD
        charmm_MMPBSA       PB and other implicit solvent-based calculations for CHARMM
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
    
    xBFreEnergy is a tool to compute Binding Free Energy with different methods
    ```
    

</div>

## Sub-commands structure
Each subcommand contains its own arguments according to the MD program and method selected. However, we decide to make 
some arguments redundant to simplify the command-line instructions and keep it similar to gmx_MMPBSA:

=== "GROMACS"
    
    <div class="termy">
        ```console
        // All flags available in **xBFreE** are shown below:
        
        $ xbfree gmx_MMPBSA -h

        usage: ...

        This is the core of gmx_MMPBSA and it will perfrom all the calculations
        
        optional arguments:
          -h, --help            show this help message and exit
          -O, --overwrite       Allow output files to be overwritten
          --create_input [{gb,pb,rism,ala,decomp,nmode,gbnsr6,all} ...]
                                Create an new input file with selected calculation type.
          --rewrite-output      Do not re-run any calculations, just parse the output files from the previous calculation and rewrite the output files.
          -s, --stability       Perform stability calculation. Only the complex parameters are required. In any other case receptor and ligand parameters will be ignored
          -ng, --nogui          No open gmx_MMPBSA_ana after all calculations finished
          -v, --version         show program's version number and exit
        
        Input and Output Files:
          These options specify the input files and optional output files.
        
          -i FILE               MM/PBSA input file.
          -xvvfile XVVFILE      XVV file for 3D-RISM.
          -o FILE               Output file with statistics for the selected method.
          -do FILE              Output file for decomposition statistics summary.
          -eo FILE              CSV-format output of all energy terms for every frame in every calculation. File name forced to end in [.csv]. This file is only written when specified on the command-line.
          -deo FILE             CSV-format output of all energy terms for each printed residue in decomposition calculations. File name forced to end in [.csv]. This file is only written when specified on the command-line.
        
        < Gromacs specific arguments >
        
        Miscellaneous Options:
          -prefix <file prefix>
                                Prefix for intermediate files.
          --input-file-help     Print all available options in the input file.
          --clean               Clean temporary files and quit.
        
        xBFreEnergy is a tool to compute Binding Free Enrgy with different methods
        ```
    </div>

=== "AMBER"
    <div class="termy">
        ```console
        // All flags available in **xBFreE** are shown below:
        
        $ xbfree amber_MMPBSA -h

        usage: ...

        This is the core of amber_MMPBSA and it will perform all the calculations
        
        optional arguments:
          -h, --help            show this help message and exit
          -O, --overwrite       Allow output files to be overwritten
          --create_input [{gb,pb,rism,ala,decomp,nmode,gbnsr6,all} ...]
                                Create an new input file with selected calculation type.
          --rewrite-output      Do not re-run any calculations, just parse the output files from the previous calculation and rewrite the output files.
          -s, --stability       Perform stability calculation. Only the complex parameters are required. In any other case receptor and ligand parameters will be ignored
          -ng, --nogui          No open gmx_MMPBSA_ana after all calculations finished
          -v, --version         show program's version number and exit
        
        Input and Output Files:
          These options specify the input files and optional output files.
        
          -i FILE               MM/PBSA input file.
          -xvvfile XVVFILE      XVV file for 3D-RISM.
          -o FILE               Output file with statistics for the selected method.
          -do FILE              Output file for decomposition statistics summary.
          -eo FILE              CSV-format output of all energy terms for every frame in every calculation. File name forced to end in [.csv]. This file is only written when specified on the command-line.
          -deo FILE             CSV-format output of all energy terms for each printed residue in decomposition calculations. File name forced to end in [.csv]. This file is only written when specified on the command-line.
        
        < Amber specific arguments >
        
        Miscellaneous Options:
          -prefix <file prefix>
                                Prefix for intermediate files.
          --input-file-help     Print all available options in the input file.
          --clean               Clean temporary files and quit.
        
        xBFreEnergy is a tool to compute Binding Free Enrgy with different methods
        ```
    </div>

=== "NAMD"
    <div class="termy">
        ```console
        // All flags available in **xBFreE** are shown below:
        
        $ xbfree namd_MMPBSA -h

        usage: ...

        This is the core of namd_MMPBSA and it will perform all the calculations
        
        optional arguments:
          -h, --help            show this help message and exit
          -O, --overwrite       Allow output files to be overwritten
          --create_input [{gb,pb,rism,ala,decomp,nmode,gbnsr6,all} ...]
                                Create an new input file with selected calculation type.
          --rewrite-output      Do not re-run any calculations, just parse the output files from the previous calculation and rewrite the output files.
          -s, --stability       Perform stability calculation. Only the complex parameters are required. In any other case receptor and ligand parameters will be ignored
          -ng, --nogui          No open gmx_MMPBSA_ana after all calculations finished
          -v, --version         show program's version number and exit
        
        Input and Output Files:
          These options specify the input files and optional output files.
        
          -i FILE               MM/PBSA input file.
          -xvvfile XVVFILE      XVV file for 3D-RISM.
          -o FILE               Output file with statistics for the selected method.
          -do FILE              Output file for decomposition statistics summary.
          -eo FILE              CSV-format output of all energy terms for every frame in every calculation. File name forced to end in [.csv]. This file is only written when specified on the command-line.
          -deo FILE             CSV-format output of all energy terms for each printed residue in decomposition calculations. File name forced to end in [.csv]. This file is only written when specified on the command-line.
        
        < NAMD specific arguments >
        
        Miscellaneous Options:
          -prefix <file prefix>
                                Prefix for intermediate files.
          --input-file-help     Print all available options in the input file.
          --clean               Clean temporary files and quit.
        
        xBFreEnergy is a tool to compute Binding Free Enrgy with different methods
        ```
    </div>

=== "CHARMM"
    <div class="termy">
        ```console
        // All flags available in **xBFreE** are shown below:
        
        $ xbfree charmm_MMPBSA -h

        usage: ...

        This is the core of charmm_MMPBSA and it will perform all the calculations
        
        optional arguments:
          -h, --help            show this help message and exit
          -O, --overwrite       Allow output files to be overwritten
          --create_input [{gb,pb,rism,ala,decomp,nmode,gbnsr6,all} ...]
                                Create an new input file with selected calculation type.
          --rewrite-output      Do not re-run any calculations, just parse the output files from the previous calculation and rewrite the output files.
          -s, --stability       Perform stability calculation. Only the complex parameters are required. In any other case receptor and ligand parameters will be ignored
          -ng, --nogui          No open gmx_MMPBSA_ana after all calculations finished
          -v, --version         show program's version number and exit
        
        Input and Output Files:
          These options specify the input files and optional output files.
        
          -i FILE               MM/PBSA input file.
          -xvvfile XVVFILE      XVV file for 3D-RISM.
          -o FILE               Output file with statistics for the selected method.
          -do FILE              Output file for decomposition statistics summary.
          -eo FILE              CSV-format output of all energy terms for every frame in every calculation. File name forced to end in [.csv]. This file is only written when specified on the command-line.
          -deo FILE             CSV-format output of all energy terms for each printed residue in decomposition calculations. File name forced to end in [.csv]. This file is only written when specified on the command-line.
        
        < CHARMM specific arguments >
        
        Miscellaneous Options:
          -prefix <file prefix>
                                Prefix for intermediate files.
          --input-file-help     Print all available options in the input file.
          --clean               Clean temporary files and quit.
        
        xBFreEnergy is a tool to compute Binding Free Enrgy with different methods
        ```
    </div>