---
template: main.html
title: The input file for xBFreE (MMPBSA)
---

# The input file

## Description

**xBFreE (MMPBSA)** input file contains all the specifications for the MMPBSA calculations. The input file is syntactically 
similar to other programs in Amber, although we have incorporated format more similar to that of GROMACS *.mdp 
files (see below). The input file contains sections called `namelist` where the variables are defined for each 
calculation. The allowed namelists are:

- [`&general`](general/#general-namelist-variables): contains variables that apply to all aspects of the 
  calculation or parameters required for building AMBER topologies.
- [`&gb`](gb/#gb-namelist-variables): unique variables to Generalized Born (GB) calculations.
- [`&gbnsr6`](gbnsr6/#gbnsr6-namelist-variables): unique variables to GBNSR6 calculations.
- [`&pb`](pb/#pb-namelist-variables): unique variables to Poisson Boltzmann (PB) calculations.
- [`&rism`](rism/#rism-namelist-variables): unique variables to 3D-RISM calculations.
- [`&alanine_scanning`](cas/#alanine_scanning-namelist-variables): unique variables to alanine scanning 
  calculations.
- [`&decomp`](decomp/#decomp-namelist-variables): unique variables to the decomposition scheme.
- [`&nmode`](nmodes/#nmode-namelist-variables): unique variables to the normal mode (NMODE) calculations used to 
  approximate vibrational entropies.

  [1]: https://pubs.acs.org/doi/10.1021/ct300418h

## Generation of input files with **xBFreE**
The input file can be created using **xBFreE** by selecting the subcommand and the calculations you want to perform.

!!! note
    Note that the command-line is basically the same for all MD programs. They only differ in the subcommand 
    selected according to the MD program.

    === "GROMACS"
        ``` title="Command-line"
        xbfree gmx_MMPBSA --create_input args
        ```
        where `args` can be:  `gb`, `gbnsr6`, `pb`, `rism`, `ala`, `decomp`, `nmode`, `all`
    
        Example:
        
        === "GB calculation"
                
                xbfree gmx_MMPBSA --create_input gb
            
        === "PB calculation"
            
                xbfree gmx_MMPBSA --create_input pb
        
        === "GB, PB and Decomposition calculations"
            
                xbfree gmx_MMPBSA --create_input gb pb decomp
        
        === "All calculations"
        
                xbfree gmx_MMPBSA --create_input
             
            or 
                
                xbfree gmx_MMPBSA --create_input all
    
    === "AMBER"
        ``` title="Command-line"
        xbfree amber_MMPBSA --create_input args
        ```
        where `args` can be:  `gb`, `gbnsr6`, `pb`, `rism`, `ala`, `decomp`, `nmode`, `all`
        
        Example:
        
        === "GB calculation"
                
                xbfree amber_MMPBSA --create_input gb
            
        === "PB calculation"
            
                xbfree amber_MMPBSA --create_input pb
        
        === "GB, PB and Decomposition calculations"
            
                xbfree amber_MMPBSA --create_input gb pb decomp
        
        === "All calculations"
        
                xbfree amber_MMPBSA --create_input
             
            or 
                
                xbfree amber_MMPBSA --create_input all
    
    === "NAMD"
        ``` title="Command-line"
        xbfree namd_MMPBSA --create_input args
        ```
        where `args` can be:  `gb`, `gbnsr6`, `pb`, `rism`, `ala`, `decomp`, `nmode`, `all`
    
        Example:
        === "GB calculation"
                
                xbfree namd_MMPBSA --create_input gb
            
        === "PB calculation"
            
                xbfree namd_MMPBSA --create_input pb
        
        === "GB, PB and Decomposition calculations"
            
                xbfree namd_MMPBSA --create_input gb pb decomp
        
        === "All calculations"
        
                xbfree namd_MMPBSA --create_input
             
            or 
                
                xbfree namd_MMPBSA --create_input all
    
    === "CHARMM"
        ``` title="Command-line"
        xbfree charmm_MMPBSA --create_input args
        ```
        where `args` can be:  `gb`, `gbnsr6`, `pb`, `rism`, `ala`, `decomp`, `nmode`, `all`
        
        Example:
        === "GB calculation"
                
                xbfree charmm_MMPBSA --create_input gb
            
        === "PB calculation"
            
                xbfree charmm_MMPBSA --create_input pb
        
        === "GB, PB and Decomposition calculations"
            
                xbfree charmm_MMPBSA --create_input gb pb decomp
        
        === "All calculations"
        
                xbfree charmm_MMPBSA --create_input
             
            or 
                
                xbfree charmm_MMPBSA --create_input all


!!! Danger 
    - Note that several variables must be explicitly defined in the input file
    - Performing all the calculations in the same run can take a long time

