
# The input file

## Description

`gmx_MMPBSA` input file contains all the specifications for the calculations. The input file is syntactically similar 
to other programs in Amber, although we incorporated a new format more similar to the one used on GROMACS *.mdp 
files (see bleow). The input file contains sections called `namelist` where the variables are defined for each 
calculation. The allowed namelists are:

- [`&general`](general/#general-namelist-variables): contains variables that apply to all aspects of the 
  calculation or parameters required for building AMBER topologies from GROMACS files.
- [`&gb`](gb/#gb-namelist-variables): unique variables to Generalized Born (GB) calculations.
- [`&gbnsr6`](gbnsr6/#gbnsr6-namelist-variables): unique variables to Generalized Born (GBNSR6) calculations.
- [`&pb`](pb/#pb-namelist-variables): unique variables to Poisson Boltzmann (PB) calculations.
- [`&delphi`](delphi/#delphi-namelist-variables): unique variables to Poisson Boltzmann (DelPhi) calculations.
- [`&rism`](rism/#rism-namelist-variables): unique variables to 3D-RISM calculations.
- [`&alanine_scanning`](cas/#alanine_scanning-namelist-variables): unique variables to alanine scanning 
  calculations.
- [`&decomp`](decomp/#decomp-namelist-variables): unique variables to the decomposition scheme.
- [`&nmode`](nmodes/#nmode-namelist-variables): unique variables to the normal mode (NMODE) calculations used to 
  approximate vibrational entropies.

  [1]: https://pubs.acs.org/doi/10.1021/ct300418h

## Generation of input files with gmx_MMPBSA
The input file can be created using `gmx_MMPBSA` by selecting the calculations you want to perform.

``` title="Command-line"
gmx_MMPBSA --create_input args

where `args` can be:  gb, gbnsr6, pb, rism, ala, decomp, nmode, all
```

Example:
=== "GB calculation"
        
        gmx_MMPBSA --create_input gb
    
=== "PB calculation"
    
        gmx_MMPBSA --create_input pb

=== "GB, PB and Decomposition calculations"
    
        gmx_MMPBSA --create_input gb pb decomp

=== "All calculations"

        gmx_MMPBSA --create_input
     
    or 
        
        gmx_MMPBSA --create_input all
        
!!! Danger 
    Note that several variables must be explicitly defined in the input file

_Implemented in v1.5.0_

## Format
All the input variables are described below according to their respective namelists. Descriptions are taken from 
original sources and modified accordingly when needed. Integers and floating point 
variables should be typed as-is while strings should be put in either single- or double-quotes. All variables should be 
set with `variable = value` and separated by commas is they appear in the same line. If the variables appear in different 
lines, the comma is no longer needed. See several [examples](#sample-input-files) are shown below. As you will see, several 
calculations can be performed in the same run (_i.e._ `&gb` and `&pb`; `&gb` and `&alanine_scanning`; `&pb` and
`&decomp`; etc). As we have mentioned, the input file can be generated using the `create_input` option of `gmx_MMPBSA`. 
This style, while retaining the same Amber format (derived from Fortran), is aesthetically more familiar to the GROMACS
style (`*.mdp`). However, it maintains the same essence, so it could be defined in any of the two format styles or even
combined. See the formats below:

=== "New format style "
    ``` title="New format style Input file example"
            
    # General namelist variables
    &general
      sys_name             = ""                      # System name
      startframe           = 1                       # First frame to analyze
      endframe             = 9999999                 # Last frame to analyze
      ...
      interval              = 1                      # The offset from which to choose frames from each trajectory file
    /
    
    # Generalized-Born namelist variables
    &gb
      igb                  = 5                       # GB model to use
      ...
      probe                = 1.4                     # Solvent probe radius for surface area calc
    /
    ```

=== "Old format style"
    ``` title="Old format style Input file example"
            
    # General namelist variables
    &general
      sys_name = "", startframe = 1, endframe = 9999999
      ...
      interval = 1
    /
    
    # Generalized-Born namelist variables
    &gb
      igb = 5, 
      ...
      probe = 1.4
    /
    ```
