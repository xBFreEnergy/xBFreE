## **`&alanine_scanning` namelist variables**

!!! note "Keep in mind"
    * A default alanine scanning input file can be created as follows:

        === "GROMACS"
            ```bash
            xbfree gmx_MMPBSA --create_input ala                # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

        === "AMBER"
            ```bash
            xbfree amber_MMPBSA --create_input ala              # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

        === "NAMD"
            ```bash
            xbfree namd_MMPBSA --create_input ala               # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

        === "CHARMM"
            ```bash
            xbfree charmm_MMPBSA --create_input ala             # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

    * A sample alanine scanning input file is shown [here](input_file.md#alanine-scanning)
    * A tutorial on alanine scanning is available [here](examples/Alanine_scanning/README.md)

[`mutant_res`](#mmpbsa_ifv_mutant_res){#mmpbsa_ifv_mutant_res} (Default = None. Must be defined)
:   Define the specific residue that is going to be mutated. Use the following format CHAIN/RESNUM (_e.g._: A/350) or 
CHAIN/RESNUM INSERTION_CODE if applicable (_e.g._: A/27B).

    !!! important
        * Only one residue can be mutated per calculation!
        * For GROMACS, we recommend using the reference structure (-cr) to ensure the perfect match between the 
        selected residue in the defined structure or topology. For AMBER, NAMD and CHARMM we recommend check 
        carefully the structure provided (-cs). 
        * When this varibale is defined, **xBFreE** performs the mutation. This way the user does not have to 
        provide the mutant topology
    

[`mutant`](#mmpbsa_ifv_mutant){#mmpbsa_ifv_mutant} (Default = ALA) 
:   Defines the residue that it is going to be mutated for. Allowed values are: 

    * `ALA` or `A`: Alanine scanning
    * `GLY` or `G`: Glycine scanning

[`mutant_only`](#mmpbsa_ifv_mutant_only){#mmpbsa_ifv_mutant_only}  (Default = 0)
:   Option to perform specified calculations only for the mutants. 

    * 0: Perform calcultion on mutant and original
    * 1: Perform calcultion on mutant only
    
    !!! note
        Note that all calculation details are controlled in the other namelists, though for alanine scanning to be 
        performed, the namelist must be included (blank if desired)

[`cas_intdiel`](#mmpbsa_ifv_cas_intdiel){#mmpbsa_ifv_cas_intdiel} (Default = 0)
:   The dielectric constant (`intdiel`(GB)/`indi`(PB)) will be modified depending on the nature of the residue to be 
mutated. 
    
    * 0: Don’t
    * 1: Adaptative `intdiel` assignation

    !!! important
        * Works with the GB and PB calculations
        * It is ignored when `intdiel`(GB)/`indi`(PB) has been explicitly defined, that is, it is ignored if 
        `intdiel != 1.0`/`indi != 1.0` (default values)
        * Dielectric constant values has been assigned according to [Yan et al., 2017][9]
    !!! warning 
        Careful. Activating this variable can cause considerable variations in the results, since the dielectric 
        constant of the solute varies. 

  [9]: https://pubs.acs.org/doi/10.1021/acs.jcim.6b00734
    

[`intdiel_nonpolar`](#mmpbsa_ifv_intdiel_nonpolar){#mmpbsa_ifv_intdiel_nonpolar} (Default = 1)
:   Define the `intdiel`(GB)/`indi`(PB) value for non-polar residues (`PHE`, `TRP`, `VAL`, `ILE`, `LEU`, `MET`, `PRO`,
`CYX`, `ALA`, `GLY`, `PRO`)
    

[`intdiel_polar`](#mmpbsa_ifv_intdiel_polar){#mmpbsa_ifv_intdiel_polar} (Default = 3)
:   Define the `intdiel`(GB)/`indi`(PB) value for polar residues (`TYR`, `SER`, `THR`, `CYM`, `CYS`, `HIE`, `HID`, 
`ASN`, `GLN`, `ASH`, `GLH`, `LYN`)
    

[`intdiel_positive`](#mmpbsa_ifv_intdiel_positive){#mmpbsa_ifv_intdiel_positive} (Default = 5)
:   Define the `intdiel`(GB)/`indi`(PB) value for positive charged residues (`LYS`, `ARG`, `HIP`)
    

[`intdiel_negative`](#mmpbsa_ifv_intdiel_negative){#mmpbsa_ifv_intdiel_negative} (Default = 5)
:   Define the `intdiel`(GB)/`indi`(PB) value for negative charged residues (`GLU`, `ASP`)
