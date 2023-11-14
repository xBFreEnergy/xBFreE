    

## **`&decomp` namelist variables**

!!! note "Keep in mind"
    * A default decomp input file can be created as follows:

        === "GROMACS"
            ```bash
            xbfree gmx_MMPBSA --create_input decomp         # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

        === "AMBER"
            ```bash
            xbfree amber_MMPBSA --create_input decomp           # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

        === "NAMD"
            ```bash
            xbfree namd_MMPBSA --create_input decomp            # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

        === "CHARMM"
            ```bash
            xbfree charmm_MMPBSA --create_input decomp          # (1)! 
            ```
            
            1.  Remember to define an energy method, for example `GB`, `PB`, etc.

    
    * A sample decomp input file is shown [here](input_file.md#decomposition-analysis)
    * A tutorial on binding free energy decomposition is available [here](examples/Decomposition_analysis/README.md)

[`idecomp`](#mmpbsa_ifv_idecomp){#mmpbsa_ifv_idecomp} (Default = None. Must be defined)
:   Energy decomposition scheme to use:
    
    * 1: Per-residue decomp with 1-4 terms added to internal potential terms
    * 2: Per-residue decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms
    * 3: Pairwise decomp with 1-4 terms added to internal potential terms
    * 4: Pairwise decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms

[`dec_verbose`](#mmpbsa_ifv_dec_verbose){#mmpbsa_ifv_dec_verbose} (Default = 0)
:   Set the level of output to print in the decomp_output file.

    * 0: DELTA energy, total contribution only
    * 1: DELTA energy, total, sidechain, and backbone contributions
    * 2: Complex, Receptor, Ligand, and DELTA energies, total contribution only
    * 3: Complex, Receptor, Ligand, and DELTA energies, total, sidechain, and backbone contributions

    !!! note
        If the values 0 or 2 are chosen, only the Total contributions are required, so only those will be printed to the
        mdout files to cut down on the size of the mdout files and the time required to parse them.

[`print_res`](#mmpbsa_ifv_print_res){#mmpbsa_ifv_print_res} (Default = "within 6")
:   Select residues whose information is going to be printed in the output file. The default selection should be 
sufficient in most cases, however we have added several additional notations
    
    !!! example "Selection schemes"

        === "By Distance (recommended)"
            Notation: [ `within` `distance` ]
            :   `within` corresponds to the keyword and `distance` to the maximum distance criterion in Å necessary to 
                select the residues from both the receptor and the ligand. In case the cutoff used is so small that 
                the number of decomp residues to print < 2, the cutoff value will be increased by 0.1 until the 
                number of decomp residues to print >= 2.
    
            !!! example
                `print_res="within 6"` Will print all residues within 6 Å between receptor and 
                ligand including both.
    
        === "Amino acid selection"
            Notation: [ `CHAIN`/(`RESNUM` or `RESNUM-RESNUM`) ]
            :   Print residues individual or ranges. This notation also supports insertion codes, in which case you must 
                define them individually
    
            !!! example
                `print_res="A/1,3-10,15,100 B/25"` This will print Chain A residues 1, 3 through 10, 15, and 100 along with 
                chain B residue 25 from the complex topology file and the corresponding residues in either the ligand and/or 
                receptor topology files.
    
                Suppost that we can have the following sequence where chain A is the receptor and B is the ligand: 
                A:LEU:5, A:GLY:6A, A:THR:6B, A:SER:6C A:ASP:6D, A:ILE:7 , B:25
                
                === "Supported notation"
                    
                    **Ranges selection**
                    :   `print_res="A/5-7 B/25` Will print all mentioned residues because all residues with insertion code are 
                        contained in the range
                    
                    **Individual selection**
                    :   `print_res="A/5,6B,6C,7 B/25` Will print all mentioned residues except the residues 6A and 
                        6D from chain A
    
                === "Wrong notation"
                    `print_res="A/5-6B,6D-7` Will end in error.
    
        === "All"
    
            Notation: `all`
            :   will print all residues. This option is often not recommended since most residues contribution is zero and 
                it is just going to be a waste of time and computational resources.
    
            !!! danger
                Using `idecomp=3 or 4` (pairwise) with a very large number of printed residues and a large number of frames 
                can quickly create very, very large temporary mdout files. Large print selections also demand a large amount 
                of memory to parse the mdout files and write decomposition output file (~500 MB for just 250 residues, since 
                that’s 62500 pairs!) It is not unusual for the output file to take a significant amount of time to print if 
                you have a lot of data. This is most applicable to pairwise decomp, since the amount of data scales as  
                O(N^2^).
     
    !!! important
        For GROMACS, we recommend using the reference structure (-cr) to ensure the perfect match between the 
        selected residue in the defined structure or topology. For AMBER, NAMD and CHARMM we recommend check 
        carefully the structure provided (-cs). 


[`csv_format`](#mmpbsa_ifv_csv_format){#mmpbsa_ifv_csv_format}  (Default = 1)
:   Print the decomposition output in a Comma-Separated-Values (CSV) file. CSV files open natively in most
spreadsheets. 

    * 0: data to be written out in the standard ASCII format.
    * 1: data to be written out in a CSV file, and standard error of the mean will be calculated and included for all 
    data.
