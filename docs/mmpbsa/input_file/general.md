## **`&general` namelist variables**

### **Basic input options**

[`sys_name`](#mmpbsa_ifv_sys_name){#mmpbsa_ifv_sys_name} (Default = None)
:   Define the System Name. This is useful when trying to analyze several systems at the same time or calculating 
the correlation between the predicted and the experimental energies. If the name is not defined, one will be assigned 
when loading the system in `gmx_MMPBSA_ana` on the loading order.

    !!! tip 
        The definition of the system name is entirely optional, however it can provide a better clarity during 
        the results analysis. All files associated with this system will be saved using its name.

[`startframe`](#mmpbsa_ifv_startframe){#mmpbsa_ifv_startframe} (Default = 1)
:   The frame from which to begin extracting snapshots from the full, concatenated trajectory comprised of
every trajectory file placed on the command-line. This is always the first frame read.

[`endframe`](#mmpbsa_ifv_endframe){#mmpbsa_ifv_endframe} (Default = 9999999)
:   The frame from which to stop extracting snapshots from the full, concatenated trajectory comprised of every
trajectory file supplied on the command-line.

[`interval`](#mmpbsa_ifv_interval){#mmpbsa_ifv_interval} (Default = 1)
:   The offset from which to choose frames from each trajectory file. For example, an interval of 2 will pull
every 2nd frame beginning at startframe and ending less than or equal to endframe.

### **Parameter options**

[`temperature`](#mmpbsa_ifv_temperature){#mmpbsa_ifv_temperature} (Default = 298.15)  
:   Specify the temperature (in K) used in the calculations.

[`PBRadii`](#mmpbsa_ifv_PBRadii){#mmpbsa_ifv_PBRadii} (Default = "mbondi2")
:   PBRadii is the parameter that defines the radius that will be assigned to each atom during the calculation of 
the solvation energy. You can combine multiple PBRadii for the same system!

    ??? warning "Effect of radii on energy calculations"
        Depending on the method selected, this parameter will have a greater or lesser impact on the computed value. 
        While in PB, this will only be used to compute the non-polar solvation component (`ENPOLAR` and `EDISPER`). In 
        GB, it is used, in addition to the non-polar solvation component, to compute the effective Born radius. 

    ??? gmx-mmpbsa "For gmx_MMPBSA users!"
        Note that notation changes from number to string. We implemented a new function to assing radii, which allow 
        customs radii defined by de user (for example, add Au radii to the `mbondi` radii) through the file path.    

    * `bondi`
        
        ??? info "`bondi` radii set"
                        
            |               | Description                                                                                                                                                                       |
            |---------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended when `igb = 7`                                                                                                                                                        |
            | Compatibility | `bondi` is a generic radii type, which means that the radii is assigned by atom properties for any type of molecule                                                               |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)                                                                         |
            | Combination   | since `bondi` is a generic type, you can combine it with specific radii set (`tyl06` or `yamagishi`). Note you can't combine it with other generic radii or complete radii types. |

    * `mbondi`
        
        ??? info "`mbondi` radii set"
                        
            |               | Description                                                                                                                                                                             |
            |---------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended when `igb = 1`                                                                                                                                                              |
            | Compatibility | `mbondi` is a generic radii type, which means that the radii is assigned by atom properties for any type of molecule                                                                    |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)                                                                               |
            | Combination   | since `mbondi` is a generic type, you can combine it with specific radii set (`tyl06` or `yamagishi`). Note you can't combine it with other generic radii or complete radii types. |
                
    * `mbondi2`
            
        ??? info "`mbondi2` radii set"

            |               | Description                                                                                                                                                                                   |
            |---------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended when `igb = 2` or `igb = 5`                                                                                                                                                       |
            | Compatibility | `mbondi2` is a generic radii type, which means that the radii is assigned by atom properties for any type of molecule                                                                         |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)                                                                                     |
            | Combination   | since `mbondi2` is a generic type, you can combine it with specific radii set (`tyl06` or `yamagishi`). Note you can't combine it with other generic radii or complete radii types. |

    * `mbondi3`

        ??? info "`mbondi3` radii set"
            
            |               | Description                                                                                                                                                                                        |
            |---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended when `igb = 8`                                                                                                                                                                         |
            | Compatibility | `mbondi3` is a generic radii type, which means that the radii is assigned by atom properties for any type of molecule                                                                              |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)                                                                                          |
            | Combination   | since `mbondi3` is a generic type, you can combine it with specific radii set (`tyl06` or `yamagishi`). Note you can't combine it with other generic radii or complete radii types. |

    * `mbondi_pb2`

        ??? info "`mbondi_pb2` radii set"
            
            |               | Description                                                                                                                                                                                           |
            |---------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended when `igb = 1`                                                                                                                                                                            |
            | Compatibility | `mbondi_pb2` is a generic radii type, which means that the radii is assigned by atom properties for any type of molecule                                                                              |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)                                                                                             |
            | Combination   | since `mbondi_pb2` is a generic type, you can combine it with specific radii set (`tyl06` or `yamagishi`). Note you can't combine it with other generic radii or complete radii types. |

            It is based on _mbondi_ radii set and contains a  new optimized set of halogen PB radii for halogenated 
            compounds (without extra point (EP) of charge) parametrized with General Amber Force Field (GAFF):

            Values from Table 3 in [§3.1 Halogen Radii Optimization Without EP][300]:

            * Cl: 1.76
            * Br: 1.97
            * I: 2.09
    
            This radii set should be used with the following PBSA setup:
    
            ```
            Sample input file for PB calculation with halogenated compounds
            
            &general
            sys_name="PB_Halogens",
            PBRadii="mbondi_pb2",
            /
            &pb
            radiopt=0, istrng=0.150, inp=1,
            /
            ```

    * `mbondi_pb3`

        ??? info "`mbondi_pb3` radii set"
            
            |               | Description                                                                                                                                                                                 |
            |---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended when `igb = 1`                                                                                                                                                                  |
            | Compatibility | `mbondi_pb3` is a generic radii type, which means that the radii is assigned by atom properties for any type of molecule                                                                    |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)                                                                                   |
            | Combination   | since `mbondi_pb3` is a generic type, you can combine it with specific radii set (`tyl06` or `yamagishi`). Note you can't combine it with other generic radii or complete radii types. |

            It is based on _mbondi_ radii set and contains a new optimized set of halogen PB radii for halogenated 
            compounds (without extra point (EP) of charge) parametrized with General Amber Force Field (GAFF):

            Values from Table 3 in [§3.1 Halogen Radii Optimization Without EP][300]:

            * Cl: 2.20
            * Br: 2.04
            * I: 2.19
                    

            This radii set should be used with the following PBSA setup:
    
            ```
            Sample input file for PB calculation with halogenated compounds
            
            &general
            sys_name="PB_Halogens",
            PBRadii="mbondi_pb3",
            /
            &pb
            radiopt=0, istrng=0.150, inp=1,
            /
            ```
            
    [300]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106

    * `charmm_radii`

        ??? info "`charmm_radii` radii set"

            |               | Description                                                                                                                         |
            |---------------|-------------------------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended for CHARMM force field only. Compatible with &pb only                                                                   |
            | Compatibility | `charmm_radii` is a complete radii type, which means that it contain radii for protein, nucleic acids, ligands and lipids molecules |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)                           |
            | Combination   | can't be combined                                                                                                                   |

            This atomic radii set for Poisson-Boltzmann calculations has been derived from average solvent 
            electrostatic charge distribution with explicit solvent. The accuracy has been tested with free energy 
            perturbation with explicit solvent. Most of the values were taken from a _*radii.str_ file used in PBEQ 
            Solver in [charmm-gui](https://www.charmm-gui.org/?doc=input/pbeqsolver).

            * Radii for protein atoms in 20 standard amino acids from [Nina, Belogv, and Roux](https://pubs.acs.org/doi/10.1021/jp970736r)
            * Radii for nucleic acid atoms (RNA and DNA) from [Banavali and Roux](https://pubs.acs.org/doi/abs/10.1021/jp025852v)
            * Halogens and other atoms from [Fortuna and Costa](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
    
    * `tyl06` 
        
        ??? info "`tyl06` radii set"
                        
            |               | Description                                                                                                      |
            |---------------|------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended for protein and nucleic acids only. Recommended when `inp=2`                                         |
            | Compatibility | `tyl06` is a semi-complete radii type, which means that it contain radii for protein and nucleic acids molecules |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)        |
            | Combination   | can be combined with any generic radii set, for example, mbondi, mbondi2, etc.                                   |

            The [Tan, Yang & Luo radii][255] are optimized for Amber atom types as in standard residues from the Amber 
            database. Please see [the original study][255] on how these radii are optimized.
        
        ??? note "For Amber users!"
            Note that this radii is the same applied when `radiopt=1`. However, as radiopt is not available in xbfree 
            due to possible fails if the system contains other molecules types in addition to proteins and nucleic 
            acids, then you have to apply it manually.  

    * `yamagishi` 
        
        ??? info "`yamagishi` radii set"
                                    
            |               | Description                                                                                                      |
            |---------------|------------------------------------------------------------------------------------------------------------------|
            | Recommended   | recommended for protein and nucleic acids only.                                                                  |
            | Compatibility | `yamagishi` is a semi-complete radii type, which means that it contain radii for protein and nucleic acids molecules |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#)        |
            | Combination   | can be combined with any generic radii set, for example, mbondi, mbondi2, etc.                                                                                                |

            The [Yamgishi radii][256] are optimized for Amber atom types as in standard residues from the Amber 
            database. Please see [the original study][256] on how these radii are optimized.

    * `custom` 
        
        ??? info "Defining a `custom` radii set"
                                                
            |               | Description                                                                                               |
            |---------------|-----------------------------------------------------------------------------------------------------------|
            | Recommended   | it dependents of radii set type.                                                                          |
            | Compatibility | it dependents of radii set type.                                                                          |
            | Modification  | you can add parameters for new atoms or modify existing ones. Please see how to do it [here](pbradii.md#) |
            | Combination   | it dependents of radii set type.                                                                          |

            You can define a new custom radii set defining the path to the radii file (*.json) as follows:
            `PBRadii = /home/user/wdir/custom_radii.json`
            This allows you to define specific parameters for new atoms, change existing ones, or define your set of 
            optimized radii. Please see ["How to create a custom radio set"](pbradii.md#)
        

    [255]: https://pubs.acs.org/doi/10.1021/jp063479b
    [256]: https://doi.org/10.1002/jcc.23728

### **Entropy options**

[`qh_entropy`](#mmpbsa_ifv_qh_entropy){#mmpbsa_ifv_qh_entropy} (Default = 0)
:    It specifies whether to perform a quasi-harmonic entropy (QH) approximation with `cpptraj` or not.
     
     * 0: Don’t
     * 1: perform QH

    !!! important "Keep in mind"
        * The number of frames used for QH analyses should be higher than 3N, N being the number of atoms in the 
        complex
        * Check this [thread](http://archive.ambermd.org/201207/0319.html) for more info on QH analysis

[`interaction_entropy`](#mmpbsa_ifv_interaction_entropy){#mmpbsa_ifv_interaction_entropy} (default = 0)
:    It specifies whether to use the [Interaction Entropy (IE)][3] approximation.
     
     * 0: Don’t
     * 1: perform IE

    !!! note "Keep in mind"
        - The Interaction Entropy can be calculated independently of the solvent model used.
        - A sample Interaction Entropy input file is shown [here](input_file.md#interaction-entropy)
        - A tutorial on the use of Interaction Entropy is 
        available [here](examples/Entropy_calculations/Interaction_Entropy/README.md)
        - The standard deviation of the interaction energy (`σIE`) should always be reported when using the Interaction 
        Entropy method.
        - The Interaction Entropy method should be avoided if `σIE > ~ 3.6 kcal/mol` because it is impossible to 
        converge the exponential average.
        - It is advisable to study how the Interaction Entropy depends on N by block averaging (which also provide an 
        estimate of the precision of the calculated entropies).
        - A sampling frequency of 10 fs, as reported in the original [IE publication][3], seems to be 3–40 times too 
        dense. A sampling frequency of 0.1 ps would be more appropriate.
        - The Interaction Entropy results may vary depending on the system flexibility or whether constraints were used 
        or not in the MD simulation. 

        Please, check this [paper][10] for further details.

  [3]: https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682
  [10]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00374


[`ie_segment`](#mmpbsa_ifv_ie_segment){#mmpbsa_ifv_ie_segment} (Default = 25)
:    Representative segment (in %), starting from the last frame, for the calculation of the
Interaction Entropy, _e.g._: `ie_segment = 25` means that the last quartile of the total number of frames
(`(endframe-startframe)/interval`) will be used to calculate the average Interaction Entropy.

[`c2_entropy`](#mmpbsa_ifv_c2_entropy){#mmpbsa_ifv_c2_entropy} (default = 0) 
:    It specifies whether to use the [C2 Entropy][11] approximation.
     
     * 0: Don’t
     * 1: perform C2

    !!! note "Keep in mind"
        - The C2 Entropy can be calculated independently of the solvent model used.
        - A tutorial on the use of C2 Entropy is 
        available [here](examples/Entropy_calculations/C2_Entropy/README.md)
        - The standard deviation of the interaction energy (`σIE`) should always be reported.
        - The C2 Entropy method should be avoided if `σIE > ~ 3.6 kcal/mol` because it gives unrealistically large 
        entropies.
        - It is advisable to study how the C2 Entropy depends on N by block averaging (which also provide an 
        estimate of the precision of the calculated entropies).
        - A sampling frequency of 10 fs, seems to be 3–40 times too dense. A sampling frequency of 0.1 ps would be more 
        appropriate.
        - The C2 Entropy results may vary depending on the system flexibility or whether constraints were used 
        or not in the MD simulation.

        Please, check this [paper][10] for further details.

  [10]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00374
  [11]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00418

### **Miscellaneous options**

[`assign_chainID`](#mmpbsa_ifv_assign_chainID){#mmpbsa_ifv_assign_chainID} (Default = 0) 
:   Defines the chains ID assignment mode. _It is ignored when defining a reference structure
(recommended)_. If `assign_chainID = 1`, **xBFreE** check if the structure has no chains ID, and it is assigned 
according to the structure[^1]. If `assign_chainID = 2`, **xBFreE** assign the chains ID, exist or not, 
according to the structure[^1] (can generate inconsistencies).

  [^1]: _The chain ID is assigned according to two criteria: **terminal amino acids** and **residue numbering**. If
        both criteria or residue numbering changes are present, we assign a new chain ID. If there are terminal 
        amino acids, but the numbering of the residue continues, we do not change the ID of the chain._

[`exp_ki`](#mmpbsa_ifv_exp_ki){#mmpbsa_ifv_exp_ki} (Default = 0.0)
:   Specify the experimental Ki (in nM) for correlations analysis. If not defined or exp_ki = 0 then this system 
will be omitted in the correlation analysis

[`full_traj`](#mmpbsa_ifv_full_traj){#mmpbsa_ifv_full_traj} (Default = 0)
:   Print trajectories

    * 0: Print only thread trajectories in *.mdcrd format
    * 1: Print a full traj and the thread trajectories in *.mdcrd format

[`exe_path`](#mmpbsa_ifv_exe_path){#mmpbsa_ifv_exe_path} 
:   Define a list of path to search for executables (gromacs, namd, delphi, etc.). This path takes precedence over the 
paths defined in the PATH variable. Please,

    !!! note "Keep in mind"
        
        * Note that if this variable is not defined, the necessary executables will be searched in the `PATH`.
        * By defining this variable you can use other versions of the same program that are in other paths than the `PATH`

[`keep_files`](#mmpbsa_ifv_keep_files){#mmpbsa_ifv_keep_files} (Default = 1)
:   Defines if temporary files will be deleted or not.

    * 0: Remove all temporary files 
    * 1: Keep all temporary files
    
    !!! note "Keep in mind"
        * Please note that temporary files may be required for compatibility with higher versions.
        * If you remove the temporary files you won't be able to do --rewrite-output to change some aspects of the 
        output like verbose
    
    ??? gmx-mmpbsa "For gmx_MMPBSA users!"
        Since we have improved the workflow to be more organized, this variable is different in gmx_MMPBSA. In this 
        case, the binary file becomes the output for the default analysis.
    

[`netcdf`](#mmpbsa_ifv_netcdf){#mmpbsa_ifv_netcdf} (Default = 0)
:   Specifies whether to use NetCDF trajectories internally rather than writing temporary ASCII trajectory
files. For very large trajectories, this could offer significant speedups, and requires less temporary space. 
However, this option is incompatible with alanine scanning.

    * 0: Do NOT use temporary NetCDF trajectories
    * 1: Use temporary NetCDF trajectories

[`process_trajectory`](#mmpbsa_ifv_process_trajectory){#mmpbsa_ifv_process_trajectory} (Default = 1)
:   Define if it is necessary to generate a clean trajectory, including remove water, ions or select molecules.
    
    * 0: Don’t
    * 1: Generate clean trajectory

    ??? gmx-mmpbsa "For gmx_MMPBSA users!"
        Replace `solvated_trajecotry`.


[`verbose`](#mmpbsa_ifv_verbose){#mmpbsa_ifv_verbose} (Default = 1)
:   Specifies how much output is printed in the output file.

    * 0: Print only difference terms
    * 1: Print all complex, receptor, ligand, and difference terms