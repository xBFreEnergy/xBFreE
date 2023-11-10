## **`&nmode` namelist variables**

!!! note "Keep in mind"
    * A default nmode input file can be created as follows:

        === "GROMACS"
            ```
            xbfree gmx_MMPBSA --create_input nmode
            ```
        === "AMBER"
            ```
            xbfree amber_MMPBSA --create_input nmode
            ```
        === "NAMD"
            ```
            xbfree namd_MMPBSA --create_input nmode
            ```
        === "CHARMM"
            ```
            xbfree charmm_MMPBSA --create_input nmode
            ```
    
    * A sample nmode input file is shown [here](input_file.md#entropy-with-nmode)
    * A tutorial on normal mode analysis is available [here](examples/Entropy_calculations/nmode/README.md)

### **Basic input options**

[`nmstartframe`](#mmpbsa_ifv_nmstartframe){#mmpbsa_ifv_nmstartframe}
:   Frame number to begin performing `nmode` calculations on 

    !!! note  
        This variable will choose a subset of the frames chosen from the variables in the `&general` namelist. Thus,
        the "trajectory" from which snapshots will be chosen for `nmode` calculations will be the collection of 
        snapshots upon which the other calculations were performed.

[`nmendframe`](#mmpbsa_ifv_nmendframe){#mmpbsa_ifv_nmendframe} (Default = 1000000)
:   Frame number to stop performing `nmode` calculations on 

    !!! note  
        This variable will choose a subset of the frames chosen from the variables in the `&general` namelist. Thus,
        the "trajectory" from which snapshots will be chosen for `nmode` calculations will be the collection of 
        snapshots upon which the other calculations were performed.

[`nminterval`](#mmpbsa_ifv_nminterval){#mmpbsa_ifv_nminterval} (Default = 1)
:   Offset from which to choose frames to perform `nmode` calculations on

    !!! note  
        This variable will choose a subset of the frames chosen from the variables in the `&general` namelist. Thus,
        the "trajectory" from which snapshots will be chosen for `nmode` calculations will be the collection of 
        snapshots upon which the other calculations were performed.

### **Parameter options**

[`nmode_igb`](#mmpbsa_ifv_nmode_igb){#mmpbsa_ifv_nmode_igb} (Default = 1)
:   Value for Generalized Born model to be used in calculations. Options are:
    
    * 0: Vacuum
    * 1: HCT GB model 

[`nmode_istrng`](#mmpbsa_ifv_nmode_istrng){#mmpbsa_ifv_nmode_istrng} (Default = 0.0)
:   Ionic strength to use in `nmode` calculations. Units are Molarity (M). Non-zero values are ignored if `nmode_igb`
is 0 above.

[`dielc`](#mmpbsa_ifv_dielc){#mmpbsa_ifv_dielc} (Default = 1.0)
:   Distance-dependent dielectric constant 

[`drms`](#mmpbsa_ifv_drms){#mmpbsa_ifv_drms} (Default = 0.001)
:   Convergence criteria for minimized energy gradient.

[`maxcyc`](#mmpbsa_ifv_maxcyc){#mmpbsa_ifv_maxcyc} (Default = 10000)
:   Maximum number of minimization cycles to use per snapshot in sander.
