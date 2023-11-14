## **`&gbnsr6` namelist variables**

!!! note "Keep in mind"
    * GBNSR6 is an implementation of the Generalized Born (GB) model in which the effective Born radii are computed 
    numerically, via the so-called "R6" integration ([ref.][222]) over molecular surface of the solute. In contrast to 
    most GB practical models, GBNSR6 model is parameter free in the same sense as the numerical PB framework is. Thus, 
    accuracy of GBNSR6 relative to the PB standard is virtually unaffected by the choice of input atomic radii. Check
    Chapter [§5](https://ambermd.org/doc12/Amber21.pdf#chapter.5) in Amber manual for a more thorough description 
    of the GBNSR6 model and its parameters.
    * A default GBNSR6 input file can be created as follows:

        === "GROMACS"
            ```
            xbfree gmx_MMPBSA --create_input gbnsr6
            ```
        === "AMBER"
            ```
            xbfree amber_MMPBSA --create_input gbnsr6
            ```
        === "NAMD"
            ```
            xbfree namd_MMPBSA --create_input gbnsr6
            ```
        === "CHARMM"
            ```
            xbfree charmm_MMPBSA --create_input gbnsr6
            ```

    * A sample GBNSR6 input file is shown [here](input_file.md#gbnsr6)
    * A tutorial on binding free energy calculation with GBNSR6 model is available 
    [here](examples/GBNSR6/README.md)

  [222]: https://pubs.acs.org/doi/abs/10.1021/ct200786m

### **Basic input options**

[`epsin`](#mmpbsa_ifv_epsin){#mmpbsa_ifv_epsin} (Default = 1.0)
:   Dielectric constant of the solute region.

[`epsout`](#mmpbsa_ifv_epsout){#mmpbsa_ifv_epsout} (Default = 78.5)
:   Implicit solvent dielectric constant for the solvent.

[`istrng`](#mmpbsa_ifv_istrng){#mmpbsa_ifv_istrng} (Default = 0.0)
:   Ionic strength in M for the GBNSR6 equation.
                           
[`dprob`](#mmpbsa_ifv_dprob){#mmpbsa_ifv_dprob} (Default = 1.4)
:   Radius of the solvent probe.

[`cavity_surften`](#mmpbsa_ifv_cavity_surften){#mmpbsa_ifv_cavity_surften} (Default = 0.005)
:   Surface tension parameter for nonpolar solvation calculation.

### **Options to select numerical procedures**

[`space`](#mmpbsa_ifv_space){#mmpbsa_ifv_space} (Default = 0.5)
:   Sets the grid spacing that determines the resolution of the solute molecular surface. 

    !!! note "Keep in mind"
        Note that memory footprint of this grid-based implementation of GBNSR6 may become large for large structures,
        e.g. the nucleosome (about 25,000 atoms) will take close to 2 GB of RAM when the default grid spacing is 
        used. For very large structures, one may consider increasing the value of space, which will reduce the 
        memory footprint and execution time; however, the accuracy will also decrease.

[`arcres`](#mmpbsa_ifv_arcres){#mmpbsa_ifv_arcres} (Default = 0.2)
:   Arc resolution used for numerical integration over molecular surface.

[`b`](#mmpbsa_ifv_b){#mmpbsa_ifv_b} (Default = 0.028)
:   Specifies the value of uniform offset to the (inverse) effective radii, the default value 0.028 gives 
better agreement with the PB model, regardless of the structure size. For best agreement with the explicit solvent 
(TIP3P) solvation energies, optimal value of B depends on the structure size: for small molecules (number of atoms 
less than 50), B=0 is recommended. With -chagb option, B is calculated automatically based on the solute size.

[`alpb`](#mmpbsa_ifv_alpb-1){#mmpbsa_ifv_alpb-1} (Default = 1)
:   Specifies if ALBP correction is to be used.

    * 0: Canonical GB is used.
    * 1: ALPB is used (default)

### **Options for CHAGB model**

[`chagb`](#mmpbsa_ifv_chagb){#mmpbsa_ifv_chagb} (Default = 0)
:   Define if CHAGB is used.

    * 0: Do not use CHAGB.
    * 1: Use CHAGB.

[`rs`](#mmpbsa_ifv_rs){#mmpbsa_ifv_rs} (Default = 0.52)
:   Dielectric boundary shift compared to the molecular surface.

[`radiopt`](#mmpbsa_ifv_radiopt){#mmpbsa_ifv_radiopt} (Default = 0)
:   Set of intrinsic atomic radii to be used.

    * 0: uses hardcoded intrisic radii optimized for small drug like molecules, and single amino acid
    dipeptides ([ref.][215])
    * 1: intrinsic radii are read from the topology file. Note that the dielectric surface defined using
    these radii is then shifted outwards by Rs relative to the molecular surface. The option is not
    recommended unless you are planning to re-optimize the input radii set for your problem.

  [215]: https://pubs.acs.org/doi/full/10.1021/ct4010917

[`roh`](#mmpbsa_ifv_roh){#mmpbsa_ifv_roh} (Default = 1)
:   Sets the value of R<sup>z</sup><sub>OH</sub> for CHAGB model, the default is 0.586Å. This parameter defines which 
explicit water model is being mimicked with respect to its propensity to cause charge hydration asymmetry. A perfectly 
tetrahedral water , which can not cause charge hydration asymmetry, would have R<sup>z</sup><sub>OH</sub> = 0. The 
options for `roh` are:

    * 1: R<sup>z</sup><sub>OH</sub> = 0.586Å corresponds to TIP3P and SPC/E. 
    * 2: R<sup>z</sup><sub>OH</sub> = 0.699Å for OPC.
    * 3: R<sup>z</sup><sub>OH</sub> = 0.734Å for TIP4P 
    * 4: R<sup>z</sup><sub>OH</sub> = 0.183Å for TIP5P/E. 

[`tau`](#mmpbsa_ifv_tau){#mmpbsa_ifv_tau} (Default = 1.47)
:   Value of τ in the CHAGB model. This dimensionless parameter controls the effective range of the neighboring 
charges (_j_) affecting the CHA of atom (_i_), see ([ref.][215]) for details.
