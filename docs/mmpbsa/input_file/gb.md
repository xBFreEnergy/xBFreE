## **`&gb` namelist variables**

!!! note "Keep in mind"
    * A default GB input file can be created as follows:

        === "GROMACS"
            ```
            xbfree gmx_MMPBSA --create_input gb
            ```
        === "AMBER"
            ```
            xbfree amber_MMPBSA --create_input gb
            ```
        === "NAMD"
            ```
            xbfree namd_MMPBSA --create_input gb
            ```
        === "CHARMM"
            ```
            xbfree charmm_MMPBSA --create_input gb
            ```
            

    * A sample GB input file is shown [here](input_file.md#gb)
    * A tutorial on binding free energy calculation with GB model is available 
    [here](examples/Protein_ligand/ST/README.md)

### **Basic input options**

[`igb`](#mmpbsa_ifv_igb){#mmpbsa_ifv_igb} (Default = 5)
:   Generalized Born method to use (see [§4](https://ambermd.org/doc12/Amber21.pdf#chapter.4) for more info).

    * 1: The [Hawkins, Cramer, Truhlar][191] pairwise GB model (GB-HCT)
    * 2: Modified GB model 1 developed by [A. Onufriev, D. Bashford and D.A. Case][188] (GB-OBC1)
    * 5: Modified GB model 2 developed by [A. Onufriev, D. Bashford and D.A. Case][188] (GB-OBC2)
    * 7: GBn model described by [Mongan, Simmerling, McCammon, Case and Onufriev][206] (GB-Neck)
    * 8: Same GB functional form as the GBn model (igb=7), but with different parameters. Developed by [Nguyen, Pérez, 
         Bermeo, and Simmerling][200] (GB-Neck2)

  [191]: https://pubs.acs.org/doi/10.1021/jp961710n
  [188]: https://onlinelibrary.wiley.com/doi/10.1002/prot.20033
  [206]: https://pubs.acs.org/doi/10.1021/ct600085e
  [200]: https://pubs.acs.org/doi/10.1021/acs.jctc.5b00271

[`alpb`](#mmpbsa_ifv_alpb){#mmpbsa_ifv_alpb} (Default = 0)
:   Use [Analytical Linearized Poisson-Boltzmann (ALPB)][209] approximation to handle electrostatic interactions 
within the implicit solvent model (see [§4.2](https://ambermd.org/doc12/Amber21.pdf#section.4.2)):

    $$
    ∆𝐺_{el} \approx ∆𝐺_{alpb} = -\frac{1}{2} (\frac{1}{ε_{in}} - \frac{1}{ε_{ex}})\frac{1}{1+αβ} \sum_{ij} q_{i}q_{j}(\frac{1}{f_{GB}} + \frac{αβ}{A})
    $$

    where $β = \frac{ε_{in}}{ε_{ex}}$ is the ratio of the internal and external dielectrics, $α=0.571412$, and A 
    is the so-called effective electrostatic size of the molecule (see `arad_method` below). The ALPB requires one 
    of the analytical GB models to be set, that is igb = 1, 2, 5, or 7, for computing the effective Born radii. It uses 
    the same sets of radii as required by the particular GB model.

    * 0: Don't
    * 1: Use ALPB

  [209]: https://aip.scitation.org/doi/10.1063/1.1857811

[`arad_method`](#mmpbsa_ifv_arad_method){#mmpbsa_ifv_arad_method} (Default = 1)
:   Method used to estimate the effective electrostatic size/radius (`A` in ALPB equation) of the molecule 
(See [Sigalov, Fenley, and Onufriev](https://aip.scitation.org/doi/10.1063/1.2177251)).

    * 1: Use structural invariants
    * 2: Use elementary functions
    * 3: Use elliptic integral (numerical)

[`intdiel`](#mmpbsa_ifv_intdiel){#mmpbsa_ifv_intdiel} (Default = 1.0)
:   Define Internal dielectric constant.

[`extdiel`](#mmpbsa_ifv_extdiel){#mmpbsa_ifv_extdiel} (Default = 78.5)
:   Define External dielectric constant.

[`saltcon`](#mmpbsa_ifv_saltcon){#mmpbsa_ifv_saltcon} (Default = 0.0)
:   Salt concentration in Molarity (M).

[`rgbmax`](#mmpbsa_ifv_rgbmax){#mmpbsa_ifv_rgbmax} (Default = 999.0)
:   Distance cutoff in Å to use when computing effective GB radii.

[`surften`](#mmpbsa_ifv_surften){#mmpbsa_ifv_surften} (Default = 0.0072)
:   Surface tension value. Units in kcal/mol/Å^2^

[`surfoff`](#mmpbsa_ifv_surfoff){#mmpbsa_ifv_surfoff} (Default = 0.0)
:   Offset to correct (by addition) the value of the non-polar contribution to the solvation free energy term.

[`molsurf`](#mmpbsa_ifv_molsurf){#mmpbsa_ifv_molsurf} (Default = 0)
:   Define the algorithm to calculate the surface area for the non-polar solvation term.
    
    * 0: LCPO (Linear Combination of Pairwise Overlaps)
    * 1: molsurf algorithm

[`msoffset`](#mmpbsa_ifv_msoffset){#mmpbsa_ifv_msoffset} (Default = 0) 
:   Offset to apply to the individual atomic radii in the system when calculating the `molsurf` surface. See the
description of the `molsurf` action command in [cpptraj][4].

[`probe`](#mmpbsa_ifv_probe){#mmpbsa_ifv_probe} (Default = 1.4)
:   Radius in Å of the probe molecule (supposed to be the size of a solvent molecule), to use when determining the 
molecular surface.
    
    !!! note
        only applicable when `molsurf` is set to 1

### **QM options**

[`ifqnt`](#mmpbsa_ifv_ifqnt){#mmpbsa_ifv_ifqnt} (Default = 0)
:   Specifies whether a part of the system is treated with quantum mechanics.
    
    * 0: Potential function is strictly classical
    * 1: Use QM/MM

    !!! note "Keep in mind"
        * Calculations where part of the system is treated with quantum mechanics can be performed only with GB
        * A sample QM/MMGBSA input file is shown [here](input_file.md#qmmmgbsa)
        * A tutorial on binding free energy calculation with QM/MMGBSA is available 
        [here](examples/QM_MMGBSA/README.md)

[`qm_theory`](#mmpbsa_ifv_qm_theory){#mmpbsa_ifv_qm_theory} 
:   Which semi-empirical Hamiltonian should be used for the quantum calculation. Options are `PM3`, `AM1`, `MNDO`, 
`PDDG-PM3`, `PM3PDDG`, `PDDG-MNDO`, `PDDGMNDO`, `PM3-CARB1`, `PM3CARB1`, `DFTB`, `SCC-DFTB`, `RM1`, `PM6`, 
`PM3-ZnB`, `PM3-MAIS`, `PM3ZNB`, `MNDO/D`, `MNDOD`. The dispersion correction can be switched on for `AM1` 
and `PM6` by choosing `AM1-D*` and `PM6-D`, respectively. The dispersion and hydrogen bond correction will be 
applied for `AM1-DH+` and `PM6-DH+`.

    !!! danger
         No `qm_theory` default, this must be specified if `ifqnt` = 1.

[`qm_residues`](#mmpbsa_ifv_qm_residues){#mmpbsa_ifv_qm_residues}
:   Complex residues to treat with quantum mechanics. All residues treated with quantum mechanics in the complex 
must be treated with quantum mechanics in the receptor or ligand to obtain meaningful results. This notation is 
the same used for `print_res` variable in `&decomp` namelist.

    !!! danger
         No `qm_residues` default, this must be specified if `ifqnt` = 1.

    !!! example "Selection schemes"

        === "By Distance (recommended)"
            Notation: [ `within` `distance` ]
            :   `within` corresponds to the keyword and `distance` to the maximum distance criterion in Å necessary to 
                select the residues from both the receptor and the ligand. In case the cutoff used is so small that 
                the number of qm_residues = 0, the cutoff value will be increased by 0.1 until the number of 
                qm_residues > 0.
    
            !!! example
                `qm_residues="within 5"` Residues within 5 Å between receptor and ligand will be treated with quantum 
                mechanic.

        === "Amino acid selection"
            Notation: [ `CHAIN`/(`RESNUM` or `RESNUM-RESNUM`) ]
            :    Treat with quantum mechanics residues individual or ranges. This notation also supports insertion 
            codes, in which case you must define them individually

            `qm_residues="A/1,3-10,15,100"` This treat with quantum mechanic Chain A residues 1, 3 through 10, 15, and 
            100 from the complex topology file and the corresponding residues in either the ligand and/or receptor 
            topology files.
    
            Let's suppose that we can have the following sequence: - A:LEU:5 - A:GLY:6:A - A:THR:6:B - A:SER:6:C - 
            A:ASP:6:D - A:ILE:7
    
            with the format `CHAIN`/`RESNUMBER` `INSERTION_CODE`
            
            === "Right notation"
                
                **Ranges selection**
                :   `qm_residues="A/5-7` Will treat with quantum mechanic all mentioned residues because all residues with 
                insertion code are contained in the range
                
                **Individual selection**
                :   `qm_residues="A/5,6B,6C,7` Will treat with quantum mechanic all mentioned residues except the 
                residues 6A and 6D from chain A
                
                **Multiple chain selection**
                :   `qm_residues="A/5-10,100 B/34,56` Will treat with quantum mechanic residues 5 through 10, and 100 from 
                chain A, and residues 34 and 56 from Chain B.
    
            === "Wrong notation"
                `qm_residues="A/5-6B,6D-7` Will end in error.


[`qmcut`](#mmpbsa_ifv_qmcut){#mmpbsa_ifv_qmcut} (Default = 9999.0)
:   The cutoff for the qm/mm charge interactions.

[`scfconv`](#mmpbsa_ifv_scfconv){#mmpbsa_ifv_scfconv} (Default = 1.0e-8)
:   Controls the convergence criteria for the SCF calculation, in kcal/mol. The tighter the 
convergence the longer the calculation will take. Values tighter than 1.0e-11 are not recommended as these can lead 
to oscillations in the SCF, due to limitations in machine precision, that can lead to convergence failures.

[`writepdb`](#mmpbsa_ifv_writepdb){#mmpbsa_ifv_writepdb} (Default = 1)
:   Write a PDB file of the selected QM region. This option is designed to act as an aid to the user to
allow easy checking of what atoms were included in the QM region. Write a PDB file of the atoms in the QM region 
on the very first step to a file named qmmm_region.pdb.

    * 0: Don't
    * 1: Write a PDB file of the selected QM region

[`peptide_corr`](#mmpbsa_ifv_peptide_corr){#mmpbsa_ifv_peptide_corr} (Default = 0)
:   Apply MM correction to peptide linkages. This correction is of the form: 

    $$
    E_{scf} = E_{scf} + h_{type}(i_{type}) * sin^{2}\phi
    $$

    where _ϕ_ is the dihedral angle of the H-N-C-O linkage and $h_{type}$ is a constant dependent on the 
    Hamiltonian used. Recommended, except for DFTB/SCC-DFTB.

    * 0: Don't
    * 1: Apply a MM correction to peptide linkages

[`verbosity`](#mmpbsa_ifv_verbosity){#mmpbsa_ifv_verbosity} (Default = 0)
:   Controls the verbosity of QM/MM related output. Values of 2 or higher will produce a lot of output.

    * 0: only minimal information is printed - Initial QM geometry and link atom positions as
    well as the SCF energy at every ntpr steps.
    * 1: Print SCF energy at every step to many more significant figures than usual. Also print the
    number of SCF cycles needed on each step.
    * 2: As 1 and also print info about memory reallocations, number of pairs per QM atom, QM core -
    QM core energy, QM core - MM atom energy, and total energy.
    * 3: As 2 and also print SCF convergence information at every step.
    * 4: As 3 and also print forces on the QM atoms due to the SCF calculation and the coordinates of
    the link atoms at every step.
    * 5: As 4 and also print all of the info in kJ/mol as well as kcal/mol.
    
  [4]: https://ambermd.org/doc12/Amber21.pdf#subsection.34.11.49

## Sample input files

!!! tip
    You can refer to the [examples](examples/README.md) to understand the input file in a practical way.

!!! warning
    These are illustrative examples, please, don't use it for production. Create a new one using 
    the instructions provides above in the section [Generation of input files with **xBFreE**](#generation-of-input-files-with-xbfree)

``` linenums="1"
Sample input file for GB calculation building the Amber topologies
from structures. Please refer to the section "How gmx_MMPBSA works"

&general
startframe=5, endframe=100, interval=5, verbose=2, 
/

&gb
igb=5, saltcon=0.150,
/
```