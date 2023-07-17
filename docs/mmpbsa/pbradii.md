
!!! warning
    Please, note that this section is in active development; thus some elements could change in future versions 

## Introduction

Previous results suggest that the PB/GB predictions are largely impacted by the radii used, as the multiple variants of 
radii give very different results (check [Wang et. al., 2019][14] for an extensive review on this and other factors). 
In our interest in providing the users with all available tools, we have implemented 
10 radii sets available in the literature that you can define with the [PBRadii](input_file.md#mmpbsa_ifv_PBRadii) 
variable in the [input file](input_file.md):

|   Radii Name   |                     Type                     |             FF Compatibility             | Compatibility |
|:--------------:|:--------------------------------------------:|:----------------------------------------:|:-------------:|
|    `amber6`    |      `atom_number`+`bonded_atom_number`      |       any, but optimized for Amber       |    generic    |
|    `bondi`     |                `atom_number`                 |       any, but optimized for Amber       |    generic    |
|    `mbondi`    |      `atom_number`+`bonded_atom_number`      |       any, but optimized for Amber       |    generic    |
|   `mbondi2`    |      `atom_number`+`bonded_atom_number`      |       any, but optimized for Amber       |    generic    |
|   `mbondi3`    | `atom_number`+`bonded_atom_number`+`residue` |       any, but optimized for Amber       |    generic    |
|  `mbondi_pb2`  |      `atom_number`+`bonded_atom_number`      |       any, but optimized for Amber       |    generic    |
|  `mbondi_pb3`  |      `atom_number`+`bonded_atom_number`      |       any, but optimized for Amber       |    generic    |
| `charmm_radii` |                `charmm_radii`                |                  charmm                  |   complete    |
|    `tyl06`     |                  `residue`                   | amber (optimized), charmm (experimental) | semi-complete |
|  `yamagishi`   |                  `residue`                   | amber (optimized), charmm (experimental) | semi-complete |

  [14]: https://pubs.acs.org/doi/abs/10.1021/acs.chemrev.9b00055

## Radii file structure

The radii data is contained in a JSON format file. It is a dictionary (key-value pair) that contains 
all the required information. The most important keys are `forcefield`, `type`, and `radii`, so please, pay attention 
to them when modifying or generating a new radii set. This format is extremely versatile and allows adding new atoms or 
modifying the current ones:

```json  title="mbondi radii"
{
  "name": "mbondi",
  "forcefield": ["all"],   // (1)!
  "implementation": "xBFreE Team",
  "RADIUS_SET": "modified Bondi radii (mbondi)",
  "reviews": [],
  "reference": "",
  "notes": "",
  "type": "atom_number+bonded_atom_number",   // (2)!
  "radii": [   // (3)!
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 6 , "mass": null, "radii": 1.3 },
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 7 , "mass": null, "radii": 1.3 },
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 8 , "mass": null, "radii": 0.8 },
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 16, "mass": null, "radii": 0.8 },
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.2 },
    {"atom_number": 6 , "atom_type": "^C1", "bonded_atom_number": null, "mass": ">13.0", "radii": 2.2 },
    {"atom_number": 6 , "atom_type": "^C2", "bonded_atom_number": null, "mass": ">14.0", "radii": 2.2 },
    {"atom_number": 6 , "atom_type": "^C3", "bonded_atom_number": null, "mass": ">15.0", "radii": 2.2 },
    {"atom_number": 6 , "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.7 },
    {"atom_number": 7 , "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.55},
    {"atom_number": 8 , "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.5 },
    {"atom_number": 9 , "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.5 },
    {"atom_number": 14, "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 2.1 },
    {"atom_number": 15, "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.85},
    {"atom_number": 16, "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.8 },
    {"atom_number": 17, "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.7 },
    {"atom_number": "*", "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.5 }
  ]
}
```

1. Please, pay attention to this key when modifying or generating new radii set
2. Please, pay attention to this key when modifying or generating new radii set
3. Please, pay attention to this key when modifying or generating new radii set

The content of the key radii corresponds to the conditions for assigning the value of radii per atom. As seen in the 
Type column of the table above, there are various combinations depending on the selected radii set. 

## Understanding the radii classification/type

In the table above, the Type column shows several combinations. Each keyword is used as a condition to assign the 
radii to any atom. Let's analyze the previous example: `mbondi`
In this case, a combination of `atom_number` and `bonded_atom_number` is highlighted. 

??? note "Why not include `atom_type` or `mass`?"
    The `atom_type` keyword is included in all generic defined radii, so it would be redundant to define it. On the 
    other hand, `mass` is evaluated for all radii set based on `atom_number` (all except `charmm_radii`, `tyl06` and 
    `yamagishi`)

Here we first check the atom_number, then if `atom_type` and `mass` or `bonded_atom_number` were defined. If the 
processed atom meets the conditions, it is assigned the specific defined value, otherwise, it is assigned a generic 
value defined in the file with `atom_number="*"`  

## Simple case to add a new atom parameter
Suppose you want to add a new atom parameter, for example, iron (Fe). For this, you must copy the radii base file 
(for example, `mbondi`) and add the new data as follows:

```json  title="modified mbondi radii"
{
  "name": "new-mbondi",    // (1)!
  "forcefield": ["all"],   
  "implementation": "Me",  // (2)!
  "RADIUS_SET": "modified Bondi radii (mbondi) for Fe",
  "reviews": [],
  "reference": "",
  "notes": "I have added the Fe parameters",    // (3)!
  "type": "atom_number+bonded_atom_number",   // (4)!
  "radii": [   
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 6 , "mass": null, "radii": 1.3 },
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 7 , "mass": null, "radii": 1.3 },
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 8 , "mass": null, "radii": 0.8 },
    {"atom_number": 1 , "atom_type": null, "bonded_atom_number": 16, "mass": null, "radii": 0.8 },
    ...
    {"atom_number": 17, "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.7 },
    {"atom_number": 26, "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.6 }, // (5)!
    {"atom_number": "*", "atom_type": null, "bonded_atom_number": null, "mass": null, "radii": 1.5 }
  ]
}
```

1. New name
2. author
3. Author's notes
4. Please, take care your modifications and the type
5. Note that we add the new atom parameters before the `"atom_number": "*"`

Now, this new file can be defined as the radii in the variable [PBRadii](input_file.md#mmpbsa_ifv_PBRadii), for example:
```yaml
Sample input file for GB calculation

&general
  PBRadii = /home/user/path/to/new_radii.json 
/

&gb
igb=5, saltcon=0.150,
/
```

!!! tip "Create, test, and share your radii"

    * This feature is experimental, so we need to test it in a broader scenario. If you find a bug, please report it.
    You will receive due credit in the "reviews" section of the radii archive.
    * If you create or optimize your radii, feel free to share it, and we will include it in the list of radii.
    * We don't show much information about the radii used at this time, but we plan to have it displayed in the output 
    file.
    