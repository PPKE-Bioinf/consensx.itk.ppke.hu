# CoNSEnsX

Source code of the [consensx.itk.ppke.hu](http://consensx.itk.ppke.hu) site.

The server is specifically designed for structural ensembles that were generated to reflect the internal dynamics of proteins at a given time scale. Such ensembles are typically generated by restrained molecular dynamics methods or selection-based approaches. CoNSEnsX+ provides an ensemble-averaged analysis of all experimental parameters recognized in the input and offers a simple greedy selection approach to identify the sub-ensemble best reflecting the parameters chosen.

For detailed usage and description visit the [CoNSEnsX description page](http://consensx.itk.ppke.hu/static/Descriptionpage.html).


## NOE restraints

The following STR fields are needed for CoNSEnsX to identify NOE data:
* in the `save` block:
```
    Gen_dist_constraint_list.Constraint_type      NOE
```

* in the `loop` block:
```
    Gen_dist_constraint.ID
    Gen_dist_constraint.PDB_residue_no_1
    Gen_dist_constraint.PDB_residue_no_2
    Gen_dist_constraint.Comp_ID_1
    Gen_dist_constraint.Comp_ID_2
    Gen_dist_constraint.Atom_ID_1
    Gen_dist_constraint.Atom_ID_2
```

> **Note:** The atomic naming convention must match the naming convention used in the input PDB file

----------


## Back-calculated parameters


#### RDC

The following STR fields are needed for CoNSEnsX to identify RDC data:
* in the `save` block:
```
    RDC_list_[num]
```
* in the `loop` block:
```
    RDC.Seq_ID_1        OR Atom_one_residue_seq_code
    RDC.Atom_type_1     OR Atom_one_atom_name
    RDC.Seq_ID_2        OR Atom_two_residue_seq_code
    RDC.Atom_type_2     OR Atom_two_atom_name
    RDC.Val             OR Residual_dipolar_coupling_value
```

#### Order parameters

The following STR fields are needed for CoNSEnsX to identify order parameter data:
* in the `save` block:
```
    order_param
```
* in the `loop` block:
```
    S2_value
    Residue_seq_code
    Atom_name
```

#### Scalar coupling

The following STR fields are needed for CoNSEnsX to identify scalar coupling data:
* in the `save` block:
```
    coupling_constants
```
* in the `loop` block:
```
    Atom_one_residue_seq_code
    Coupling_constant_code
```

#### Chemical shifts

The following STR fields are needed for CoNSEnsX to identify chemical shifts data:
* in the `save` block:
```
    chem_shift_list_[num] OR assigned_chem_shift_list_[num]
```
* in the `loop` block:
```
    Atom_chem_shift.Seq_ID        OR Residue_seq_code
    Atom_chem_shift.Comp_ID     OR Residue_label
    RDC.Seq_ID_2        OR Atom_name
    RDC.Atom_type_2     OR Chem_shift_value
```
