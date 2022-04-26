abYpack V1.0
============

(c) 2022, UCL, Prof. Andrew C.R. Martin
---------------------------------------

abYpack takes a specified VH/VL packing angle and a PDB file and packs
the domains together at the specified angle (as calculated by
`abpackingangle`, Abhinandan, K.R. and Martin, A.C.R. (2010) Analysis
and prediction of VH/VL packing in antibodies. *Protein Engineering
Design and Selection* **23**:689-697. [doi:10.1093/protein/gzq043] [PMID:
20591902]).

### Installation

The following symbolic links must be made:

- `bin` - the abYmod binary directory containing `abpackinganghle`,
   `fitlhpdb`, and `pdbgetchain`
- `DATA` - the abYmod DATA directory containing `abpdblib` (library of
   antibody PDB files) and `angles.pa`, (a two-column list of packing
   angles and filenames in `abpdblib`)

If these links are created, then `config.pm`, is just a symbolic link
to `config.pm.dist`. Alternatively the above data/programs can be
placed elsewhere and the config file adapted accordingly.

