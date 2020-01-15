# lPDBmngr

local PDB manager is a python tool I wrote to deal with a local PDB copy files location and computed derivative data. The main motivation was to facilitate how to deal with writing and loading data computed for PDB entry files on different data mining projects. Everything is centred on the _lPDB class_, defined at 'lPDB.py'. 


PS: I wrote this with my own selfish desires on mind, but if you find anything useful feel free to use it. =)

---

## HOW TO?

### Generate metadata '.csv'

To obtain a '.csv' file with file paths, bc groups associated and pre-computed data[1] for each entry, edit the

```{bash}
>$ python3.6 mount_lpdb_mtdt.py /path/to/lPDB_dir/
```

### Compute derivative data

To run some software on entries of a local PDB copy, you can use one of the methods available[1].

```{bash}
>$ python3.6 compute_drvtv_data /path/to/lPDB_dir/ 4
```
the second argument sets the number of cpus to use.


[1] Currently, only DSSP and FleXgeo are supported.

---
