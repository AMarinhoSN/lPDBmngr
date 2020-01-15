import lPDB
from sys import argv
'''
This script generates a dataframe with local PDB copy metadata.
Currently:
    1 - get pdb files location;
    2 - associate each entry with different bc groups
    3 - dssp out files
    4 - xgeo files
    5 - export metadata as a csv file
    6 - save lPDB object state as a pickle file.

USAGE: python3.6 mount_lpdb_mtdt.py /path/to/lPDB_dir/
'''
# set lPDB source dir
src_dir = argv[1]

print('@ get entries info...')
lPDB_obj = lPDB.lPDB(src_dir)
print('@ get pdbid chains location...')
lPDB_obj.add_chains_coordfls()
print('@ add bc groups columns...')
lPDB_obj.add_bcgrps_col()
print('@ add dssp filepaths...')
lPDB_obj.add_dssp_fls()
print('@ add xgeo filepaths...')
lPDB_obj.add_xgeo_fls()
print('@ saving csv file...')
lPDB_obj.export_to_csv()
print('@ saving lPDB object state...')
lPDB_obj.save_state()
print(':: DONE ::')
