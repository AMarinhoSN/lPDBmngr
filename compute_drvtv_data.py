import lPDB
import pickle as pck
import subprocess

# set input variables #
source_dir = '/path/to/lPDB_dir/'
ncpus = 10
#######################
'''
        >> compute_drvtv_data.py <<

This script generate a lPDB object load coordinate files location and metadata
for chain files, compute derivative data[1] and associate to the inputs. This
metadata of the local pdb copy is stored as a pandas dataframe and writen as
.csv file as well.

[1] Currently, DSSP and FleXgeo.
'''
print('@ generating min metadata csv...')
# 1 - mount minimal metadata
lPDB_obj = lPDB.lPDB(source_dir)
# add chain files path to metadata
lPDB_obj.add_chains_coordfls()
# export metadata to csv
lPDB_obj.export_to_csv()

# 2 - compute derivative data
print('@ computing derivative data...')
print('  -> dssp...')
lPDB_obj.run_dssp(ncpus=ncpus)
print('  -> flexgeo...')
lPDB_obj.run_flexgeo(ncpus=ncpus)

# 3 - add info to lPDB metadata dataframe
print('@ adding files path to metadata dataframe...')
lPDB_obj.add_xgeo_fls()
lPDB_obj.add_dssp_fls()
lPDB_obj.export_to_csv()
print(':: DONE ::')
