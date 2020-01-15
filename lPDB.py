import pandas as pd
import glob
import multiprocessing as mp
import subprocess
import pickle as pck

'''
Here the main class of local PDB copy is defined and methods to compute and load
information should be added here.
'''

## FUNCTIONS ###################################################################
# ~~ Load metadata files provided on PDB server and return dataframe ~~
def load_resolu_data_df(resolu_data_flpth):
    '''load data from resolu.idx'''
    f = open(resolu_data_flpth, 'r')
    on_data = False
    dct_lst = []
    for line in f:
        if line.startswith('----'):
            on_data = True
            continue
        if on_data == True:
            data_set = line.replace('\t', '').replace('\n','').split(';')
            pdbid = data_set[0].lower()
            if data_set[1] == '':
                res = None
            else:
                res = float(data_set[1])
            dct = {'pdbid': pdbid, 'res': res}
            dct_lst.append(dct)
    return pd.DataFrame(dct_lst)

def load_entry_type_df(entry_type_flpth):
    '''load data from pdb_entry_type.txt'''
    col_nms = ['pdbid', 'mol_type', 'method']
    return pd.read_csv(entry_type_flpth, header=None, names=col_nms,sep='\t')

# ~~ pandas handy functions ~~~
def get_dct(flspth, newcol_prfx):
    '''
    get dictionary to associate pdbid and chain for arbitrary list of files.

    This method assumes first 6 characteres of files names follow the
    <pdbid>_<chain> scheme (e.g. 1r9w_A.dssp).

    Parameters
    ---
    :flspth:<str>, file path
    :newcol_name:<str>, new column name

    Return
    ---
    :<dct>, a dictionary {'pdbid', 'chain', 'colname'}
    '''
    flnm = flspth.split('/')[-1]
    pdbid = flnm[0:4]
    chain = flnm[5]
    return {'pdbid':pdbid, 'chain':chain, newcol_prfx+'_flspth':flspth}

def get_ncpus(ncpus):
    if ncpus == 'All':
        ncpus = mp.cpu_count()
    else:
        assert(type(ncpus)==int), 'ncpus must be int'
        return ncpus

def run_command_on(cmd, colnm):
    ''' '''
    cmd = cmd.split(' ')
    print(cmd)

def load_bc_df(bc_flpth, bc_name):
    '''Parser bc-X0.out files data and load a dataframe'''
    # load file
    bc_f = open(bc_flpth, 'r')

    # each line contains a set of entries which belongs to the same cluster
    dct_lst = []
    c = 1
    for line in bc_f:
        grp_entries = line.split(' ')
        grp_n=c
        for entry in grp_entries[0:20]:
            pdbid = entry.split('_')[0].lower()
            chain = entry.split('_')[1][0]
            dct = {'pdbid':pdbid, 'chain':chain, bc_name+'_grp':grp_n}
            dct_lst.append(dct)
        c+=1
    return pd.DataFrame(dct_lst)

# ~~ command handy functions ~~~
def create_derivative_dir(newdir_nm, lPDB_obj):
    # create dir
    newdir_path = lPDB_obj.source_dir+'pdb/'+newdir_nm+'/'
    subprocess.run(['mkdir', newdir_path])

    # create subdir structure
    subdir_lst = glob.glob(lPDB_obj.source_dir+'pdb/pdb_chain/*')
    for dir in subdir_lst:
        dirnm = dir.split('/')[-1]
        subprocess.run(['mkdir', newdir_path+dirnm])

def get_outflpth(in_f, format, out_dir):
    subdir = in_f.split('/')[-2]+'/'
    flnm = in_f.split('/')[-1].replace('.pdb', format)
    return out_dir+subdir+flnm

# ~~~ misc functions ~~~
def get_fls_mtd_df(search_str, newcol_nm):
    found_fls_lst = glob.glob(search_str)
    assert(len(found_fls_lst) > 0), 'No files found using '+search_str

    # get chains dataframe
    dct_lst = []
    for fl in found_fls_lst:
        dct_lst.append(get_dct(fl, newcol_nm))

    return pd.DataFrame(dct_lst)

def run_shell_cmd_PP(cmd_lst, ncpus):
    '''
    run a list of terminal comands in parallel.
    PS: This function do not return anything.
    '''
    ncpus = get_ncpus(ncpus)
    workers = mp.Pool(processes=ncpus)
    workers.map(subprocess.run, cmd_lst)
    workers.close()

################################################################################

class lPDB:
    '''
    Local PDB class.
    This class aim to provide handy methods to deal with a local copy of the
    PDB.

    USAGE:
        >>> import lPDB
        >>> lPDB = lPDB.lPDB('/path/to/locaPDBcopy_dir/')
    '''

    def __init__(self, source_dir):
        # set source directory
        self.source_dir = source_dir
        self.metadata = None
        # load minimum metadata
        self.load_barebones_mtdt()

    def load_barebones_mtdt(self):
        '''
        load information provided by PDB server.
        Will load data from the following files at './pdb/derived_data/':
        1 - pdb_entry_type.txt
        2 - resol.idx

        and generate a single dataframe with the combined information of those
        files.
        '''
        # set drived data dir path
        derived_dirpth = self.source_dir+'pdb/derived_data/'
        # load entry type data
        entry_type_data_flpth = derived_dirpth+'pdb_entry_type.txt'
        entry_type_df = load_entry_type_df(entry_type_data_flpth)
        # load resol.idx
        resolu_data_flpth = derived_dirpth+'resolu.idx'
        resol_df = load_resolu_data_df(resolu_data_flpth)

        # merge dataframes
        self.metadata = pd.merge(resol_df, entry_type_df, on='pdbid')

    def export_to_csv(self):
        self.metadata.to_csv(self.source_dir+'lPDB_metadata.csv')

    def save_state(self, prfx=''):
        '''save state as pickle file at source dir'''
        pck.dump(self, open(self.source_dir+'lPDB.p', 'wb'))
    # ~~~~~~~~~~ methods to add new columns to metadata ~~~~~~~~~~~~~~~~~~~~~~~~
    def add_chains_coordfls(self):
        # get files location
        newcol_prfx = 'chain'
        chains_srch_str = self.source_dir+'pdb/pdb_chain/*/*.pdb'
        chains_mtdt_df = get_fls_mtd_df(chains_srch_str, newcol_prfx)
        # add to metadata
        self.metadata = pd.merge(self.metadata, chains_mtdt_df, on='pdbid')

    def add_dssp_fls(self, dssp_srch_str=None):
        '''add dssp files path to metadata '''
        # sanity check
        if dssp_srch_str == None:
            dssp_srch_str = self.source_dir+'pdb/pdb_dssp/*/*'
        else:
            assert(type(dssp_srch_str)==str), 'dssp_srch_str must be a string'
        # get files location
        dssp_mtdt_df = get_fls_mtd_df(dssp_srch_str, 'dssp')
        # add to metadata
        self.metadata = pd.merge(self.metadata, dssp_mtdt_df, on=['pdbid', 'chain'])

    def add_xgeo_fls(self, xgeo_srch_str=None):
        '''add xgeo files path to metadata '''
        # sanity check
        if xgeo_srch_str == None:
            xgeo_srch_str = self.source_dir+'pdb/pdb_xgeo/*/*_xgeo.csv'
        else:
            assert(type(xgeo_srch_str)==str), 'xgeo_srch_str must be a string'
        # get files location
        xgeo_mtdt_df = get_fls_mtd_df(xgeo_srch_str, 'xgeo')
        # add to metadata
        self.metadata = pd.merge(self.metadata, xgeo_mtdt_df,
                                 on=['pdbid', 'chain'])

    def add_bcgrps_col(self, bc_srch_str=None):
        '''
        Load clustered pdb data from 'bc-X0.out' files and add bcX_groups
        column to metadata
        '''
        # get list of .out files
        if bc_srch_str == None:
            bc_dir = self.source_dir+'pdb/bc/'
            bc_src_str = bc_dir+'bc-*.out'

        bc_fls_found = glob.glob(bc_src_str)

        # load data and store new columns to be added on self.metadata
        df = self.metadata[['pdbid','chain']]
        for bc_fl in bc_fls_found:
            bc_name = bc_fl.split('/')[-1].replace('.out','')
            bc_df = load_bc_df(bc_fl, bc_name)
            df = df.merge(bc_df, on=['pdbid', 'chain'])
        self.metadata = self.metadata.merge(df, on=['pdbid', 'chain'], how='outer')

    # ~~~~~~~ methods to generate derivative data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def run_dssp(self, coord_col='chain_flspth', ncpus='All'):
        '''run dssp on all entries of  '''
        # create dir structure
        create_derivative_dir('pdb_dssp', self)


        # get list of commands to run
        prot_only = self.metadata['mol_type'] == 'prot'
        input_fls = self.metadata[coord_col].loc[prot_only].values
        format = '.dssp'
        out_dir = self.source_dir+'pdb/pdb_dssp/'
        cmd_lst = []
        for in_fl in input_fls:
            out_fl = get_outflpth(in_fl, format, out_dir)
            cmd_lst.append(['dssp', '-i', in_fl, '-o', out_fl])

        # run commands
        ncpus = get_ncpus(ncpus)
        workers = mp.Pool(processes=ncpus)
        workers.map(subprocess.run, cmd_lst)
        workers.close()

    def run_flexgeo(self, coord_col='chain_flspth', ncpus='All'):
        ''' '''
        ## TEMPORARY SOLUTION##
        flexgeo_str = '/path/to/FleXgeo/bin/FleXgeo_LINUX'
        ################################
        # create dir structure
        create_derivative_dir('pdb_xgeo', self)
        # get list of commands
        cmd_lst =[]
        prot_only = self.metadata['mol_type'] == 'prot'
        # get list of commands
        input_fls = self.metadata[coord_col].loc[prot_only].values
        format = ''
        out_dir = self.source_dir+'pdb/pdb_xgeo/'

        cmd_lst = []
        for in_fl in input_fls:
            out_fl = get_outflpth(in_fl, format, out_dir)
            cmd=[flexgeo_str,'-pdb='+in_fl,'-isSingle','-ncpu=1',
                '-outprfx='+out_fl]
            cmd_lst.append(cmd)

        # run commands
        run_shell_cmd_PP(cmd_lst, ncpus)
        ## WARNING TEMPORARY SOLUTION ##
        # on next flexgeo update, those files should not be written
        ####
        subprocess.run(['/usr/bin/find', out_dir+'*/*_NORM.csv', '-delete'])
        subprocess.run(['/usr/bin/find', out_dir+'*/*_MEAN.csv', '-delete'])
        subprocess.run(['/usr/bin/find', out_dir+'*/*_VAR.csv', '-delete'])
        subprocess.run(['/usr/bin/find', out_dir+'*/*_STD.csv', '-delete'])
        #################################
