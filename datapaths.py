from pathlib import Path

data_directory_name = Path('kids_acquisition_data')
working_directory = Path('/media/federico/DATA/' / data_directory_name)
print('Setting up the local data directory at "{:s}"...'.format(working_directory.as_posix()))
if not working_directory.exists():
    working_directory.mkdir()
    print('Local data directory not found. Created.')
    

target = working_directory / Path('target')
vna = working_directory / Path('vna')
target_processed = working_directory / Path('target_processed')
vna_processed = working_directory / Path('vna_processed')
dirfile = working_directory / Path('data_logger')

remote_username = 'admin'
remote_hostname = 'nasg31.roma1.infn.it'
remote_directory_path = Path('/share/kids_acquisition_data')


def check_if_dirs_exist(mkdir = False):
    '''
    Function that checks if data directories exist.

    Parameters
    ----------
    mkdir : boolean, optional
        If True the data directories will be created. The default is False.

    Returns
    -------
    None.

    '''
    
    error_msg = ' directory does not exist'
    add_dir_msg = ' directory added'
    
    if not target.exists():
        print(target.as_posix() + error_msg)
        if mkdir: 
            target.mkdir()
            print(target.as_posix() + add_dir_msg)
    if not vna.exists():
        print(vna.as_posix() + error_msg)
        if mkdir: 
            vna.mkdir()
            print(vna.as_posix() + add_dir_msg)
    if not target_processed.exists():
        print(target_processed.as_posix() + error_msg)
        if mkdir: 
            target_processed.mkdir()
            print(target_processed.as_posix() + add_dir_msg)
    if not vna_processed.exists():
        print(vna_processed.as_posix() + error_msg)
        if mkdir: 
            vna_processed.mkdir()
            print(vna_processed.as_posix() + add_dir_msg)
    if not dirfile.exists():
        print(dirfile.as_posix() + error_msg)
        if mkdir: 
            dirfile.mkdir()
            print(dirfile.as_posix() + add_dir_msg)
            
            
# check if data directories exist
check_if_dirs_exist(mkdir=True)
