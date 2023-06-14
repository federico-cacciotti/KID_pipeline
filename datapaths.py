from pathlib import Path

print('Looking for remote data folder...')
working_directory = Path('/Volumes/homes')
if not working_directory.exists():
    print('Remote data folder not found.')
    print('Switching to /home/cacciotf/data directory')
    working_directory = Path('/home/cacciotf/data')

target = working_directory / Path('target')
vna = working_directory / Path('vna')
target_S21 = working_directory / Path('target_S21')
vna_S21 = working_directory / Path('vna_S21')
dirfile = working_directory / Path('data_logger')
output_noise = working_directory / Path('noise')
cosmic_rays = working_directory / Path('cosmic_rays')
picolog = working_directory / Path('picolog')
anritsuMS2034B = working_directory / Path('anritsu_MS2034B')


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
    if not target_S21.exists():
        print(target_S21.as_posix() + error_msg)
        if mkdir: 
            target_S21.mkdir()
            print(target_S21.as_posix() + add_dir_msg)
    if not vna_S21.exists():
        print(vna_S21.as_posix() + error_msg)
        if mkdir: 
            vna_S21.mkdir()
            print(vna_S21.as_posix() + add_dir_msg)
    if not dirfile.exists():
        print(dirfile.as_posix() + error_msg)
        if mkdir: 
            dirfile.mkdir()
            print(dirfile.as_posix() + add_dir_msg)
    if not output_noise.exists():
        print(output_noise.as_posix() + error_msg)
        if mkdir: 
            output_noise.mkdir()
            print(output_noise.as_posix() + add_dir_msg)
    if not cosmic_rays.exists():
        print(cosmic_rays.as_posix() + error_msg)
        if mkdir: 
            cosmic_rays.mkdir()
            print(cosmic_rays.as_posix() + add_dir_msg)
    if not anritsuMS2034B.exists():
        print(anritsuMS2034B.as_posix() + error_msg)
        if mkdir: 
            anritsuMS2034B.mkdir()
            print(anritsuMS2034B.as_posix() + add_dir_msg)
    if not picolog.exists():
        print(picolog.as_posix() + error_msg)
        if mkdir: 
            picolog.mkdir()
            print(picolog.as_posix() + add_dir_msg)
            
            
            
# check if data directories exist
check_if_dirs_exist(mkdir=True)
