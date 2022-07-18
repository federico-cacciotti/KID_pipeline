from pathlib import Path


target = Path('./target')
vna = Path('./vna')
target_S21 = Path('./target_S21')
vna_S21 = Path('./vna_S21')
dirfile = Path('./data_logger')
output_noise = Path('./noise')
cosmic_rays = Path('./cosmic_rays')
picolog = Path('./picolog')
anritsuMS2034B = Path('./anritsu_MS2034B')

# check if directories already exists
error_msg = ' directory does not exists'
add_dir_msg = ' directory added'
def check_if_dirs_exist(mkdir = False):
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