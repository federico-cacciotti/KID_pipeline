from . G31_KID_pipeline import *

__version__ = '1.1.0'
__author__ = 'Federico Cacciotti'




'''
            Data_paths class
'''
from pathlib import Path
class Data_paths():
    
    def __init__(self):
        self.target = Path('target')
        self.vna = Path('vna')
        self.target_S21 = Path('target_S21')
        self.vna_S21 = Path('vna_S21')
        self.dirfile = Path('data_logger')
        self.output_noise = Path('noise')
        self.cosmic_rays = Path('cosmic_rays')
        self.picolog = Path('picolog')
        self.anritsuMS2034B = Path('anritsu_MS2034B')
    
    def __call__(self, working_directory=''):
        self.target = Path(working_directory) / Path('target')
        self.vna = Path(working_directory) / Path('vna')
        self.target_S21 = Path(working_directory) / Path('target_S21')
        self.vna_S21 = Path(working_directory) / Path('vna_S21')
        self.dirfile = Path(working_directory) / Path('data_logger')
        self.output_noise = Path(working_directory) / Path('noise')
        self.cosmic_rays = Path(working_directory) / Path('cosmic_rays')
        self.picolog = Path(working_directory) / Path('picolog')
        self.anritsuMS2034B = Path(working_directory) / Path('anritsu_MS2034B')
        
        # check if data directories exist
        self.check_if_dirs_exist(mkdir=True)


    def check_if_dirs_exist(self, mkdir = False):
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
        
        error_msg = ' directory does not exists'
        add_dir_msg = ' directory added'
        
        if not self.target.exists():
            print(self.target.as_posix() + error_msg)
            if mkdir: 
                self.target.mkdir()
                print(self.target.as_posix() + add_dir_msg)
        if not self.vna.exists():
            print(self.vna.as_posix() + error_msg)
            if mkdir: 
                self.vna.mkdir()
                print(self.vna.as_posix() + add_dir_msg)
        if not self.target_S21.exists():
            print(self.target_S21.as_posix() + error_msg)
            if mkdir: 
                self.target_S21.mkdir()
                print(self.target_S21.as_posix() + add_dir_msg)
        if not self.vna_S21.exists():
            print(self.vna_S21.as_posix() + error_msg)
            if mkdir: 
                self.vna_S21.mkdir()
                print(self.vna_S21.as_posix() + add_dir_msg)
        if not self.dirfile.exists():
            print(self.dirfile.as_posix() + error_msg)
            if mkdir: 
                self.dirfile.mkdir()
                print(self.dirfile.as_posix() + add_dir_msg)
        if not self.output_noise.exists():
            print(self.output_noise.as_posix() + error_msg)
            if mkdir: 
                self.output_noise.mkdir()
                print(self.output_noise.as_posix() + add_dir_msg)
        if not self.cosmic_rays.exists():
            print(self.cosmic_rays.as_posix() + error_msg)
            if mkdir: 
                self.cosmic_rays.mkdir()
                print(self.cosmic_rays.as_posix() + add_dir_msg)
        if not self.anritsuMS2034B.exists():
            print(self.anritsuMS2034B.as_posix() + error_msg)
            if mkdir: 
                self.anritsuMS2034B.mkdir()
                print(self.anritsuMS2034B.as_posix() + add_dir_msg)
        if not self.picolog.exists():
            print(self.picolog.as_posix() + error_msg)
            if mkdir: 
                self.picolog.mkdir()
                print(self.picolog.as_posix() + add_dir_msg)
                

paths = Data_paths()