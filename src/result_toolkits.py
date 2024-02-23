# import package
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import warnings
import nums_from_string as nfs
import re

# Define class
class phantom_analysis:
    def __init__(self,analysis_type , time, data_dict, theta, radius, column_names) -> None:
        self.analysis_type = analysis_type
        self.time = time
        self.radius = radius
        self.data_dict = data_dict
        self.theta = theta
        self.column_names = column_names
        self.type_check()
        
    
    def type_check(self) -> None:
        if self.analysis_type == 'Pitch_analysis':
            assert type(self.radius) == np.ndarray, f"SelfCheckError: Radius in the type {self.analysis_type} should be 'ndarray':, rather then {type(self.radius)}"
        elif self.analysis_type == 'Spiral_analysis':
            assert type(self.radius) == float, f"SelfCheckError: Radius in the type {self.analysis_type} should be 'float', rather then {type(self.radius)}"
        else:
            warnings.warn("SelfCheckWarn: Unknown type of analysis. Please check before further operation.")
    
    @staticmethod
    def read_H5DF(filepath):
        def read_dict(f, name, K=str):
            dictionary = {}
            g = f[name]
            for key,val in g.items():
                if type(val[()]) == np.bytes_:
                    dictionary[K(key)] = val[()].decode('utf-8')
                elif type(val[()]) == np.ndarray:
                    dictionary[K(key)] = val[()].copy().T # Trasepose because of the stroge system
                else:
                    dictionary[K(key)] = val[()]
            return dictionary 
            
        with h5py.File(filepath, 'r') as f:
            analysis_type = f['struct_type'][()].decode('utf-8')
            if analysis_type == 'Pitch_analysis':
                radius = f['radius'][:].copy()
            elif analysis_type == 'Spiral_analysis':
                radius = f['radius'][()]
            else:
                warnings.warn("ReadFileWarn: Unknown type of analysis. Please check before further operation.")
                radius = f['radius']
            time = f['time'][()]
            theta = f['theta'][:]
            
            data_dict = read_dict(f,'data_dict',int)
            column_names = read_dict(f, "column_names", int)
            
            data = phantom_analysis(analysis_type , time, data_dict, theta, radius, column_names)
            return data
        
    @staticmethod
    def read_dat(filepath):
        def drop_column_names(line):
            column_pattern = r'(?<=\[)([^\]]{1,20})(?=\])'
            column_names =  sorted(list(set(re.findall(column_pattern,line))))
            if not column_names:
                return False
            else:
                column_names = ['['+ name + ']' for name in column_names]
                return column_names
        
        def drop_time(line):
            if r'# Analysis data at t =' in line:
                multi_part, expo_part = nfs.get_nums(line)
                decimal = multi_part * (10**expo_part)
                return decimal
            else:
                return False
            
        def drop_radius(line):
            if r'# theta_rad =' in line:
                multi_part, expo_part = nfs.get_nums(line)
                decimal = multi_part * (10**expo_part)
                return decimal
            else:
                return False
        
        def generate_column_name_dict(column_name_list):
            result_dict = {}
            for i,name in enumerate(column_name_list):
                result_dict[i+1] = name
            return result_dict
                
        def Spiral_mode(filepath):
            data_array = np.genfromtxt(filepath,invalid_raise=False)
            time = None
            column_names_list = None
            radius = None
            with open(filepath,'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#'):
                        time_temp = drop_time(line)
                        if not time_temp:
                            pass
                        else:
                            time = time_temp
                            continue
                        
                        radius_temp = drop_radius(line)
                        if not radius_temp:
                            pass
                        else:
                            radius = radius_temp
                            continue
                        
                        column_names_list_temp = drop_column_names(line)
                        if not column_names_list_temp:
                            pass
                        else:
                            column_names_list = column_names_list_temp.copy()
                            continue
            column_names = generate_column_name_dict(column_names_list)
            
            theta = data_array[:,0]
            data_dict = {}
            for i in range(1,len(column_names_list)):
                data_dict[column_names[i+1]] = data_array[:,i].copy()
                
            return phantom_analysis('Spiral_analysis',time,data_dict,theta,radius,column_names)
        
        def Pitch_mode(filepath):
            data_array = np.genfromtxt(filepath,invalid_raise=False)
            time = None
            column_names_list = None
            radius_list = []
            with open(filepath,'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#'):
                        if not time:
                            time_temp = drop_time(line)
                            if not time_temp:
                                pass
                            else:
                                time = time_temp
                                continue
                        if not column_names_list:
                            column_names_list_temp = drop_column_names(line)
                            if not column_names_list_temp:
                                pass
                            else:
                                column_names_list = column_names_list_temp.copy()
                                continue
                        
                        radius_temp = drop_radius(line)
                        if not radius_temp:
                            pass
                        else:
                            radius_list.append(radius_temp)
                            continue
            column_names = generate_column_name_dict(column_names_list)
            radius = np.array(radius_list)
            block_num = len(radius)
            block_length = int(round(data_array.shape[0]/block_num))
            theta = data_array[:block_length,0]
            data_dict = {}
            for i in range(1,len(column_names_list)):
                data_dict[column_names[i+1]] = data_array[:,i].reshape(block_num,block_length)
                
            return phantom_analysis('Pitch_analysis',time,data_dict,theta,radius,column_names)
  
        if 'pitch_' in filepath:
            return Pitch_mode(filepath)
        elif 'spiral_' in filepath:
            return Spiral_mode(filepath)
        else:
            raise AttributeError
            # warnings.warn("ReadFileWarn: Unknown type of analysis. Please check before further operation.")
            # analysis_type = 'Unknown_analysis'
            
    def transfer_cgs(self,umass=1.9891E+33, udist=1.496E+13, utime=5.0227287E+06, year=True,au=True):
        try:
            _ = self._cgs
        except AttributeError:
            regex = r"(cm\$^{-)(\d+)(\}\])"
            def replace_grident_exp(match):
                exponent = int(match.group(2)) - 1
                new_exponent_str = f"{match.group(1)}{exponent}{match.group(3)}"
                return new_exponent_str
            usigma = umass/(udist**2)
            urho = umass/(udist**3)
            uv = udist/utime
            
            if year:
                self.time *= (utime/31536000)
            else:
                self.time *= utime
            if not au:
                self.radius *= udist
                
            self.column_unit = {}
            for key,value in self.column_names.items():
                if ('Sigma' in value) or (('sigma' in value)):
                    self.data_dict[key] *= usigma
                    if '_g' in value:
                        self.column_unit[key] = r'$\Sigma_g$\,[g cm$^{-2}]$'
                    elif '_d' in value:
                        self.column_unit[key] = r'$\Sigma_d$\,[g cm$^{-2}]$'
                    else:
                        self.column_unit[key] = r'$\Sigma$\,[g cm$^{-2}]$'
                elif 'rho' in value:
                    self.data_dict[key] *= urho
                    if '_g' in value:
                        self.column_unit[key] = r'$\rho_g$\,[g cm$^{-3}]$'
                    elif '_d' in value:
                        self.column_unit[key] = r'$\rho_d$\,[g cm$^{-3}]$'
                    else:
                        self.column_unit[key] = r'$\rho$\,[g cm$^{-3}]$'
                elif 'vr' in value:
                    self.data_dict[key] *= uv
                    if '_g' in value:
                        self.column_unit[key] = r'$v_{r,g}$\,[cm s$^{-1}]$'
                    elif '_d' in value:
                        self.column_unit[key] = r'$v_{r,d}$\,[cm s$^{-1}]$'
                    else:
                        self.column_unit[key] = r'$v_{r}$\,[cm s$^{-1}]$'
                elif 'vphi' in value:
                    self.data_dict[key] *= uv
                    if '_g' in value:
                        self.column_unit[key] = r'$v_{\phi,g}$\,[cm s$^{-1}]$'
                    elif '_d' in value:
                        self.column_unit[key] = r'$v_{\phi,d}$\,[cm s$^{-1}]$'
                    else:
                        self.column_unit[key] = r'$v_{\phi}$\,[cm s$^{-1}]$'
                else:
                    self.column_unit[key] = ''
                
                if '∇' in value:
                    self.data_dict[key] /= udist
                    self.column_unit[key] = r'$\nabla$ ' + re.sub(regex,replace_grident_exp,self.column_unit[key])
        else:
            pass
        self._cgs = True

    def add_more_label(self,column_index,label):
        try:
            _ = self.column_unit
        except AttributeError:
            self.column_unit = {}
        else:
            try:
                _ = self.column_unit[column_index]
            except KeyError:
                pass
            else:
                if not self.column_unit[column_index] == '': 
                    while True:
                        print(f'Old column unit {self.column_unit[column_index]} has found. Are you sure you want to replace it? [y/n]')
                        input_val = str(input())
                        if input_val == 'y':
                            break
                        elif input_val == 'n':
                            return
                        else:
                            pass
        self.column_unit[column_index] = label
            