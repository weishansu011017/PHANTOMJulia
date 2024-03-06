# import package
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import h5py
import warnings
import nums_from_string as nfs
import re
import cv2
from scipy import ndimage, optimize

#Include LaTeX rendering
plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 13,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsfonts}"
})

# Useful function
def value2closestvalueindex(array, value):
        return np.int64(np.nanargmin(np.abs(array - value)))

def openinteractive():
        if not plt.isinteractive():
            plt.ion()
            
def Get_vmaxmin(array,minzero = False):
            median = np.nanmedian(array)
            std = np.nanstd(array)
            vmax = median + 3*std
            vmin = median - 3*std
            if minzero and vmin<0.0:
                vmin = 0.0
            return vmax,vmin

def Theoretical_sigma(Sigma_type="gas"):
    if Sigma_type=="gas":
        Sigma0 = 2.21E-01
    elif Sigma_type=="dust":
        Sigma0 = 2.21E-03
    r = np.linspace(10.0,100.0,301)
    Sigma = Sigma0*(r/100.0)**(-1.5)*(1-np.sqrt(10.0/r))
    return Sigma,r
    

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
                    dictionary[K(key)] = val[()].copy().T # Trasepose because of the storage system
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
                elif 'vtheta' in value:
                    self.data_dict[key] *= uv
                    if '_g' in value:
                        self.column_unit[key] = r'$v_{\theta,g}$\,[cm s$^{-1}]$'
                    elif '_d' in value:
                        self.column_unit[key] = r'$v_{\theta,d}$\,[cm s$^{-1}]$'
                    else:
                        self.column_unit[key] = r'$v_{\theta}$\,[cm s$^{-1}]$'
                elif 'e_' in value:
                    if '_g' in value:
                        self.column_unit[key] = r'$e_g$'
                    elif '_d' in value:
                        self.column_unit[key] = r'$e_d$'
                    else:
                        self.column_unit[key] = r'$e$'
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
        
    def setup_fig(self,nrows,ncols,index):
        openinteractive()
        if not hasattr(self,"fig") or not hasattr(self,"ax"):
            self.fig= plt.figure(figsize=(10,6))
            self.ax = self.fig.add_subplot(nrows,ncols,index)

    def extract_radius_average(self,column_index,draw=False,Log_scale=True,Theoretical_value=True):
        self.transfer_cgs()
        try:
            full_image = self.data_dict[column_index]
        except KeyError:
            raise KeyError("Invalid column index.")
        radius_array = self.radius
        radius_length = len(radius_array)
        result_array = np.zeros(radius_length,dtype=float)
        for i in range(radius_length):
            result_array[i] = np.average(full_image[i,:])
        if draw:
            openinteractive()
            self.setup_fig(1,1,1)
            self.ax.cla()
            self.ax.plot(radius_array,result_array,c="orange",label=self.column_names[column_index])
            if Theoretical_value:
                if column_index==2:
                    T_sigma, T_r = Theoretical_sigma("gas")
                    self.ax.plot(T_r,T_sigma,c="red",label="Theoretical value")
                elif column_index==3:
                    T_sigma, T_r = Theoretical_sigma("dust")
                    self.ax.plot(T_r,T_sigma,c="red",label="Theoretical value")
                else:
                    pass
            self.ax.set_xlabel(r"$r$ (au)")
            self.ax.set_xlim((radius_array[0],radius_array[-1]))
            self.ax.set_ylabel(self.column_unit[column_index])
            self.ax.legend()
            if Log_scale:
                self.ax.set_yscale('log')
            plt.draw()
        return result_array
    
    def close_fig(self):
        try:
            verify = self.fig
            verify = self.ax
        except AttributeError:
            pass
        else:
            del self.fig,self.ax

    def setup_polar(self):
        openinteractive()
        if not hasattr(self,"polar_fig") or not hasattr(self,"polar_ax"):
            self.polar_fig= plt.figure(figsize=(10,6))
            self.polar_ax = self.polar_fig.add_subplot(projection='polar')
        
    
    def polar_plot(self,column_index,label=r'',Log_mode=False):
        def Extend_array_grid(array):
            s = array[0]
            e = array[-1]
            new_array = np.linspace(s,e,len(array)+1)
            return new_array
        if not self.analysis_type == 'Pitch_analysis':
            print(f'PlotError: The analysis type {self.analysis_type} is not allowed for plotting the polar plot.')
        self.transfer_cgs()
        radius_array = self.radius
        theta_array = self.theta
        time = self.time
        z = self.data_dict[column_index]
        vmax,vmin = Get_vmaxmin(z,np.nanmin(z)>-1e-7)
        if label == r'':
            colorbar_unit = self.column_unit[column_index].replace('∇',r'$\nabla$')
        else:
            colorbar_unit = label.replace('∇',r'$\nabla$')
        if Log_mode:
            color_norm = mcolors.LogNorm(vmin,vmax)
        else:
            color_norm = mcolors.Normalize(vmin,vmax)
        extendx,extendy = np.meshgrid(Extend_array_grid(theta_array),Extend_array_grid(radius_array))
        props = dict(boxstyle='round', facecolor='black', alpha=0.5)
        self.setup_polar()
        self.polar_ax.cla()
        
        cont = self.polar_ax.pcolormesh(extendx,extendy,z,cmap = 'inferno',norm = color_norm,rasterized=True)
        
        colarbar_ax = None
        for ax in self.polar_fig.axes:
            if ax.get_label() == '<colorbar>':
                colarbar_ax = ax
        if colarbar_ax is None:
            colorbar = plt.colorbar(cont)
        else:
            colorbar = plt.colorbar(cont,cax = colarbar_ax)
        colorbar.set_label(colorbar_unit)
        
        self.polar_ax.set_rorigin(-10)
        self.polar_ax.text(0,1.05,'Time = %d yrs' % (time),transform=self.polar_ax.transAxes,c='white', fontsize=14, verticalalignment='top', bbox=props)
        self.polar_ax.text(0.7,1.05,self.column_names[column_index].replace('∇',r'$\nabla$'),transform=self.polar_ax.transAxes,c='white', fontsize=14, verticalalignment='top', bbox=props)
        plt.draw()
        
    def close_polar(self):
        try:
            verify = self.polar_fig
            verify = self.polar_ax
        except AttributeError:
            pass
        else:
            del self.polar_fig,self.polar_ax
                
            
    
    def spiral_detection(self,image_index = 2,gaussian_sigma = 2.2,radius_range=[65.0,100.0], fitting_model="Logarithmic", initial_p=[100.0,-0.5],draw=False):
        def Logarithmic_spiral(r,a,k_0):
            k = k_0
            spiral = (1/k)*np.log(r/a)%2*np.pi
            return spiral
        def quasi_Logarithmic_spiral(r,a,k_0,kappa):
            k = k_0 + kappa*(r)
            spiral = (1/k)*np.log(r/a)%2*np.pi
            return spiral
        
        def dilate_erose(two_values_image):
            kernel = np.ones((5,5),np.uint8)
            dilated_image = cv2.dilate(two_values_image, kernel, iterations=2)
            eroded_image = cv2.erode(dilated_image, kernel, iterations=2)
            binary_image = eroded_image/255
            y_coords, x_coords = np.where(binary_image == 1)
            points = np.column_stack((x_coords, y_coords))
            return points
        
        def sort_points(detected_points):
            sorted_indices = np.argsort(detected_points[1,:])
            result = detected_points[:,sorted_indices]
            return result
        
        def two_spiral_classification(detected_points):
            sorted_detected_points = sort_points(detected_points)
            sorted_detected_points_t = sorted_detected_points.transpose()
            bukket1 = []
            bukket2 = []
            for i in range(sorted_detected_points_t.shape[0]):
                point = sorted_detected_points_t[i,:]
                if not bukket1:
                    bukket1.append(point)
                else:
                    if not bukket2:
                        previous_point = bukket1[i-1]
                        if np.abs(point[0]-previous_point[0])>np.pi/2:
                            bukket2.append(point)
                        else:
                            bukket1.append(point)
                    else:
                        previous_point1 = bukket1[-1]
                        previous_point2 = bukket2[-1]
                        sub_azi1 = np.abs(point[0]-previous_point1[0])
                        sub_azi2 = np.abs(point[0]-previous_point2[0])
                        if sub_azi1 < sub_azi2:
                            bukket1.append(point)
                        else:
                            bukket2.append(point)
            bukket1_array = np.array(bukket1).transpose()
            bukket2_array = np.array(bukket2).transpose()
            return bukket1_array,bukket2_array
        
        def spiral_fitting(points,spiral_model,initial_p,bounds):
            '''
            Notice that points should act as:
             
            points
            >>> array(
                [[theta1,theta2...],
                [r1,r2....]]
            )
            '''
            
            print(f"Input parameters: {initial_params}")
            theta = points[0,:]
            r = points[1,:]
            params, covariance = optimize.curve_fit(spiral_model, r, theta, p0=initial_p,bounds=bounds,method='trf',maxfev=1e7)
            print(f"Best-fit parameters: {params}")
            return params,covariance
            
        def pitch_angle(k0,kappa=0.0,r=50.0):
            k = k0 + kappa*r
            beta = np.arctan(np.abs(k))*(180/np.pi)
            return beta
        
        if self.analysis_type != 'Pitch_analysis':
            raise ValueError(f"Spiral Dection only valid when analysis_type is 'Pitch_analysis', but we get {self.analysis_type}!")
        if not type(radius_range) in [np.ndarray,list,tuple]:
            raise TypeError(f"Type of radius_range should be either numpy.ndarray, list, tuple, but we get {type(radius_range)}")
        
        # Getting image
        full_image = self.data_dict[image_index]
        self.transfer_cgs()
        
        # Preparation of spiral detection
        radius_array = self.radius
        theta_array = self.theta
        radius_index_range = np.zeros(2)
        for i in range(2):
            radius_index_range[i] = value2closestvalueindex(radius_array,radius_range[i])
        radius_index_range = np.int64(radius_index_range)
        try:
            image = full_image[radius_index_range[0]:radius_index_range[1],:]
            cutoff_radius_array = radius_array[radius_index_range[0]:radius_index_range[1]]
        except KeyError:
            raise ValueError(f"Invalid radius range! The value of radius only valid in {radius_array[0]} to {radius_array[-1]}, but {radius_range} was gotten!")
        
        # Spiral detection
        image = ndimage.gaussian_filter(image,sigma=gaussian_sigma)
        image_level = (255*(image / image.max())).astype(np.uint8)
        edges = cv2.Canny(image=image_level, threshold1=100, threshold2=200)
        detected_points_index = dilate_erose(edges).transpose()
        detected_points_array = np.zeros_like(detected_points_index,dtype=float)
        detected_points_array[0,:] = theta_array[detected_points_index[0,:]]
        detected_points_array[1,:] = cutoff_radius_array[detected_points_index[1,:]]
        
        # Point classification
        spiral1,spiral2 = two_spiral_classification(detected_points_array)
    
        # Spiral fitting
        if fitting_model=="Logarithmic":
            spiral_model = Logarithmic_spiral
            lower_bound = [50.0,-np.pi/4]
            upper_bound = [250.0,0.0]
            initial_params = initial_p.copy()
        elif fitting_model=="quasi-Logarithmic":
            spiral_model = quasi_Logarithmic_spiral
            lower_bound = [50.0,-np.pi/4,-0.05]
            upper_bound = [250.0,0.0,0.05]
            initial_params = initial_p.copy()
            initial_params.append(0.0)
        else:
            raise ValueError("Invalid model of spiral.")

        k_list = []
        kappa_list = []
        
        if spiral1.size > 0:
            self.spiral1_params, self.spiral1_covariance= spiral_fitting(spiral1,spiral_model,initial_params,(lower_bound,upper_bound))
            spiral1_model = spiral_model(radius_array,*self.spiral1_params)
            k_list.append(self.spiral1_params[1])
            if fitting_model=="quasi-Logarithmic":
                kappa_list.append(self.spiral1_params[2])
        if spiral2.size > 0:
            self.spiral2_params,self.spiral2_covariance = spiral_fitting(spiral2,spiral_model,initial_params,(lower_bound,upper_bound))
            spiral2_model = spiral_model(radius_array,*self.spiral2_params)
            k_list.append(self.spiral2_params[1])
            if fitting_model=="quasi-Logarithmic":
                kappa_list.append(self.spiral2_params[2])
        
        k_average = np.mean(np.array(k_list))
        if np.array(kappa_list).size >0:
            kappa_average = np.mean(np.array(kappa_list))
        else:
            kappa_average = 0.0
            
        self.pitch_angle = pitch_angle(k_average,kappa_average,75)
        print(f"The pitch angle of spiral is {self.pitch_angle}.")
        
        
        # Plotting
        if draw:
            self.polar_plot(image_index)
            for artist in self.polar_ax.get_children():
                if isinstance(artist, matplotlib.collections.PathCollection):
                    artist.remove()
            for line in self.polar_ax.lines:
                line.remove()
            if spiral1.size > 0:
                self.polar_ax.scatter(spiral1[0,:],spiral1[1,:],label="Traced_spiral",c='cyan')
                self.polar_ax.plot(spiral1_model,radius_array,label="Traced_spiral",c='b')
            if spiral2.size > 0:
                self.polar_ax.scatter(spiral2[0,:],spiral2[1,:],label="Traced_spiral",c='lime')
                self.polar_ax.plot(spiral2_model,radius_array,label="Traced_spiral",c='g')

            plt.draw()

        
        

    
    
        
        
        
            