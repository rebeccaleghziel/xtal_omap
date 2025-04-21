# Dataset class
# Rebecca Leghziel, Weizmann Institute of Science, February 2025

import numpy as np
import matplotlib.pyplot as plt
import py4DSTEM
import math
from ..stereofunctions import *
from .CrystalGroup import CrystalGroup

class Dataset:
    
    def __init__(self, crystal_class, sweep_size, coordinates, orientation_matrix, correlation, bragg_peaks):
        
        self.crystal_class = crystal_class
        self.sweep_size = sweep_size
        self.coordinates = coordinates # these are the original x,y coordinates
        self.orientation_matrix = orientation_matrix
        self.correlation = correlation
        self.bragg_peaks = bragg_peaks


    def run(self, corr_threshold = 0.6, dist_threshold = 25, theta_th = 3.5, min_groupsize = 50, filter = 'ON'):
        
        self.set_thresholds(corr_threshold, dist_threshold, theta_th, min_groupsize)
        self.run_clustering()
        
        if filter == 'ON':
            self.run_filtering()

        self.arrange_cluster_clockwise()

        if filter == 'ON':
             self.define_singlecrystals(option = 'filtered')
        elif filter == 'OFF':
            self.define_singlecrystals(option = 'other')
            
######################################################################################################################################################################            
# Methods in run() 

    
    def set_thresholds(self, corr_threshold = 0.6, dist_threshold = 25, theta_th = 3.5, min_groupsize = 50):
    
        self.corr_threshold = corr_threshold
        self.dist_threshold = dist_threshold
        self.theta_th = theta_th
        self.min_groupsize = min_groupsize

    def run_clustering(self):
        
        count = self.correlation_filter(self.corr_threshold)
        _, _, Z_corr = self.get_latticevectors(self.orientation_matrix[self.idx_corr])
        data_corr = self.array_forclustering(Z_corr, self.pixl_corr)
        self.cluster(data_corr, Z_corr)

    def run_filtering(self):
        
        groups_to_flip, idcs_to_flip = self.groups_to_rot()
        self.rotate_orientation(self.orientation_matrix[self.idx_corr], idcs_to_flip) #this should be a global variable
        _, _, Z_f = self.get_latticevectors(self.rotated_orientation)
        data_secondfilter = self.array_forclustering(Z_f, self.pixl_corr)
        self.cluster(data_secondfilter, Z_f)

    def arrange_cluster_clockwise(self):
        
        data_to_arrange = []
        
        for group in self.group_data:
            
            CoM_shift = self.shift_center(group[1])
            theta_CoM = self.theta_azimuthal_360(CoM_shift, [-1, 0])
            data_to_arrange.append([group[0], theta_CoM])
        
        sorted_data = sorted(data_to_arrange, key=lambda x: x[1])
        self.i_arranged = [x[0] for x in sorted_data]

    def define_singlecrystals(self, option = 'filtered'):
        
        self.crystal_unit = []
        i = 0

        for ind in self.i_arranged:
            
            indices_in = self.idx_corr[self.groups[ind]]
            indices_wrt_f = self.groups[ind]
            
            coordinates_in = self.pixl_corr[indices_wrt_f]
            correlation_in = self.correlation[indices_in]
            bragg_peaks_in = self.bragg_peaks[indices_in]

            if option == 'filtered':
                orientation_matrix_in = self.rotated_orientation[indices_wrt_f].squeeze()
            else:
                orientation_matrix_in = self.orientation_matrix[indices_in].squeeze() #
        
            single_crystal = CrystalGroup(
                    crystal_class = self.crystal_class,
                    sweep_size = self.sweep_size,
                    segment_number = i,
                    indices = indices_in, #with respect to the main dataset
                    coordinates = coordinates_in,
                    orientation_matrix = orientation_matrix_in,
                    correlation_number = correlation_in,
                    bragg_peaks = bragg_peaks_in
                )
        
            single_crystal.run()
        
            self.crystal_unit.append(single_crystal)
            
            i = i+1

######################################################################################################################################################################
    
    # plot_correlation creates a plot of all the pixels that are above correlation threshold
    # useful to understand which threshold to choose
    def plot_correlation(self, corr_threshold, fig_handle = 'None'):
        
        count = self.correlation_filter(corr_threshold)
        
        if fig_handle == 'None':
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        else:
            fig = fig_handle[0]
            ax = fig_handle[1]

        # plot based on correlation
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(0, self.sweep_size[1])
        ax.set_ylim(0, self.sweep_size[3])
        fig_width, fig_height = fig.get_size_inches() * fig.dpi  # Total figure size in pixels
        
        pixel_width = fig_width / self.sweep_size[1]
        pixel_height = fig_height / self.sweep_size[3]
        pixel_size = min(pixel_width, pixel_height)
        
        # pixl_corr = np.array(pixl_corr)
        ax.scatter(self.pixl_corr[:, 1], self.pixl_corr[:, 0], c='lightsteelblue', marker='s',  s=pixel_size , alpha=1)
        ax.set_title(f"Pixels with correlation > {corr_threshold}")
        print(f'Number of pixels with correlation > {corr_threshold} is {count}')

    def plot_map(self, colormap_type = 'Spectral', CoM = 'ON'):
        
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(20, 5))

        # Possible colormaps: viridis, inferno, Blues, RdYlBu, Spectral, tab20c
        colormap = plt.cm.get_cmap(colormap_type, len(self.crystal_unit))
        
        ax2 = initialize_stereogram(ax = ax2, style = 'default')
        ax3 = initialize_stereogram(ax = ax3, style = 'default')
        ax4 = initialize_stereogram(ax = ax4, style = 'default')
        
        ax2.set_title('c-axis projections')
        ax3.set_title('a-axis projections')
        ax4.set_title('b-axis projections')
        
        j = 0
        for i in range(0, len(self.crystal_unit)):
            color = colormap(j)
            j = j+1
        
            self.crystal_unit[i].plot_single_crystal(legend = 'false', ax = ax1, color_in = color)
            self.crystal_unit[i].add_projection(ax=ax2, lattice='Z', color = color)
            self.crystal_unit[i].add_projection(ax=ax3, lattice='X', color = color)
            self.crystal_unit[i].add_projection(ax=ax4, lattice='Y', color = color)

            if CoM == "ON":
                ax1.scatter(self.crystal_unit[i].CoM[0], self.crystal_unit[i].CoM[1], color='black', marker='+', s = 40)

######################################################################################################################################################################


    # correlation_filter returns the pixel coordinates, the indices and the number of pixels that are above a user defined correlation threshold
    def correlation_filter(self, corr_threshold):
        
        self.pixl_corr =[]
        self.idx_corr = []
        
        count = 0
        
        for i in range(0, len(self.coordinates)):
            if self.correlation[i] > corr_threshold:
                self.pixl_corr.append(self.coordinates[i]) #this is redundant here
                self.idx_corr.append(i)
                count+=1
                
        self.idx_corr = np.array(self.idx_corr)
        self.pixl_corr = np.array(self.pixl_corr)
        
        return count
    
    # get_latticevectors transforms the orientation matrix outputted from py4dstem into a, b, c of calcite
    # it's specific for calcite but can be modified for other crystal systems
    def get_latticevectors(self, orientation_matrix):
    
            a = 4.99
            b = 4.99
            c = 17.0615
        
            lattice_vectors = np.array([
                [a, 0, 0],
                [-a*np.sin(np.radians(30)), a*np.sqrt(3)/2, 0],
                [0, 0, c]
            ])
    
            result_matrix = np.dot(lattice_vectors, orientation_matrix)
    
            Xlattice_vec = result_matrix[0, :]
            Ylattice_vec = result_matrix[1, :]
            Zlattice_vec = result_matrix[2, :]
    
            return Xlattice_vec, Ylattice_vec, Zlattice_vec
    
    #angle between two vectors - taken from stereographic projections module
    
    def theta_azimuthal_360(self, vector, x):
        # x = [1, 0, 0]
        dot_product = np.dot(vector, x)
        magnitude_v = np.linalg.norm(vector)
        magnitude_x = np.linalg.norm(x)
        cos_theta = dot_product / (magnitude_v * magnitude_x)
        theta_radians = np.arccos(np.clip(cos_theta, -1.0, 1.0))
        cross_product = vector[0] * x[1] - vector[1] * x[0]
    
        # Adjust the angle to 0-360 degrees
        if cross_product < 0:  # Clockwise
            theta_radians = 2 * np.pi - theta_radians
    
        theta_d = np.degrees(theta_radians)  # Convert to degrees
        return theta_d
    
    def array_forclustering(self, Z_cf, pxl_cf):
    
        theta = []
        xytheta = []
        x = [1, 0]
        
        for Z in Z_cf:
    
            Z_in = [Z[0][0], Z[0][1]]
            theta_deg = self.theta_azimuthal_360(Z_in, x)
            theta.append(theta_deg)
    
        for i in range(0, len(pxl_cf)):
            joined = pxl_cf[i], theta[i]
            xytheta.append(joined)
    
        flattened_data = [(int(pair[0]), int(pair[1]), value) for pair, value in xytheta]
        flattened_data = np.array(flattened_data)
    
        return flattened_data
    
    def cluster(self, data_corr, Z_cf):
        
        self.groups = []
        self.group_data = []
        
        idcs_to_analyse = np.array(range(0, len(data_corr)))
        code = 0
        
        while len(idcs_to_analyse)>0:
            
            first_idx = idcs_to_analyse[0]
            other_idx = idcs_to_analyse[1:-1]
            
            cluster = []
            cluster.extend([idcs_to_analyse[0]])
            
            avg_az = data_corr[first_idx][2]
            avg_Z = Z_cf[first_idx]
            sum_V = avg_Z
            avg_CoM = [data_corr[first_idx][0], data_corr[first_idx][1]]
    
            N_divide = 1
            
            for i in other_idx:
                
                distance = math.sqrt((data_corr[i][0] - avg_CoM[0])**2 + (data_corr[i][1] - avg_CoM[1])**2) # element average
                theta_Z = angle_between_vectors(Z_cf[i][0], avg_Z[0]) # from stereofunctions
                theta_d = abs(data_corr[i][2] - avg_az)
                
                if theta_Z < self.theta_th and distance < self.dist_threshold:
                    N_divide = N_divide + 1
                    
                    cluster.extend([i])
                    
                    # Update avg_Z
                    sum_V = sum_V + Z_cf[i]
                    avg_Z = sum_V/N_divide
                    avg_az = np.mean(data_corr[cluster][:, 2])
                    avg_CoM = [np.mean(data_corr[cluster][:, 0]), np.mean(data_corr[cluster][:, 1])]
        
            if len(cluster)>self.min_groupsize:
                self.groups.append(cluster)
                group_data_append = [code, avg_CoM, avg_az, avg_Z]
                self.group_data.append(group_data_append)
                code = code+1
                
            indices_to_delete = np.where(np.isin(idcs_to_analyse, cluster))[0]
            idcs_to_analyse = np.delete(idcs_to_analyse, indices_to_delete)

    def groups_to_rot(self):

        group_code_flip =[]
        idcs_to_flip = []
        
        for group in self.group_data:
            
            CoM_shift = self.shift_center(group[1])
            avg_Z = group[3]
            
            theta = angle_between_vectors(CoM_shift, avg_Z[0, :2])
            
            if theta > 90:
                group_code = group[0]
                group_code_flip.append(group_code)
                idcs_to_flip.extend([self.groups[group_code]])
    
        idcs_to_flip_flattened = [item for sublist in idcs_to_flip for item in sublist]
        
        return group_code_flip, idcs_to_flip_flattened

    def reflect_180yx(self, orientation_matrix):
        theta = np.radians(180)
        Rxy_180 = np.array([
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, 1]
        ])
    
        OM_rot_tot = np.dot(orientation_matrix, Rxy_180)
    
        return OM_rot_tot

    def shift_center(self, CoM):

        lx = (self.sweep_size[1] - self.sweep_size[0])/2
        ly = (self.sweep_size[3] - self.sweep_size[2])/2
        
        CoM_shift = [CoM[1] - lx, CoM[0] - ly]
    
        return CoM_shift #x,y
    
    def rotate_orientation(self, orientation_f, idcs_to_flip_flattened):
    
        self.rotated_orientation = []
    
        for i, orientation in enumerate(orientation_f):
            
            if i in idcs_to_flip_flattened:
                om = self.reflect_180yx(orientation)
            else:
                om = orientation
    
            self.rotated_orientation.append(om)
        
        self.rotated_orientation = np.array(self.rotated_orientation)
        