# CrystalGroup class
# Rebecca Leghziel

# This class deals with a single crystal group and has methods for:
# 1. Plotting the crystal
# 2. Obtaining a histogram for the correlation of the matches for each pixel in the group
# 3. Averaging the correlation and orientation matrix, as well as obtaining standard deviations for the group
# 4. Computing angular parameters which are important for the coccolith analysis
# 6. Computing the stereogram for the single crystal


# TO DO: Add text to plot_stereo function, add V/R classification

import numpy as np
import matplotlib.pyplot as plt
from ..stereofunctions import *
import py4DSTEM


class CrystalGroup:
    
    def __init__(self, crystal_class, sweep_size, segment_number, indices, coordinates, orientation_matrix, correlation_number, bragg_peaks):
        self.crystal_class = crystal_class
        self.sweep_size = sweep_size
        self.segment_number = segment_number
        self.indices = indices
        self.coordinates = coordinates # these are the original x,y coordinates
        self.orientation_matrix = orientation_matrix
        self.correlation_number = correlation_number
        self.bragg_peaks = bragg_peaks

    def run(self):
        self.run_orientation()
        # self.calculate_angles()
        self.corr_avg = np.mean(self.correlation_number)
        self.compute_CoM()

    def run_orientation(self):
        self.get_lattice_vectors()
        self.lattice_vect_avg()
        self.compute_projections()

    # Methods used when running the crystal group to define the class attributes
    
    def get_lattice_vectors(self):

        a = 4.99
        b = 4.99
        c = 17.0615
    
        lattice_vectors = np.array([
            [a, 0, 0],
            [-a*np.sin(np.radians(30)), a*np.sqrt(3)/2, 0],
            [0, 0, c]
        ])

        result_matrix = np.dot(lattice_vectors, self.orientation_matrix)

        self.Xlattice_vec = result_matrix[0, :]
        self.Ylattice_vec = result_matrix[1, :]
        self.Zlattice_vec = result_matrix[2, :]

        return self.Xlattice_vec, self.Ylattice_vec, self.Zlattice_vec

    def lattice_vect_avg(self):
        
        self.Xlattice_avg = np.mean(self.Xlattice_vec, axis=0)
        self.Ylattice_avg = np.mean(self.Ylattice_vec, axis=0)
        self.Zlattice_avg = np.mean(self.Zlattice_vec, axis=0)

        return self.Xlattice_avg, self.Ylattice_avg, self.Zlattice_avg

    def flip_lattice(self):

        self.Xlattice_avg = self.rot_180x_180y(self.Xlattice_avg)
        self.Ylattice_avg = self.rot_180x_180y(self.Ylattice_avg)
        self.Zlattice_avg = self.rot_180x_180y(self.Zlattice_avg)

        return self.Xlattice_avg, self.Ylattice_avg, self.Zlattice_avg
    
    def compute_projections(self):

        self.proj_a = stereoproj(self.Xlattice_avg)
        self.proj_b = stereoproj(self.Ylattice_avg)
        self.proj_c = stereoproj(self.Zlattice_avg)

    def compute_CoM(self):
        y_coords = self.coordinates[:, 0]
        x_coords = self.coordinates[:, 1]

        mean_x = np.mean(x_coords)
        mean_y = np.mean(y_coords)

        self.CoM = (mean_x, mean_y)

    def apply_rotation_and_shift(self, angle, center):

        self.shift_crystal(center)
        self.rotate_pixels(angle)
        self.Xlattice_avg = self.rotate_xy(self.Xlattice_avg, -angle)
        self.Ylattice_avg = self.rotate_xy(self.Ylattice_avg, -angle)
        self.Zlattice_avg = self.rotate_xy(self.Zlattice_avg, -angle)
        self.compute_CoM()

    # Functions used in main methods

    # FIX: POLAR and AZIMUTHAL
    
    def polar(self, lattice_vec):
        
        z_axis = np.array([0, 0, 1])
        vectors = lattice_vec
        vectors = vectors.reshape(-1, 3)
        
        # Calculate the magnitude of each vector
        vector_magnitudes = np.linalg.norm(vectors, axis=1) # these are the lattice magnitude of z = 17.0615
        z_axis_magnitude = np.linalg.norm(z_axis) # Magnitude of the z-axis unit vector (which is 1)
        dot_product = np.dot(vectors, z_axis)
        cos_theta = dot_product / (vector_magnitudes * z_axis_magnitude)
        theta_rad = np.arccos(cos_theta)
        theta_V = np.degrees(theta_rad)

        average_theta = np.mean(theta_V)
        std_theta = np.std(theta_V)
        
        return average_theta, std_theta

    def azimuthal(self, lattice_vec):
        
        vectors = lattice_vec.reshape(-1, 3)
        projected_vector = np.column_stack((vectors[:, 0], vectors[:, 1]))
    
        # Compute angles using arctan2 to get values in [-180, 180]
        theta_rad = np.arctan2(projected_vector[:, 1], projected_vector[:, 0])
        theta_deg = np.degrees(theta_rad)
    
        # Convert negative angles to [0, 360] range
        theta_deg = np.mod(theta_deg, 360)
    
        average_theta = np.mean(theta_deg)
        std_theta = np.std(theta_deg)

        return average_theta, std_theta

    def rot_180x_180y(self, input_V):
        
        theta = np.radians(180)
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)]
        ])
            
        Ry =  np.array([
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)]
        ])
        
        Rxy_180 = np.dot(Rx, Ry)
        Vrot = np.dot(Rxy_180, input_V)
        
        return Vrot

    def rotate_xy(self, vector, phi):
            
        # Rotation matrix in xy-plane
        R = np.array([
            [np.cos(phi), -np.sin(phi), 0],
            [np.sin(phi), np.cos(phi), 0],
            [0, 0, 1]
        ])
        
        # Perform matrix-vector multiplication
        rotated_v = R @ vector
        
        return rotated_v

    def shift_crystal(self, center):
        
        coordinates_new = []
        xc = center[0]
        yc = center[1]
        
        for coord in self.coordinates:
            cx = coord[1] - xc
            cy = coord[0] - yc
            c = [cy, cx]
    
            coordinates_new.append(c) 

        coordinates_new = np.array(coordinates_new)
        self.coordinates = coordinates_new
       
        deltaX1 = abs(self.sweep_size[0] - xc)
        deltaX2 = abs(self.sweep_size[1] - xc)
        deltaY1 = abs(self.sweep_size[2] - yc)
        deltaY2 = abs(self.sweep_size[3] - yc)

        
        self.sweep_size = [-deltaX1, deltaX2, -deltaY1, deltaY2]
        # self.compute_CoM()

    def rotate_pixels(self, angle): #along 0, 0
        
        theta = angle
        R = np.array([[np.cos(theta), -np.sin(theta)],
            [np.sin(theta),  np.cos(theta)]])

        rotated_coordinates = []

        for coord in self.coordinates:
            
            cy = coord[0]
            cx = coord[1]
            coord_to_rot = [cy, cx]
            
            rot_coord = R @ coord_to_rot
            
            rotated_coordinates.append([rot_coord[0], rot_coord[1]])
        
        rotated_coordinates = np.array(rotated_coordinates)
        self.coordinates = rotated_coordinates
    
    # METHODS FOR DATA VISUALIZATION
    
    def plot_single_crystal(self, ax = None, legend = 'true', color_in = 'blue'):
        """Plot the x, y coordinates of the crystal group."""
        if ax is None:
            ax = plt.gca()

        pixel_size = 1
        point_size = pixel_size * (plt.rcParams['figure.dpi'] / 72.0)**2
        
        ax.scatter(self.coordinates[:, 1], self.coordinates[:, 0], label=f'Cluster {self.segment_number}', s=point_size, marker='s', alpha=1, color = color_in)
        
        # Maintain the original scale
        ax.set_xlim(self.sweep_size[0], self.sweep_size[1])
        ax.set_ylim(self.sweep_size[2], self.sweep_size[3])
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        # ax.set_title(f'Crystal Group {self.segment_number} Coordinates')
        
        if legend == 'true':
            ax.legend()
    
    def histogram_of_correlation(self, ax = None, bins=10):
        """Plot a histogram of the correlation numbers."""
        if ax is None:
            ax = plt.gca()
        
        ax.hist(self.correlation_number, bins=bins)
        ax.set_xlabel('Correlation Number')
        ax.set_ylabel('Frequency')
        ax.set_title(f'Correlation - Cluster {self.segment_number}')


    def plot_stereo(self, axis, xcolor = 'steelblue', ycolor = 'limegreen', zcolor = 'red'):
        
        axis = initialize_stereogram(ax = axis, style = 'default')
        self.add_projection(ax=axis, lattice='X', color = xcolor)
        self.add_projection(ax=axis, lattice='Y', color = ycolor)
        self.add_projection(ax=axis, lattice='Z', color = zcolor)

        
    def add_projection(self, ax, lattice='Z', color = 'blue'):
        # color code
        if lattice == 'X':
            selected_lattice = self.Xlattice_avg
        elif lattice == 'Y':
            selected_lattice = self.Ylattice_avg
        elif lattice == 'Z':
            selected_lattice = self.Zlattice_avg

        plot_pointproj(selected_lattice, ax, color)

    def add_trace(self, ax, lattice='Z', color = 'blue'):
        # color code
        if lattice == 'X':
            selected_lattice = self.Xlattice_avg
        elif lattice == 'Y':
            selected_lattice = self.Ylattice_avg
        elif lattice == 'Z':
            selected_lattice = self.Zlattice_avg

        plot_trace(selected_lattice, ax, color)

    # ACCESS DIFFRACTION PATTERNS

    def diffraction_match(self, fig, ax, n_diff, k_max = 1.2):
        bragg_peaks_fit = self.crystal_class.generate_diffraction_pattern(
            orientation_matrix=self.orientation_matrix[n_diff],
            sigma_excitation_error=0.01
        )

        py4DSTEM.process.diffraction.plot_diffraction_pattern(
            bragg_peaks_fit,
            bragg_peaks_compare=self.bragg_peaks[n_diff],
            scale_markers=3000,
            scale_markers_compare=1000,
            plot_range_kx_ky=[k_max, k_max],
            figsize=(4, 4),
            returnfig=False,  # Set returnfig to False because we're using input_fig_handle
            input_fig_handle=(fig, [ax]) 
        )
