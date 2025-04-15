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
        self.CoM()

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

    def compute_projections(self):

        self.proj_a = stereoproj(self.Xlattice_avg)
        self.proj_b = stereoproj(self.Ylattice_avg)
        self.proj_c = stereoproj(self.Zlattice_avg)

    def CoM(self):
        y_coords = self.coordinates[:, 0]
        x_coords = self.coordinates[:, 1]

        mean_x = np.mean(x_coords)
        mean_y = np.mean(y_coords)

        self.CoM = (mean_x, mean_y)

    # Functions used in main methods
    
    def angle_with_z_axis(self, lattice_vec):
        
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

    def angle_with_x_axis(self, lattice_vec):
        
        vectors = lattice_vec
        vectors = vectors.reshape(-1, 3)
        projected_vector = np.column_stack((vectors[:, 0], vectors[:, 1]))
        x_axis = np.array([1, 0])
        
        # Magnitude of the projected vector and x-vector (=1)
        projected_vector_magnitude = np.linalg.norm(projected_vector, axis=1)
        x_axis_magnitude = np.linalg.norm(x_axis)
        
        # Dot product between the projected vector and the x-axis unit vector
        dot_product = np.dot(projected_vector, x_axis)
        
        # Obtain angles
        cos_theta = dot_product / (projected_vector_magnitude * x_axis_magnitude)
        theta_rad = np.arccos(cos_theta)
        theta_deg = np.degrees(theta_rad)
        
        average_theta = np.mean(theta_deg)
        std_theta = np.std(theta_deg)

        return average_theta, std_theta
 
    # COMPUTE ANGLES TO PUT INTO CLASS VARIABLES
    
    # def calculate_angles(self):
    #     self.average_thetaZz, self.std_thetaZz = self.angle_with_z_axis(self.Zlattice_vec)
    #     self.average_thetaYz, self.std_thetaYz = self.angle_with_z_axis(self.Ylattice_vec)
    #     self.average_thetaXz, self.std_thetaXz = self.angle_with_z_axis(self.Xlattice_vec)
    #     self.average_thetaZx, self.std_thetaZx = self.angle_with_x_axis(self.Zlattice_vec)
    #     self.average_thetaYx, self.std_thetaYx = self.angle_with_x_axis(self.Ylattice_vec)
    #     self.average_thetaXx, self.std_thetaXx = self.angle_with_x_axis(self.Xlattice_vec)
        
    
    # def reflect_180(self):
    #     self.Xlattice_avg = rot_180x_180y(self.Xlattice_avg)
    #     self.Ylattice_avg = rot_180x_180y(self.Ylattice_avg)
    #     self.Zlattice_avg = rot_180x_180y(self.Zlattice_avg)
    
    
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
