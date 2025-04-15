# Written by R. Leghziel
# This module is used to obtain stereographic projections of vectors
# 1. Construct the Wulff net which consists in 3D circles delimiting a sphere at an interval angle of 10 deg and projected in 2D
# 2. Project the vectors intersecting a unit sphere
# 3. Project the trace of a plane given it's normal vector

import numpy as np
import matplotlib.pyplot as plt

# Wulff Net Construction

def initialize_stereogram(ax, style = 'default'):

    plt.style.use([style])
   
    if style == 'dark_background':
        color_maincircle = 'white'
        color_incircles = 'white'
    else:
        color_maincircle = 'black'
        color_incircles = 'gray'

    r = 1
    res_circles = 75
    res_smallcircles = 25

    # Define figure

    # DEFINE OUTSIDE
    # fig = plt.figure()
    # plt.figure(figsize=(6, 6))
    
    # LONGITUDES
    
    # Main circle
    x = np.linspace(0, 2*np.pi, res_circles)
    ax.plot(r*np.cos(x), r*np.sin(x), color = color_maincircle, alpha = 1, linewidth = 0.9)
    
    # Projection of central circle - a line
    xc_long = [0] * 2 # need just 2 points for a line
    yc_long = (np.linspace(-r, r, 2))
    ax.plot(xc_long, yc_long, color=color_incircles, alpha=1, linewidth=0.5, linestyle = "--")
    
    # PLOT GREAT CIRCLES
    range1 = np.linspace(0, 2*np.pi, res_circles)
    deg_long = 10 # degrees between 3D circles to project
    
    # start from 10 deg because do not project main circle, and finish with 80 deg because don't project central circle
    increment = np.linspace(np.radians(deg_long), np.pi/2 - np.radians(deg_long), 8)
    great_circles = np.empty((0, 2))  # Initialize an empty numpy array with shape (0, 2)
    
    i = 0
    for theta in increment:
        for t in range1:
            # Rotation matrix across y-axis
            Ry = np.array([
                [np.cos(theta), 0, np.sin(theta)],
                [0, 1, 0],
                [-np.sin(theta), 0, np.cos(theta)],
            ])
    
            circle_vector_long = np.array([r*np.cos(t), r*np.sin(t), 0]) # lying in xy plane
            point_rot = np.dot(Ry, circle_vector_long)
    
            # Projection
            if point_rot[2] == 0:
                point = np.array([point_rot[0], point_rot[1]])        
            else:
                point = np.array([point_rot[0] / (np.abs(point_rot[2]) + 1), point_rot[1] / (np.abs(point_rot[2]) + 1)])
            
            great_circles = np.vstack((great_circles, point))  # Append the point to great_circles
    
        ax.plot(great_circles[i*res_circles:(i+1)*res_circles, 0], great_circles[i*res_circles:(i+1)*res_circles, 1], color=color_incircles, alpha=1, linewidth=0.5, linestyle = "--")
        i += 1
    
    small_circles = np.empty((0, 2))
    
    psi_range = np.linspace(-np.pi/2+np.radians(10), np.pi/2-np.radians(10), 17)
    range2 = np.linspace(0, np.pi, res_smallcircles)
    
    j=0
    for psi in psi_range:
        for t in range2:
    
            r_sc = r*np.cos(psi) #radius of small circle
            cy = r*np.sin(psi)
            
            x = r_sc*np.cos(t)
            y = cy
            z = r_sc*np.sin(t)
            smallcircle_vector = np.array([x, y, z]) # lying in xy plane
    
            # Projection
            point = np.array([smallcircle_vector[0] / (np.abs(smallcircle_vector[2]) + 1), smallcircle_vector[1] / (np.abs(smallcircle_vector[2]) + 1)])        
            small_circles = np.vstack((small_circles, point))  # Append the point to great_circles
    
        ax.plot(small_circles[j*res_smallcircles:(j+1)*res_smallcircles, 0], small_circles[j*res_smallcircles:(j+1)*res_smallcircles, 1], color=color_incircles , alpha=1, linewidth=0.5, linestyle = '--')
        j += 1
    
    # # Add labels
    
    ax.text(np.cos(np.radians(60))-0.01, np.sin(np.radians(60))+0.05, str(30))
    ax.text(np.cos(np.radians(60))-0.01, -np.sin(np.radians(60))-0.1, str(30))
    ax.text(-np.cos(np.radians(60))-0.1, -np.sin(np.radians(60))-0.1, str(30))
    ax.text(-np.cos(np.radians(60))-0.1, np.sin(np.radians(60))+0.03, str(30))
    
    
    ax.text(np.cos(np.radians(30))+0.03, np.sin(np.radians(30))+0.03, str(60))
    ax.text(-np.cos(np.radians(30))-0.1, np.sin(np.radians(30))+0.03, str(60))
    ax.text(np.cos(np.radians(30))+0.03, -np.sin(np.radians(30)+0.1), str(60))
    ax.text(-np.cos(np.radians(30))-0.1, -np.sin(np.radians(30)+0.1), str(60))
    
    ax.text(0-0.025, 1+0.05, str(0))
    ax.text(0-0.025, -1-0.1, str(0))
    ax.text(1+0.025, 0-0.025, str(90))
    ax.text(-1-0.12, 0-0.025, str(90))
    
    # PLOT
    ax.set_aspect('equal')
    ax.axis('off')

    return ax

# Compute Stereographic Projection of a vector on the Wulff Net

def stereoproj(vector):

    V = scalevector(vector, 1)
    
    if V[2] < 0:
        vect_to_proj = np.array([V[0], V[1], -V[-1]])
        point = np.array([vect_to_proj[0] / (vect_to_proj[2] + 1), vect_to_proj[1] / (vect_to_proj[2] + 1)])
        hemisphere = False

    elif V[2] == 0:
        point = np.array([V[0], V[1]])
        hemisphere = True

    else:
        vect_to_proj = np.array([V[0], V[1], V[-1]])
        point = np.array([vect_to_proj[0] / (vect_to_proj[2] + 1), vect_to_proj[1] / (vect_to_proj[2] + 1)])
        hemisphere = True
    
    return point, hemisphere

# Compute Trace of a plane given it's normal vector

def compute_trace(normal_vector):
      
    theta_d, phi_d = angles_theta_phi(normal_vector)
    
    phi = np.radians(phi_d)
    theta = np.radians(theta_d)

    Rz = np.array([
        [np.cos(phi), np.sin(phi), 0],
        [-np.sin(phi), np.cos(phi), 0],
        [0, 0, 1]
    ])

    # make sure the - sin is correct / changed it to fit the z-axis 

    Ry = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)],
    ])

    Rzy = np.dot(Rz, Ry)

    range = np.linspace(0, 2*np.pi, 1000) # change number to 100
    trace_p = []
    trace_n = []

    for t in range:
        circle_vector = [np.cos(t), np.sin(t), 0]
        point_rot = np.dot(Rzy, circle_vector) # Apply rotation vector

        if point_rot[2]> 0:
            point = [point_rot[0] / (point_rot[2] + 1), point_rot[1] / (point_rot[2] + 1)]
            trace_p.append(point)

        elif point_rot[2]< 0:
            point = [point_rot[0] / (-point_rot[2] + 1), point_rot[1] / (-point_rot[2] + 1)]
            trace_n.append(point)

        else: # Case where ==0
            point = [point_rot[0], point_rot[1]]
            # vect_to_proj = np.array([point_rot[0], point_rot[1], point_rot[-1]])
            # point = [vect_to_proj[0] / (vect_to_proj[2] + 1), vect_to_proj[1] / (vect_to_proj[2] + 1)]
            # trace_p.append(point)
            # trace_n.append(point)


    trace_p_array = np.array(trace_p)
    trace_n_array = np.array(trace_n)

    # sort the arrays
    
    P_sorted_indices = np.argsort(trace_p_array[:, 1])
    P_sorted = trace_p_array[P_sorted_indices]

    N_sorted_indices = np.argsort(trace_n_array[:, 1])
    N_sorted = trace_n_array[N_sorted_indices]
    
    return P_sorted, N_sorted

# function to obtain polar and azimuthal angles

def angles_theta_phi(vector):
    
    z = np.array([0, 0, 1])
    x = np.array([1, 0, 0])

    vector_xy = np.array([vector[0], vector[1], 0])
    
    theta_p = angle_between_vectors(vector, z) # polar angle
    phi_a = angle_vectors_360(vector_xy, x) # azimuthal angle clockwise

    return theta_p, phi_a
    
# Plotting Functions: plot proj, plot trace

def plot_pointproj(vector, ax, color):
    
    point, hemisphere = stereoproj(vector)
    
    if hemisphere is True:
        ax.scatter(point[0], point[1], s=50, marker='s', facecolors=color)
    else:
        ax.scatter(point[0], point[1], s=50,  marker='s', facecolors='none', edgecolors=color, linewidths=3)

def plot_trace(normal_vector, ax, color):

    P, N = compute_trace(normal_vector)
    ax.scatter(P[:, 0], P[:, 1], color = color, alpha = 0.2, s=1)
    ax.scatter(N[:, 0], N[:, 1], color = color, alpha = 0.1, s=1)


# Computation functions

def angle_between_vectors(a, b):
    dot_product = np.dot(a, b)
    magnitude_a = np.linalg.norm(a)
    magnitude_b = np.linalg.norm(b)
    cos_theta = dot_product / (magnitude_a * magnitude_b)
    theta_radians = np.arccos(cos_theta) # Calculate the angle in radians
    theta_degrees = np.degrees(theta_radians) # Convert the angle to degrees
    return theta_degrees

# clockwise angle between two vectors
def angle_vectors_360(a, b):
    dot_product = np.dot(a, b)
    magnitude_a = np.linalg.norm(a)
    magnitude_b = np.linalg.norm(b)
    cos_theta = dot_product / (magnitude_a * magnitude_b)
    theta_radians = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    cross_product = a[0] * b[1] - a[1] * b[0]

    # Adjust the angle to 0-360 degrees
    if cross_product < 0:  # Clockwise
        theta_radians = 2 * np.pi - theta_radians

    theta_d = np.degrees(theta_radians)  # Convert to degrees
    return theta_d

def scalevector(vector, length): # returns a vector scaled to length but maintaining the same direction
    scale_factor = length/np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    scaled_vector = scale_factor*vector
    return scaled_vector
