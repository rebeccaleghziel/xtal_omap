import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt

def generate_ellipse(x0, y0, a, b, theta_deg):
    
    theta = np.deg2rad(theta_deg)
    t = np.linspace(0, 2 * np.pi, 300)
    
    # Parametric equations with rotation
    x = x0 + a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)
    y = y0 + a * np.cos(t) * np.sin(theta) + b * np.sin(t) * np.cos(theta)

    return x, y

def generate_mask(params_inner, params_outer, bf):
    
    height, width = bf.shape
    Y, X = np.ogrid[:height, :width]

    xi = X - params_inner[0]
    yi = Y - params_inner[1]
    cos_i = np.cos(np.deg2rad(params_inner[4]))
    sin_i = np.sin(np.deg2rad(params_inner[4]))
    inner_eq = ((xi * cos_i + yi * sin_i) / params_inner[2])**2 + ((xi * sin_i - yi * cos_i) / params_inner[3])**2
    inner_mask = inner_eq <= 1

    xo = X - params_outer[0]
    yo = Y - params_outer[1]
    cos_o = np.cos(np.deg2rad(params_outer[4]))
    sin_o = np.sin(np.deg2rad(params_outer[4]))
    outer_eq = ((xo * cos_o + yo * sin_o) / params_outer[2])**2 + ((xo * sin_o - yo * cos_o) / params_outer[3])**2
    outer_mask = outer_eq <= 1

    ring_mask = outer_mask & ~inner_mask
    y_coords, x_coords = np.where(ring_mask)


    return y_coords, x_coords
