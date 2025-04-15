# header hooks

_PROGRAM_NAME = 'xtal_omap'
_USER_NAME = ''

# helper functions

from .stereofunctions import (
    initialize_stereogram,   
    stereoproj,
    compute_trace,
    angles_theta_phi,
    plot_pointproj,
    plot_trace,
    angle_between_vectors,
    angle_vectors_360,
    scalevector
)

# classes                                                                                             

from xtal_omap.classes import (
    CrystalGroup,
    Dataset
)


                                                                                      
