# xtal_omap - Python classes for automated orientation mapping of 4D-STEM nanodiffraction data

xtal_omap provides processing classes for mapping crystal orientations from scanning nanodiffraction data, clustering connected
orientation components and reprojecting orientations into stereographic projections for 3D visualisation in real space. 



<bf>
<bf>
    
## Installation

You can use pip to install xtal_omap into your preferred environment.

- First open your preferred shell (Windows users may be using something like gitbash) 

- Activate the Python environment that you wish to use.

- Run

      python -m pip install git+https://github.com/LotharHouben/xtal_omap.git@main



## Test Your Installation

type “python” at the command prompt in your chosen terminal to start a Python session in your active Python environment.

You can now import , create an instance of the nbed class and dosplay the docstring for the LoadFile method:


    ➜ python
    Python 3.11.4 | packaged by conda-forge'
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import xtal_omap
    >>> a=Dataset.Dataset()
    >>> 

## Documentation

Please refer to the jupyter notebooks with examples in the directory 'examples' at 

    https://github.com/LotharHouben/xtal_omap/scripts/


The notebooks use a demonstration data set that is available under https://doi.org/10.5281/zenodo.15212905 
The nbed module is required for preprocessing the raw data and obtain a vector list of diffraction peak. nbed is available from

https://github.com/LotharHouben/nbed/