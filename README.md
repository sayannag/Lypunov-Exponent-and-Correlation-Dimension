# Lypunov-Exponent-and-Correlation-Dimension
Calculates the largest Lyapunov exponent and the correlation dimension.


Python implementation of the Objective C code (https://www.physionet.org/content/lyapunov/1.0.0/l1d2/l1d2.c)


Reference: M.T. Rosenstein, J.J. Collins, C.J. De Luca, A practical method for calculating largest Lyapunov exponents from small data sets, Physica D 65:117-134, 1993.

If you use this reporsitory please cite as:

APA format:

Sayan, N. (2022). Lypunov-Exponent-and-Correlation-Dimension (Version 1.0.0) [Computer software]. https://github.com/sayannag/Lypunov-Exponent-and-Correlation-Dimension

BibTeX:

@software{Sayan_Lypunov-Exponent-and-Correlation-Dimension_2022,
author = {Sayan, Nag},
month = {3},
title = {{Lypunov-Exponent-and-Correlation-Dimension}},
url = {https://github.com/sayannag/Lypunov-Exponent-and-Correlation-Dimension},
version = {1.0.0},
year = {2022}
}

Usage (on Collab):

1. Write your own test script as:

```

!git clone https://github.com/sayannag/Lypunov-Exponent-and-Correlation-Dimension.git

%cd /content/Lypunov-Exponent-and-Correlation-Dimension/Lypunov-Exponent-and-Correlation-Dimension

import numpy as np

from math import log

import matplotlib.pyplot as plt

from lyapunov_Dc_utils import *

def main():

    t = np.linspace(1,100,1000)
    
    series_data = (1/t**2)*np.sin(t)

    LLE, Dc = Lyapunov().compute(series_data, debug_lyap_plot=False, debug_Dc_plot=False)

    print('Largest Lyapunov:', LLE, ', Dc: ', Dc)

if __name__ == "__main__":

    main()

```
2. Use the given test script as:

```

!git clone https://github.com/sayannag/Lypunov-Exponent-and-Correlation-Dimension.git

%cd /content/Lypunov-Exponent-and-Correlation-Dimension/Lypunov-Exponent-and-Correlation-Dimension

!python test_lyapunov_Dc.py

```
