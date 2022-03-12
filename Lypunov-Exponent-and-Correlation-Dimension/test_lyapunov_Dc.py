import numpy as np
from math import log
import matplotlib.pyplot as plt
from lyapunov_Dc_utils import *

def main():
    # series_data = [0.01]
    # for i in range(1000):
    #   series_data.append(series_data[i-1]*3.5*(1-series_data[i-1]))

    t = np.linspace(1,100,1000)
    series_data = (1/t**2)*np.sin(t)

    LLE, Dc = Lyapunov().compute(series_data, debug_lyap_plot=False, debug_Dc_plot=False)

    print('Largest Lyapunov:', LLE, ', Dc: ', Dc)

if __name__ == "__main__":
    main()