import numpy as np
from math import log
import matplotlib.pyplot as plt
from lyapunov_Dc_utils import *

import argparse

def main(args):
    # series_data = [0.01]
    # for i in range(1000):
    #   series_data.append(series_data[i-1]*3.5*(1-series_data[i-1]))
    
    t = np.linspace(1,100,1000)
    
    if args.series_type == 'oscillating':
        series_data = np.sin(t)
        
    elif args.series_type == 'converging':
        series_data = (1/t**2)*np.sin(t)
        
    elif args.series_type == 'diverging':
        series_data = (t**2)*np.sin(t)

    LLE, Dc = Lyapunov().compute(series_data, debug_lyap_plot=False, debug_Dc_plot=False)

    print('Series Type: '+args.series_type,', Largest Lyapunov:', LLE, ', Dc: ', Dc)

if __name__ == "__main__":
    parser = argparse.ArgumentParser('Lyapunov_Dc')
    parser.add_argument('--series_type', type=str, default='converging', choices=['converging', 'diverging', 'oscillating'])
    args = parser.parse_args()
    main(args)
