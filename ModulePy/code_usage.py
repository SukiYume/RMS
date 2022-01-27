import numpy as np
import pandas as pd
from extract_pulse import get_pulse
from rm_qu_fitting import mcmc_fit, plot_mcmc_samp
from rm_synthesis import synthesis_rm, get_rm, plot_rm_synthesis

if __name__ == '__main__':
    
    ##### Read File; Format File
    data = pd.read_csv('../CalData/RM-190417-7.txt', sep='\s+', header=None)
    data = pd.read_csv('../CalData/RM-190303-1.txt', sep='\s+', header=None)
    freq = np.linspace(1000, 1500, data.loc[:, 1].max()+1)
    
    ##### I, Q, U, V - 2D numpy array; axis0 is frequency channel, axis1 is time sample
    I = data.iloc[:, 3].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    Q = data.iloc[:, 4].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    U = data.iloc[:, 5].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    V = data.iloc[:, 6].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    
    ##### Extract pulse; 
    I, Q, U, _, freq, snr = get_pulse(I, Q, U, V, freq)
    ##### RM synthesis
    rm_list, Linear = synthesis_rm(I, Q, U, freq, rm_left=-10000, rm_right=10000)
    RM, RM_error_left, RM_error_right = get_rm(rm_list, Linear, snr)
    print('RM: {:.0f} +{:.0f} -{:.0f}'.format(RM, RM_error_right, RM_error_left))
    plot_rm_synthesis(rm_list, Linear, save=False)
    
    ##### Squeeze the time dimension and turn the data into one dimension
    Q, U = np.mean(Q, axis=1), np.mean(U, axis=1)
    ##### QU fitting with MCMC, Q and U - 1D numpy array; axis0 is frequency channel
    result, RM, RM_error_left, RM_error_right = mcmc_fit(Q, U, freq, rm_left=-10000, rm_right=10000)
    print('RM: {:.0f} +{:.0f} -{:.0f}'.format(RM, RM_error_right, RM_error_left))
    plot_mcmc_samp(result, save=False)