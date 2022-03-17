import numpy as np
from scipy import interpolate
from astropy import constants as const

import seaborn as sns
import matplotlib.pyplot as plt

plt.style.use('default')
sns.set_color_codes()

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.cal'] = 'Arial'
rcParams['mathtext.it'] = 'Arial'
rcParams['mathtext.rm'] = 'Arial'

from extract_pulse import get_time_index_new

def synthesis_rm(I, Q, U, freq, rm_left=-10000, rm_right=10000):
    rm_num     = int((rm_right - rm_left) / 20)
    wave       = const.c.value / freq / 1e6
    wave       = wave.reshape(-1, 1)
    Linear, rm_list = [], np.linspace(rm_left, rm_right, rm_num)
    for RM in rm_list:
        PA     =   2 * RM * wave**2
        Q_C    =   np.cos(PA) * Q + np.sin(PA) * U
        U_C    = - np.sin(PA) * Q + np.cos(PA) * U
        Lsum   =   np.sum(np.sqrt(np.sum(Q_C, axis=0)**2 + np.sum(U_C, axis=0)**2))
        Linear =   np.append(Linear, Lsum/np.sum(I))
    return rm_list, Linear

# def get_rm(rm_list, Linear, snr):
#     rm_peak_range_index = get_time_index_new(Linear, squeeze_frac=10)
#     x, y = rm_list[rm_peak_range_index!=0], Linear[rm_peak_range_index!=0]
#     f = interpolate.interp1d(x, y, kind='cubic')
#     x_new = np.linspace(x.min(), x.max(), 500)
#     y_new = f(x_new)
#     rm_max = x_new[np.argmax(y_new)]
#     half_rm = x_new[y_new > (y_new.max() - y_new.min()) / 2 + y_new.min()]
#     rm_error_left, rm_error_right = rm_max - half_rm.min(), half_rm.max() - rm_max
#     return rm_max, rm_error_left / snr, rm_error_right / snr

def get_rm(rm_list, Linear, snr):
    rm_max         = rm_list[np.argmax(Linear)]
    rm_index       = np.argmax(Linear) - np.where(Linear<(np.max(Linear)+np.min(Linear))/2)
    rm_error_left  = rm_max - rm_list[np.argmax(Linear) - rm_index[rm_index>0].min()]
    rm_error_right = rm_list[np.argmax(Linear) - rm_index[rm_index<0].max()] - rm_max
    return rm_max, rm_error_left / snr, rm_error_right / snr

def plot_rm_synthesis(rm_list, Linear, save):
    
    plt.figure(figsize=(5, 4))
    plt.plot(rm_list, Linear*100, color='royalblue', alpha=0.7)
    plt.xlabel('RM (rad/m$^2$)')
    plt.ylabel('Degree of Linear Polarization (%)')
    
    if save:
        plt.savefig('../Figure/RM-Synthesis-{}.png'.format(save), format='png', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()