import os, re
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import constants as const

import seaborn as sns
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
clist = ['white', 'royalblue']
new_cmap = LinearSegmentedColormap.from_list('wb', clist)

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

def cal_linear_circular(profile_I, L, profile_V):
    I_index      = get_time_index_new(profile_I)
    I, L, V      = profile_I[I_index==1].sum(), L[I_index==1].sum(), profile_V[I_index==1].sum()
    rms          = profile_I[:20].std()
    snr          = I / np.sqrt(len(profile_I[I_index==1])) / rms
    linear       = L / I * 100
    linear_err   = np.sqrt(1 + L**2 / I**2) / snr * 100
    circular     = V / I * 100
    circular_err = np.sqrt(1 + V**2 / I**2) / snr * 100
    return linear, linear_err, circular, circular_err

def translate_file(file_name):

    with fits.open('../Data/' + file_name) as f:
        tbin     = f[2].header['TBIN']
        nbin     = f[2].header['NBIN']
        freq     = f[2].data['DAT_FREQ'][0, :]
    
    if not os.path.exists('../CalData/' + file_name + '.txt'):
        with open('../Tfile/' + file_name + '.txt', 'r') as f:
            data = f.read()
        rfi_mask = re.findall(r'(subint=0 chan=(\d+) zapped)', data)
        for i in rfi_mask:
            data = re.sub(i[0], '\n'.join(['0 {} {} 0 0 0 0 0 0'.format(i[1], j) for j in range(nbin)]), data)
        with open('../CalData/' + file_name + '.txt', 'w') as f:
            f.write(data)
    
    data         = pd.read_csv('../CalData/' + file_name + '.txt', sep='\s+', header=None, skiprows=1)
    I            = data.iloc[:, 3].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    Q            = data.iloc[:, 4].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    U            = data.iloc[:, 5].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    V            = data.iloc[:, 6].values.reshape(data.loc[:, 1].max()+1, data.loc[:, 2].max()+1)
    
    tbin         = tbin * 4
    return I, Q, U, V, freq, tbin

def fix_rm(I, Q, U, V, freq, rm):
    wave       = const.c.value/freq/1e6
    wave       = wave.reshape(-1, 1)
    PA         = 2 * rm * wave**2
    Q_C        = np.cos(PA)*Q + np.sin(PA)*U
    U_C        = -np.sin(PA)*Q + np.cos(PA)*U
    I, Q, U, V = I, Q_C, U_C, V
    return I, Q, U, V

def plot_spec(file_name, rm_max, save=False):
    
    ##### Plot-Spectrum
    I, Q, U, V, freq, tbin = translate_file(file_name)
    I, Q, U, V             = fix_rm(I, Q, U, V, freq, rm=rm_max)
    
    fig = plt.figure(figsize=(5, 6))
    gs  = gridspec.GridSpec(11, 1)
    plt.subplots_adjust(hspace=0)
    
    profile_I, profile_Q, profile_U, profile_V = np.mean(I, axis=0), np.mean(Q, axis=0), np.mean(U, axis=0), np.mean(V, axis=0)
    normfactor, rms        = profile_I.max(), profile_I[:20].std()
    
    ax = plt.subplot(gs[0:3, 0])
    sigma_Q, sigma_U     = profile_Q[:10].std(), profile_U[:10].std()
    PA                   = np.degrees(np.arctan2(profile_U, profile_Q)/2)
    PAE                  = np.degrees(np.sqrt((profile_Q**2 * sigma_U**2 + profile_U**2 * sigma_Q**2)/(4 * (profile_Q**2 + profile_U**2)**2)))
    PAT                  = np.arange(len(PA)).astype(np.float64)
    pulse_index          = get_time_index_new(profile_I, squeeze_frac=3)
    thres                = np.std(PAE[pulse_index==1])
    PA[pulse_index==0]   = np.nan
    PAE[pulse_index==0]  = np.nan
    PA[PAE>thres*3]      = np.nan
    PAE[PAE>thres*3]     = np.nan
    plt.errorbar(PAT, PA, PAE, color='r', fmt='.', capsize=3, lw=1, ms=1)
    plt.hlines(90, 0, len(profile_I), ls='--', color='gray', alpha=0.5, lw=0.5)
    plt.hlines(-90, 0, len(profile_I), ls='--', color='gray', alpha=0.5, lw=0.5)
    plt.xticks([])
    plt.xlim(0, len(profile_I))
    plt.ylim(-120, 120)
    plt.ylabel('PA (degree)')
    ax = plt.gca()
    plt.text(0.95, 0.85, 'A', weight='bold', transform=ax.transAxes)
    
    plt.subplot(gs[3:6, 0])
    plt.plot(profile_I / normfactor, color='gray', alpha=0.8)
    L = np.sqrt(np.mean(Q, axis=0)**2 + np.mean(U, axis=0)**2)
    L[L/rms <= 1.57]     = 0
    L[L/rms > 1.57]      = np.sqrt(L[L/rms > 1.57]**2 - rms**2)
    plt.plot(L/normfactor, color='r', alpha=0.8)
    plt.plot(profile_V / normfactor, color='royalblue', alpha=0.8)
    plt.xlim(0, len(L))
    plt.xticks([])
    plt.ylabel('Normalized Intensity', labelpad=13)
    ax = plt.gca()
    plt.text(0.95, 0.85, 'B', weight='bold', transform=ax.transAxes)
    
    ax = plt.subplot(gs[6:, 0])
    x, y                 = np.meshgrid(np.arange(I.shape[1]) * tbin * 1e3, freq)
    vmin, vmax           = np.sort(I.flatten())[int(I.shape[0]*I.shape[1]/50)], np.sort(I.flatten())[int(I.shape[0]*I.shape[1]/50*49)]
    c = plt.pcolor(x, y, I, cmap='mako', vmin=vmin, vmax=vmax)
    plt.xlabel('Time (ms)')
    plt.ylabel('Frequence (MHz)', labelpad=6)
    plt.colorbar(cax=plt.axes([0.92, 0.11, 0.025, 0.35]), label='Intensity (arbitrary)')
    ax = plt.gca()
    plt.text(-2.35, 0.9, 'C', weight='bold', color='white', fontsize=11, transform=ax.transAxes)
    fig.align_labels()
    
    linear, linear_err, circular, circular_err = cal_linear_circular(profile_I, L, profile_V)
    print('Linear: {:.2f}_{:.2f}    Circular: {:.2f}_{:.2f}'.format(linear, linear_err, circular, circular_err))
    
    if save:
        plt.savefig('../Figure/Spec-{}-Linear_{:.2f}_{:.2f}-Circular_{:.2f}_{:.2f}.png'.format(file_name, linear, linear_err, circular, circular_err),
                   format='png', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()