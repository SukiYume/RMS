import emcee
import numpy as np
import pandas as pd
from astropy import constants as const

import seaborn as sns
from scipy import interpolate
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

def gridcal_rm_pa(Q, U, freq, rm_left=-20000, rm_right=20000):
    ########## 粗搜 ##########
    rm_num      = int((rm_right - rm_left) / 20)
    RM          = np.linspace(rm_left, rm_right, rm_num)
    PAI         = np.deg2rad(np.linspace(-90, 90, 100))
    xx, yy      = np.meshgrid(RM, PAI)
    wave        = const.c.value / freq / 1e6
    L           = np.sqrt(Q**2 + U**2)
    a           = 1 / np.sum((Q / L - np.cos(2 * yy[:,:,np.newaxis] + 2 * xx[:,:,np.newaxis] * wave**2))**2 + \
                             (U / L - np.sin(2 * yy[:,:,np.newaxis] + 2 * xx[:,:,np.newaxis] * wave**2))**2, axis=-1)
    
    xindex      = np.where(a==np.max(a))[0][0]
    yindex      = np.where(a==np.max(a))[1][0]
    
    ########## 细搜 ##########
    rm_tmp_max  = RM[yindex]
    rm_right, rm_left, rm_num = rm_tmp_max + 2000, rm_tmp_max - 2000, 1000
    RM          = np.linspace(rm_left, rm_right, rm_num)
    PAI         = np.deg2rad(np.linspace(-90, 90, 100))
    xx, yy      = np.meshgrid(RM, PAI)
    a           = 1 / np.sum((Q / L - np.cos(2 * yy[:,:,np.newaxis] + 2 * xx[:,:,np.newaxis] * wave**2))**2 + \
                             (U / L - np.sin(2 * yy[:,:,np.newaxis] + 2 * xx[:,:,np.newaxis] * wave**2))**2, axis=-1)
    
    xindex      = np.where(a==np.max(a))[0][0]
    yindex      = np.where(a==np.max(a))[1][0]
    
    rm_max_line = a[xindex]
    intindex    = (np.diff(np.sign(np.diff(rm_max_line))) > 0).nonzero()[0] + 1
    intindex    = intindex[np.abs(np.argmax(rm_max_line) - intindex) > 5]
    intindex    = np.sort(intindex[np.argsort(np.abs(np.argmax(rm_max_line) - intindex))[:2]])
    x, y        = RM[intindex[0]: intindex[1]], rm_max_line[intindex[0]: intindex[1]]
    f           = interpolate.interp1d(x, y, kind='cubic')
    x_new       = np.linspace(x.min(), x.max(), 500)
    y_new       = f(x_new)
    rm_max      = x_new[np.argmax(y_new)]
    
    return rm_max, np.degrees(yy[xindex][yindex])

def log_likelihood(par, wave, y1, y2):
    pa, rm       = par
    p            = 2 * np.deg2rad(pa) + 2 * rm * wave**2
    calul, calql = np.sin(p), np.cos(p)
    return -0.5 * np.sum((calul-y1)**2 + (calql-y2)**2)

def log_probability(par, wave, y1, y2):
    pa, rm = par
    if -90 < pa < 90 and -20000 < rm < 20000:
        lp =  0
    else:
        lp = -np.inf
    lk     = log_likelihood(par, wave, y1, y2)
    if np.isnan(lk):
        return -np.inf
    return lp + lk

def mcmc_fit(Q, U, freq, rm_left=-20000, rm_right=20000):
    wave                  = const.c.value / freq / 1e6
    L                     = np.sqrt(Q**2 + U**2)
    approx_rm, approx_pa  = gridcal_rm_pa(Q, U, freq, rm_left, rm_right)
    ndim, nwalkers        = 2, 10
    p0 = np.zeros((nwalkers, ndim))
    p0[:, 0]              = np.random.rand(nwalkers) * 20 + approx_pa
    p0[:, 1]              = np.random.rand(nwalkers) * 20 + approx_rm
    
    sampler               = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=[wave, U/L, Q/L])
    pos                   = sampler.run_mcmc(p0, 3000, progress=True)
    sampler.reset()
    sampler.run_mcmc(pos, 1000, progress=True)
    result                = sampler.chain[:, :, :].reshape((-1, ndim))
    
    lim_value             = np.sort(result[:, 1])[[int(result.shape[0]*0.16), int(result.shape[0]*0.5), int(result.shape[0]*0.84)]]
    lowlim, uplim         = np.diff(lim_value)
    return result, np.round(lim_value[1], 0), np.round(lowlim, 0), np.round(uplim, 0)

def plot_mcmc_samp(result, save=False):
    
    plt.figure(figsize=(6, 5))
    plt.subplots_adjust(wspace=0, hspace=0)
    gs = gridspec.GridSpec(6, 7)

    lim_value = np.sort(result[:, 1])[[int(result.shape[0]*0.16), int(result.shape[0]*0.5), int(result.shape[0]*0.84)]]
    lowlim, uplim = np.diff(lim_value)
    samp_data = pd.DataFrame({'PA': result[:, 0], 'RM': result[:, 1]})
    pa_min, pa_max = np.round(np.sort(result[:, 0])[[int(len(result)*0.01), int(len(result)*0.99)]], 0)
    rm_min, rm_max = np.round(np.sort(result[:, 1])[[int(len(result)*0.01), int(len(result)*0.99)]], 0)

    ax1 = plt.subplot(gs[0:1, 0:5])
    g = sns.kdeplot(data=samp_data, x='PA', fill=False, ax=ax1, color='royalblue')
    sns.despine(right=True, left=True, top=True, ax=ax1)
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks([])
    plt.yticks([])
    plt.xlim(pa_min, pa_max)
    ax = plt.gca()
    plt.text(0.93, 0.3, 'A', weight='bold', transform=ax.transAxes)

    ax1 = plt.subplot(gs[1:6, 0:5])
    g = sns.kdeplot(data=samp_data, x='PA', y='RM', levels=6, cmap=new_cmap, fill=True, ax=ax1, 
                    label='RM: {:.0f}+{:.0f}-{:.0f}'.format(np.round(lim_value[1], 0), np.round(uplim, 0), np.round(lowlim, 0)),
                    cbar=True, cbar_kws={'label': 'Number Density'})
    plt.xlim(pa_min, pa_max)
    plt.ylim(rm_min, rm_max)
    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
    pos_joint_ax = g.get_position()
    g.set_position([pos_joint_ax.x0, pos_joint_ax.y0, 0.5, pos_joint_ax.height])
    g.figure.axes[-1].set_position([.73, pos_joint_ax.y0, .07, pos_joint_ax.height])

    cbar_ticks = g.figure.axes[-1].get_yticks()
    _, cbar_max = g.figure.axes[-1].get_ylim()
    g.figure.axes[-1].set_yticks(cbar_ticks, ['{} %'.format(np.int64(np.round(t/cbar_max*100, 0))) for t in cbar_ticks])

    g.set_xlabel('PA (degree)')
    g.set_ylabel('RM (rad/m$^2$)')
    ax = plt.gca()
    plt.text(0.93, 0.93, 'B', weight='bold', transform=ax.transAxes)

    ax1 = plt.subplot(gs[1:6, 5:6])
    g = sns.kdeplot(data=samp_data, y='RM', fill=False, ax=ax1, color='royalblue')
    sns.despine(right=True, bottom=True, top=True, ax=ax1)
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks([])
    plt.yticks([])
    plt.ylim(rm_min, rm_max)
    ax = plt.gca()
    plt.text(0.3, 0.93, 'C', weight='bold', transform=ax.transAxes)
    if save:
        plt.savefig('./Figure/RM-QU-Fitting-{}.png'.format(save), format='png', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()