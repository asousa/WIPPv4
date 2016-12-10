import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm


# --------------- Latex Plot Beautification --------------------------
fig_width_pt = 650.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width+1,fig_height+1]
params = {'backend': 'ps',
          'axes.labelsize': 14,
          'font.size': 14,
          'legend.fontsize': 10,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'text.usetex': False,
          'figure.figsize': fig_size}
plt.rcParams.update(params)
# --------------- Latex Plot Beautification --------------------------


# ----------------------------------------------------------------------------
# Plot a single time-energy deflection matrix
# ----------------------------------------------------------------------------

def plot_pN_pS(pN, pS, sc):
    tvec = np.linspace(sc.T_STEP,sc.T_MAX,sc.NUM_STEPS)

    clims = [-7, 0]


    pN_P = np.log10(pN)
    np.clip(pN_P,clims[0],clims[1],out=pN_P)

    pS_P = np.log10(pS)
    np.clip(pS_P,clims[0],clims[1],out=pS_P)

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    # The goods
    p1 = ax1.imshow(pN_P,origin='lower',aspect='auto')
    p1.set_clim(clims)
    p2 = ax2.imshow(pS_P,origin='lower',aspect='auto')
    p2.set_clim(clims)

    # Colorbar
#    cax = fig.add_axes([0.77, 0.12, 0.025, 0.78])
    cax = fig.add_axes([0.92, 0.12, 0.025, 0.78])
    fig.colorbar(p1,cax=cax)

    # Label axes
    
    tlabels = np.arange(1,np.floor(sc.T_MAX))
    tinds =  [np.argmax(tt <= tvec) for tt in tlabels]
    ax2.set_xticks(tinds)
    tlabel_strings = ['%d'%k for k in tlabels]
    ax2.set_xticklabels(tlabel_strings)
    

    ax1.get_xaxis().set_visible(False)
#     ax2.set_xticks(np.floor(np.linspace(0,sc.NUM_STEPS-1,sc.T_MAX + 1)))
#     ax2.set_xticklabels(np.floor(tvec[ax2.get_xticks().astype(int)]))

    # Label each power of 10
    logvals = np.arange(np.log10(sc.E_MIN), np.log10(sc.E_MAX)+1)
    einds =  [np.argmax(lv <= np.log10(sc.E_tot_arr)) for lv in logvals]
    einds[-1] = sc.NUM_E - 1 

    ylabel_strings = ['$10^%d$'%k for k in logvals]
    ax1.set_yticks(einds)
    ax1.set_yticklabels(ylabel_strings)
    ax2.set_yticks(einds)
    ax2.set_yticklabels(ylabel_strings)
    ax2.set_xlabel('Time (sec)')
    ax1.set_ylabel('Energy (eV)')
    ax2.set_ylabel('Energy (eV)')

#     fig.subplots_adjust(hspace=0.03, wspace=0.05)
    fig.canvas.draw()
    # plt.show()

    # print np.max(pN_P)
    # print np.min(pN_P)
    # print np.max(pS)
    # print sum(sum((pN!=0)))



