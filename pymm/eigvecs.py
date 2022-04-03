import logging
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def get_ext(filename):
    '''
    '''

    p = re.compile('\.([a-zA-Z]{3})')

    match = p.search(filename)

    try:
        ext = match.group()[1:]
    except:
        logging.info('File extention was not recognized. Falling back to png.')
        ext = 'png'

    return ext

def eig_corr(eig, l, m, state, save_fn, dpi, bins):
    '''
    '''

    fig, ax = plt.subplots(1, 1)

    df = pd.DataFrame(eig[:,:,state], columns = ['Unperturbed state {}'.format(i) for i in range(eig.shape[1])])

    sns.histplot(df, x='Unperturbed state {}'.format(l),
                 y='Unperturbed state {}'.format(m), bins=bins,
                 binrange=(-1,1), cmap='flare',
                 ax=ax, cbar=True, cbar_kws={'label': 'Number of frames'})

    # cb = fig.colorbar(cset, ax=axs.ravel().tolist(), fraction=0.05)

    ax.set_title('Perturbed state {}'.format(state))

    ax.set_xlabel = r'$c_{{}}^{{}}$'.format(l, state)
    ax.set_ylabel = r'$c_{{}}^{{}}$'.format(m, state)

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    plt.tight_layout()

    ext = get_ext(save_fn)
    plt.savefig(save_fn[:-4] + '.' + ext, dpi=dpi, format=ext)
    logging.info('\n * Correlation plot saved as {}.\n'.format(save_fn[:-4] + '.' + ext))

    plt.show()

def eig_corr_tot(eig, state, save_fn, dpi, bins):
    '''
    '''

    plt.rcParams['font.size'] = 22
    plt.rcParams['legend.fontsize'] = 20
    plt.rcParams['axes.linewidth'] = 1.5


    df = pd.DataFrame(eig[:,:,state], columns = ['state {}'.format(i) for i in range(eig.shape[1])])
    g = sns.PairGrid(df, despine=False)
    g.map_diag(sns.kdeplot,
               color='k'
               )
    g.map_offdiag(sns.histplot,
                  bins=bins,
                  binrange=(-1,1),
                  #sns.kdeplot,
                  #fill=True,
                  cmap='flare'
                  )

    g.set(xlim=[-1,1], ylim=[-1,1])
    g.fig.subplots_adjust(top=0.95)
    g.fig.suptitle('Perturbed state {}'.format(state))

    ext = get_ext(save_fn)


    plt.savefig(save_fn[:-4] + '.' + ext, dpi=dpi, format=ext)
    logging.info('\n * Correlation plot between all states'
                 ' saved as {}.\n'.format(save_fn[:-4] + '.' + ext))

    plt.show()

def eig_hist(eig, n_states, save_fn, dpi):
    '''
    '''
    
    plt.rcParams['font.size'] = 22
    plt.rcParams['legend.fontsize'] = 20
    plt.rcParams['axes.linewidth'] = 1.5

    plt.style.use('seaborn-colorblind')

    xsize =  6 + n_states*0.7
    fig, ax = plt.subplots(1, 1, figsize=(xsize, 7))

    tot_states = eig.shape[1]
    if not n_states:
        n_states = tot_states

    ys = [[] for i in range(n_states)]
    # x = [i+1 for i in range(n_states)]

    for i in range(n_states):
        for j in range(n_states):
            ys[i].append((eig[:,j,i]**2).mean())

        other_sum = 0
        for j in range(n_states, tot_states):
            other_sum += (eig[:,j,i]**2).mean()
        ys[i].append(other_sum)
    ys = np.array(ys)
    df = pd.DataFrame(ys, columns=['ground state'] + [str(i) for i in range(1,n_states)] + ['>{}'.format(n_states - 1)])
    df.plot(ax=ax, kind='bar', stacked=True, legend=False)

    ax.set_ylabel('${(c_l^i)}^2$')

    ax.legend(title="$l$-th Unperturbed state", frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlabel("$i$-th Perturbed state")
    ax.set_ylim(0., 1.)

    plt.tight_layout()

    ext = get_ext(save_fn)
    plt.savefig(save_fn[:-4] + '.' + ext, dpi=dpi, format=ext)
    logging.info('\n * Cumulative histogram'
                 ' saved as {}.\n'.format(save_fn[:-4] + '.' + ext))

    plt.show()

if __name__ == '__main__':
    pass
