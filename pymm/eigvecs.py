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

def eig_corr(eig, l, m, state, save_fn, dpi):
    '''
    '''

    plt.style.use('seaborn-colorblind')
    fig, ax = plt.subplots(1, 1)

    df = pd.DataFrame(eig[:,:,state], columns = [str(i) for i in range(eig.shape[1])])

    sns.histplot(df, x=str(l), y=str(m), bins=30, cmap='viridis')

    ext = get_ext(save_fn)
    plt.savefig(save_fn[:-4] + '.' + ext, dpi=dpi, format=ext)

    plt.show()

def eig_corr_tot(eig, state, save_fn, dpi):
    '''
    '''

    plt.style.use('seaborn-colorblind')

    df = pd.DataFrame(eig[:,:,state], columns = [str(i) for i in range(eig.shape[1])])
    g = sns.PairGrid(df)
    g.map_diag(sns.kdeplot,
               color='k'
               )
    g.map_offdiag(sns.histplot,
                  bins=30,
                  #sns.kdeplot,
                  #fill=True,
                  cmap='viridis',
                  #sns.scatterplot,
                  #s=0.5
                  )

    ext = get_ext(save_fn)
    plt.savefig(save_fn[:-4] + '.' + ext, dpi=dpi, format=ext)
    
    plt.show()



def eig_hist(eig, n_states, save_fn, dpi):
    '''
    '''

    plt.style.use('seaborn-colorblind')
    fig, ax = plt.subplots(1, 1)

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

    ax.set_ylabel('$(c_l^i)^2$')

    ax.legend(title="$l$-th Unperturbed state", frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlabel("$i$-th Perturbed state")
    ax.set_ylim(0., 1.)

    plt.tight_layout()

    ext = get_ext(save_fn)
    plt.savefig(save_fn[:-4] + '.' + ext, dpi=dpi, format=ext)

    plt.show()

if __name__ == '__main__':
    pass