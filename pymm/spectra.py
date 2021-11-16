import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def calc_pert_matrix(unpert_matrixs: np.ndarray,
                eigvecs: np.ndarray) -> np.ndarray:
    '''Calculate trajectory of the perturbed matrix of a generic property.

    Parameters:
        unpert_matrixs (np.ndarray): (n_states, n_states, 3) matrix of the
            unperturbed matrix.
        eigvecs (np.ndarray): (n_frames, n_states, n_states) matrix
            of the perturbed eigenvectors along the MD trajectory.

    Returns:
        pert_matrixs (np.ndarray): (n_frames, n_states_n_states, 3) perturbed
            matrix trajectory.
    '''

    pert_matrixs = np.einsum('ijk,mjn->mink', unpert_matrixs, eigvecs)
    pert_matrixs = np.einsum('ijk,ijmn->ikmn', eigvecs, pert_matrixs)
    return pert_matrixs

def calc_abs(energies: np.ndarray, pert_matrix: np.ndarray,
             output_fn: str, sigma: float, extra_range) -> np.ndarray:
    '''Calculate UV spectrum from the PMM trajectory.

    Parameters:
        energies (np.ndarray): (n_frames, n_states) array of the perturbed
            energies trajectory.
        pert_matrixs (np.ndarray): (n_frames, n_states_n_states, 3) perturbed
            matrix trajectory.
        sigma (float): TODO: #5 add description to sigma
    Returns:
    '''

    plt.style.use('seaborn-colorblind')

    fig, ax = plt.subplots(1, 1, figsize=(8,5))

    exc_ens = (energies[:, 1:] - np.expand_dims(energies[:, 0], axis=1)) / (2
              * np.pi)
    mu_squareds = (pert_matrix[:, 0, 1:, :]**2).sum(axis=2)
    n_frames = energies.shape[0]
    # histos = []

    vmin = exc_ens.min() - extra_range
    if vmin < 0.:
        vmin = 0.
    vmax = exc_ens.max() + extra_range
    n_bins = int(round((vmax - vmin) / (0.0016 / (2 * np.pi * 10))))
    bin_edges = np.histogram_bin_edges(exc_ens, range=(vmin,
                                                       vmax),
                                       bins=n_bins)
    bin_step = bin_edges[1] - bin_edges[0]
    bin_centers = bin_edges[:-1] + bin_step / 2
    spectra = []
    for i in range(exc_ens.shape[1]):
        histo, bin_edges_tmp = np.histogram(exc_ens[:, i], bins=bin_edges,
                                            weights=mu_squareds[:, i])
        # histos.append(histo)
        histo = histo * 4 * np.pi / (6 * n_frames)
        spectrum_tmp = np.zeros_like(bin_centers)
        for j, intensity in enumerate(histo):
            if intensity:
                spectrum_tmp += bin_centers * intensity *\
                                np.exp(-(bin_centers - bin_centers[j]) ** 2
                                       / (2 * sigma ** 2))
        spectra.append(spectrum_tmp * np.sqrt(2*np.pi) / (137 * sigma))

    spectra = np.array(spectra)
    tot_spectrum = np.sum(spectra, axis=0)
    lambdas = 137 * 0.0529 / bin_centers
    tot_spectrum = np.vstack((lambdas, tot_spectrum)).transpose((1,0))
    spectra = np.vstack((lambdas, spectra)).transpose((1,0))
    np.savetxt(f'{output_fn}_tot.dat', tot_spectrum)
    np.savetxt(f'{output_fn}_transitions.dat', spectra)

    logging.info(' * The total UV-Vis spectrum has been saved as {}'.format(f'{output_fn}_tot.dat'))
    logging.info(' * All the transitions considered on their own have been saved'
                 ' as {}'.format(f'{output_fn}_transitions.dat'))

    for i in range(1, spectra.shape[1]):
        ax.plot(lambdas, spectra[:, i] * 16863 / 2.3, label=r'0 $\longrightarrow$ {}'.format(i))
    ax.plot(lambdas, tot_spectrum[:, 1] * 16863 / 2.3, color='k',
            linestyle='--', label='Total spectrum')

    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('$\epsilon$ (M$^{-1}$cm$^{-1}$)')

    ax.legend(frameon=False)

    plt.show()

    return 0

if __name__ == '__main__':
    pass