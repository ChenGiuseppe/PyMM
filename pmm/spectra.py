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

def calc_uv(energies: np.ndarray, pert_matrix: np.ndarray, sigma=0.0003) -> np.ndarray:
    '''Calculate UV spectrum from the PMM trajectory.

    Parameters:
        energies (np.ndarray): (n_frames, n_states) array of the perturbed
            energies trajectory.
        pert_matrixs (np.ndarray): (n_frames, n_states_n_states, 3) perturbed
            matrix trajectory.
        sigma (float): TODO: #5 add description to sigma
    Returns:
    '''
    fig, ax = plt.subplots(1, 1, figsize=(8,5))
    exc_ens = (energies[:, 1:] - np.expand_dims(energies[:, 0], axis=1)) / (2
              * np.pi)
    mu_squareds = (pert_matrix[:, 0, 1:, :]**2).sum(axis=2)
    n_frames = energies.shape[0]
    # histos = []
    extra_range = 0.03
    n_bins = int(round((exc_ens.max() + extra_range - (exc_ens.min() 
             - extra_range)) / (0.0016 / (2 * np.pi * 10))))
    bin_edges = np.histogram_bin_edges(exc_ens, range=(exc_ens.min() 
                                                       - extra_range,
                                                       exc_ens.max()
                                                       + extra_range),
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
    np.savetxt('tot_UV_spectrum.txt', tot_spectrum)
    np.savetxt('UV_transitions.txt', spectra)
    ax.plot(lambdas, tot_spectrum[:, 1] * 16863 / 2.3)
    for i in range(1, spectra.shape[1]):
        ax.plot(lambdas, spectra[:, i] * 16863 / 2.3)
    plt.show()

if __name__ == '__main__':
    def read_raw_matrix(matrix_fn: str) -> np.ndarray:
        '''Read a (N, N, 3) matrix directly from a text file
    (for each matrix element i j: "i j x y z\n").

    Parameters:
        matrix_fn (str): matrix filename.

    Returns:
        matrix (np.ndarray): matrix as a numpy array.
    '''
        with open(matrix_fn, 'r') as matrix_file:
            matrix_x = []
            matrix_y = []
            matrix_z = []
            for i, line in enumerate(matrix_file):
                content = line.split()
                # print(i, content)
                # Check if the first line is to be read.
                if i == 0 and not (len(content) == 3 or len(content) == 5):
                    continue
                # NOTE: skip empty line. It was added having in mind the
                # possibility that the last line could be empty, as it is now it
                # accounts for all empty lines.
                elif content == []:
                    continue
                elif not (len(content) == 3 or len(content) == 5):
                    raise ValueError("There's an error in the matrix file " +
                                     f'at line {i}: there are ' +
                                     f'{len(content)} elements.\n' +
                                     f'{line}')
                if len(content) == 3:
                    matrix_x.append(float(content[0]))
                    matrix_y.append(float(content[1]))
                    matrix_z.append(float(content[2]))
                elif len(content) == 5:
                    matrix_x.append(float(content[2]))
                    matrix_y.append(float(content[3]))
                    matrix_z.append(float(content[4]))
        n_rows = round(np.sqrt(len(matrix_x)))
        # Workaround to the maximum dimension of npdarray (32).
        matrix_tmp = np.zeros_like(matrix_x)
        # print(len(matrix_x), np.sqrt(len(matrix_x)), n_rows)
        matrix_tmp[:] = matrix_x
        matrix_x = np.copy(matrix_tmp.reshape((n_rows, n_rows)))
        matrix_tmp[:] = matrix_y
        matrix_y = np.copy(matrix_tmp.reshape((n_rows, n_rows)))
        matrix_tmp[:] = matrix_z
        matrix_z = np.copy(matrix_tmp.reshape((n_rows, n_rows)))
        # print(matrix_x, matrix_y, matrix_z)
        matrix = np.stack((matrix_x, matrix_y, matrix_z), axis=2)
        return matrix
    matrix = read_raw_matrix('/mnt/d/Dottorato/Programs/prova/thioindigo/matrix_pbe0')
    energies = np.loadtxt('/mnt/d/Dottorato/Programs/prova/thioindigo/eigenvals_last.txt')
    eigvecs = np.load('/mnt/d/Dottorato/Programs/prova/thioindigo/eigenvecs.npy')
    pert_matrix = calc_pert_matrix(matrix, eigvecs)
    calc_uv(energies, pert_matrix)
    #for i in range(100):
    #    print(i, energies[i,1]-energies[i,0], (pert_matrix[i,0,1,:]**2).sum())