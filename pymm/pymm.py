from pymm.conversions import au2eV
from pymm.eigvecs import get_ext, eig_corr, eig_corr_tot, eig_hist
from pymm.inputs import read_raw_matrix
from pymm.pmm import pmm
from pymm.spectra import calc_pert_matrix, calc_abs
from pymm.free_en import calc_dA, calc_dA_mean
from datetime import datetime
import getpass
import logging
import sys
from timeit import default_timer as timer
from argparse import ArgumentParser, FileType
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def main():
    '''Program with CLI to perform MD-PMM calculations.
    TODO: #3 add documentation.
    '''
    parser = ArgumentParser(prog='PyMM',
                            description='Perform MD-PMM calculation.')
    subparsers = parser.add_subparsers(dest='command')

    # Perform MD-PMM calculation.
    parser_pmm = subparsers.add_parser('run_pmm',
                                       help='Perform MD-PMM calculation.')
    parser_pmm.add_argument('-g', '--ref-geom', action='store', type=str,
                        help='QC reference (QM) geometry filename')
    parser_pmm.add_argument('-gu', '--geom-units', choices=['angstrom', 'bohr',
                        'nm'], help='specify units used in the reference (default: angstrom)'
                        'geometry', default='angstrom')
    parser_pmm.add_argument('-dm', '--dip-matrix', action='store', type=str,
                        help='QC unperturbed electric dipole moment matrix '
                        'filename')
    parser_pmm.add_argument('-e', '--energies', action='store', type=str,
                        help='QC unperturbed electronic states energies '
                        'filename')
    parser_pmm.add_argument('-ch', '--charges', action='store', type=str,
                        default=False, help='file with the QC atomic '
                        'charges for each unperturbed electronic state '
                        '(default=False)')
    parser_pmm.add_argument('-traj', '--trajectory-path', action='store',
                        type=str, help='XTC file of the MD simulation '
                        'trajectory')
    parser_pmm.add_argument('-top', '--topology-path', action='store', type=str,
                        help='TPR file of the MD simulation')
    parser_pmm.add_argument('-q', '--qc-charge', action='store', default=0,
                        type=int, help='total QC charge (default=0)')
    parser_pmm.add_argument('-nm', '--mm-indexes', action='store', type=str,
                        help='indexes of the QC in the MM trajectory')
    parser_pmm.add_argument('-nq', '--qm-indexes', action='store', type=str,
                        default=False, help='indexes of the portion of the '
                        'reference (QM) geometry to be considered in the '
                        'MD-PMM calculation')
    parser_pmm.add_argument('--match', action='store_true',
                            help='Reorder the QC reference geometry '
                            'to match the atoms order in the MD simulation')
    parser_pmm.add_argument('-o', '--output', action='store', type=str,
                        default='pymm', help='job name used to name the '
                        'output files (default: pymm)')

    # Calculate absorption spectra.
    parser_abs = subparsers.add_parser('calc_abs',
                                       help='Calculate absorption spectra')
    parser_abs.add_argument('-dm', '--dip-matrix', action='store', type=str,
                            help='QC unperturbed electric dipole moment matrix '
                            'filename')
    parser_abs.add_argument('-el', '--eigvals', action='store', type=str,
                            default='pymm_eigvals.dat', help='Perturbed '
                            'eigenvalues trajectory (default: pymm_eigvals.dat)')
    parser_abs.add_argument('-ev', '--eigvecs', action='store', type=str,
                            default='pymm_eigvecs.npy', help='Perturbed '
                            'eigenvectors trajectory (default: pymm_eigvecs.npy)')
    parser_abs.add_argument('-sigma', action='store', type=float,
                            default=0.034, help='Square root of the variance of the gaussian '
                            'broadening applied to the signal (expressed in eV). Default: 0.034 eV.')
    parser_abs.add_argument('-xtr', '--extra-range', action='store', type=float,
                            default=0.005, help='Additional value to add to the'
                            'frequency range considered for the spectrum '
                            'calculation')
    parser_abs.add_argument('-ot', '--output', action='store', type=str,
                            default='abs_spectrum',
                            help='Calculated absorption spectra names '
                            '(default: abs_spectrum)')

    # Calculate free energy
    parser_dA = subparsers.add_parser('free_en',
                                     help='Calculate free energy')
    parser_dA.add_argument('-T', '--temperature', action='store', type=float,
                           default=298., help='Temperature of the system '
                           'during the simulation')
    parser_dA.add_argument('-eii', '--en_in_in', action='store', type=str,
                           help='Filename of the MD-PMM energies trajectory '
                           'for the initial state in the initial ensemble')
    parser_dA.add_argument('-efi', '--en_fin_in', action='store', type=str,
                           help='Filename of the MD-PMM energies trajectory '
                           'for the final state in the initial ensemble')
    parser_dA.add_argument('-eif', '--en_in_fin', action='store', type=str,
                           help='Filename of the MD-PMM energies trajectory '
                           'for the initial state in the final ensemble')
    parser_dA.add_argument('-eff', '--en_fin_fin', action='store', type=str,
                           help='Filename of the MD-PMM energies trajectory '
                           'for the final state in the final ensemble')
    parser_dA.add_argument('-col', action='store', type=int, default=0,
                           help='Column (or electronic state) to consider '
                           'from the energies trajectory file') 
    parser_dA.add_argument('-o', '--output', action='store', type=str,
                            default='dA_mean.dat',
                            help='Calculated free energy differences'
                            '(default: dA_mean.dat)')                    

    # Eigenvectors analysis 
    parser_eig = subparsers.add_parser('eig', help='Analyse the perturbed '
                                       'eigenvectors')
    parser_eig.add_argument('-i', '--input', action='store', type=str,
                            default='pymm_eigvecs.npy',
                            help='Perturbed eigenvectors trajectory provided '
                            'as an npy file (n_frames, n_unp_states, n_per_states;'
                            'default=pymm_eigvecs.npy).')
    parser_eig.add_argument('-first', action='store', type=int, default=0,
                            help='First unperturbed state to consider (ignored'
                            ' when calculating the cumulative histograms; '
                            'default=0, i.e. the electronic ground state).')
    parser_eig.add_argument('-last', action='store', type=int,
                            help='Last unperturbed state to consider.')
    parser_eig.add_argument('-state', action='store', type=int,
                            help='Perturbed state to consider.')
    parser_eig.add_argument('-dpi', action='store', type=int, default=300,
                            help='Dpi of the saved image (default=300).')
    parser_eig.add_argument('-bins', action='store', type=int, default=20,
                            help='Number of bins used to calculate '
                            'the histograms (default=20).')
    parser_eig.add_argument('-oc', '--output_corr', action='store', type=str,
                            const='eig_corr.png', nargs='?',
                            help='Calculate the correlation between a pair of'
                            ' coefficients of the l-th and m-th unperturbed '
                            'states (provided by -first and -last) of the i-th '
                            'perturbed state (-state). Choose the preferred file '
                            'extention (default=eig_corr.png).')
    parser_eig.add_argument('-ot', '--output_tot', action='store', type=str,
                            const='eig_corr_tot.png', nargs='?',
                            help='Calculate the correlation between all '
                            'the pairs of coefficients of the l-th and m-th '
                            'unperturbed states of the i-th perturbed state '
                            '(provided by -state). Choose the preferred file '
                            'extention (default=eig_corr_tot.png).')
    parser_eig.add_argument('-oh', '--output_hist', action='store', type=str,
                            const='eig_hist.png', nargs='?',
                            help='Calculate the mean squared coefficients with which each '
                            'unperturbed state contributes to each perturbed '
                            'state. Choose the preferred file extention '
                            '(default=eig_hist.png).')

    cmdline = parser.parse_args()

    logging.basicConfig(format='%(message)s',
                        level=logging.INFO)
    logging.info('==========================================================================\n'
                 '|    PyMM: A computational package for PMM-MD simulations in Python.     |\n'
                 '==========================================================================\n\n'
                 'User: {}\n'.format(getpass.getuser()) +
                 'Date: {}'.format(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')))
    logging.info('\nJob launched:')
    logging.info('{}\n\n\n'.format(' '.join(sys.argv)))

    if cmdline.command == 'run_pmm':
        start = timer()

        logging.info('==========================================================================\n'
                     '|                                                                        |\n'
                     '|                  Launching PMM-MD calculation:                         |\n'
                     '|                                                                        |\n'
                     '==========================================================================\n')

        pmm(cmdline)
    
        end = timer()

        logging.info('\n'
                     '==========================================================================\n'
                     '|                                                                        |\n'
                     '|                                 FINISHED                               |\n'
                     '|                                                                        |\n'
                     '==========================================================================\n'
                     '\n'
                     'Date: {}.'.format(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')) + 
                     '\n'
                     'The calculation took: {}'.format(end - start) +
                     '\n\n'
                     '--------------------------------------------------------------------------\n')

    elif cmdline.command == 'calc_abs':

        logging.info('==========================================================================\n'
                     '|                                                                        |\n'
                     '|                     Calculate UV-Vis spectrum:                         |\n'
                     '|                                                                        |\n'
                     '==========================================================================\n')

        # Change matplotlib settings:
        plt.rcParams['axes.linewidth'] = 1.3
        plt.rcParams['font.size'] = 18
        plt.rcParams['legend.fontsize'] = 15

        dip_matrix = read_raw_matrix(cmdline.dip_matrix)
        eigvals = np.loadtxt(cmdline.eigvals)
        eigvecs = np.load(cmdline.eigvecs)
        pert_matrix = calc_pert_matrix(dip_matrix, eigvecs)

        logging.info(' * Number of frames = {}'.format(eigvals.shape[0]))
        logging.info(' * Number of transitions = {}'.format(eigvals.shape[1] - 1))
        logging.info(f' * sigma = {cmdline.sigma} eV.')

        calc_abs(eigvals, pert_matrix, cmdline.output, cmdline.sigma / (2*np.pi*au2eV),
                 cmdline.extra_range)
    
        logging.info('\n'
                 '==========================================================================\n'
                 '|                                                                        |\n'
                 '|                                 FINISHED                               |\n'
                 '|                                                                        |\n'
                 '==========================================================================\n'
                 '\n'
                 'Date: {}.'.format(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')) + 
                 '\n\n'
                 '--------------------------------------------------------------------------\n')

    elif cmdline.command == 'free_en':

        logging.info('==========================================================================\n'
                     '|                                                                        |\n'
                     '|                     Calculate Free Energy Change:                      |\n'
                     '|                                                                        |\n'
                     '==========================================================================\n')

        col = 0

        if cmdline.en_in_in and cmdline.en_fin_in and cmdline.en_in_fin and cmdline.en_fin_fin:

            logging.info(' * Free energy will be calculated considering two '
                         'ensembles.')
            en_in_in = np.loadtxt(cmdline.en_in_in)[:,col]
            en_fin_in = np.loadtxt(cmdline.en_fin_in)[:,col]
            en_in_fin = np.loadtxt(cmdline.en_in_fin)[:,col]
            en_fin_fin = np.loadtxt(cmdline.en_fin_fin)[:,col]

            if en_in_in.shape[0] == en_in_fin.shape[0] == en_fin_in.shape[0] == en_fin_fin.shape[0]:

                logging.info(' * Number of frames = {}\n'.format(en_in_in.shape[0]))
                dA = calc_dA_mean(cmdline.temperature, en_in_in, en_fin_in,
                                  en_in_fin, en_fin_fin)

            else:
                logging.error(' ! The trajectories are not of the same length\n'
                              f'   * eii: {en_in_in.shape[0]} \n'
                              f'   * efi: {en_fin_in.shape[0]} \n'
                              f'   * eif: {en_in_fin.shape[0]} \n'
                              f'   * eff: {en_fin_fin.shape[0]} \n\n')
                raise IOError

        elif (cmdline.en_in_in and cmdline.en_fin_in) and not cmdline.en_in_fin and not cmdline.en_fin_fin:

            logging.info(' * Free energy will be calculated considering only the '
                         'ensemble in the initial state.')
            logging.warning(' ! For a more rigorous estimation of the free energy '
                            'consider using both the initial and final ensembles '
                            '(see the documentation).')
            en_in_in = np.loadtxt(cmdline.en_in_in)[:,col]
            en_fin_in = np.loadtxt(cmdline.en_fin_in)[:,col]

            if en_in_in.shape[0] == en_fin_in.shape[0]:

                logging.info(' * Number of frames = {}\n'.format(en_in_in.shape[0]))
                dA = calc_dA(cmdline.temperature, en_in_in, en_fin_in, state='ini')

            else:
                logging.error(' ! The trajectories are not of the same length\n'
                              f'   * eii: {en_in_in.shape[0]} \n'
                              f'   * efi: {en_fin_in.shape[0]} \n\n')
                raise IOError

        elif not cmdline.en_in_in and not cmdline.en_fin_in and (cmdline.en_in_fin and cmdline.en_fin_fin):

            logging.info(' * Free energy will be calculated considering only the '
                         'ensemble in the final state.')
            logging.warning(' ! For a more rigorous estimation of the free energy '
                            'consider using both the initial and final ensembles '
                            '(see the documentation).')

            en_in_fin = np.loadtxt(cmdline.en_in_fin)[:,col]
            en_fin_fin = np.loadtxt(cmdline.en_fin_fin)[:,col]

            if en_in_fin.shape[0] == en_fin_fin.shape[0]:

                logging.info(' * Number of frames = {}\n'.format(en_in_fin.shape[0]))
                dA = calc_dA(cmdline.temperature, en_in_fin, en_fin_fin, state='ini')
    
            else:
                logging.error(' ! The trajectories are not of the same length\n'
                              f'   * eii: {en_in_fin.shape[0]} \n'
                              f'   * efi: {en_fin_fin.shape[0]} \n\n')
                raise IOError

        else:
            raise IOError('Input files provided incorrectly. See documentation.')

        logging.info('================================= RESULTS ================================\n'
                     '                     * Calculated Free Energy Change:\n'
                     '                       {} J/mol\n'.format(dA) +
                     '==========================================================================\n')

        logging.info('\n'
                     '==========================================================================\n'
                     '|                                                                        |\n'
                     '|                                 FINISHED                               |\n'
                     '|                                                                        |\n'
                     '==========================================================================\n'
                     '\n'
                     'Date: {}.'.format(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')) + 
                     '\n\n'
                     '--------------------------------------------------------------------------\n')

    elif cmdline.command == 'eig':

        logging.info('==========================================================================\n'
                     '|                                                                        |\n'
                     '|                     Eigenvectors Analysis:                             |\n'
                     '|                                                                        |\n'
                     '==========================================================================\n')

        eigvecs = np.load(cmdline.input)

        logging.info(' * Number of frames = {}'.format(eigvecs.shape[0]))
        logging.info(' * Number of electronic states = {}\n\n'.format(eigvecs.shape[1]))

        # Change matplotlib settings:
        plt.rcParams['axes.linewidth'] = 1.3
        plt.rcParams['font.size'] = 18
        plt.rcParams['legend.fontsize'] = 15

        if cmdline.output_corr is not None:

            logging.info(' * Correlation plot between the contributions arising from the unperturbed\n'
                         f'   states {cmdline.first} and {cmdline.last} to the perturbed state {cmdline.state}.\n'
                         )

            eig_corr(eigvecs, cmdline.first, cmdline.last, cmdline.state,
                     cmdline.output_corr, cmdline.dpi, cmdline.bins)
        if cmdline.output_tot is not None:

            logging.info(' * Correlation plot between all the pairs obtained between the contributions\n'
                         f'   arising from the unperturbed states to the perturbed state {cmdline.state}.\n'
                         )

            eig_corr_tot(eigvecs, cmdline.state, cmdline.output_tot,
                         cmdline.dpi, cmdline.bins)
        if cmdline.output_hist is not None:

            logging.info(' * Cumulative histogram of the mean contribution arising from each \n'
                         f'   unperturbed state to each the perturbed state.\n'
                         )

            logging.info(f' * Only the first {cmdline.state} states (including the ground state) have\n'
                         '   been explicitly considered.\n')

            eig_hist(eigvecs, cmdline.state, cmdline.output_hist,
                     cmdline.dpi)

        logging.info('\n'
                     '==========================================================================\n'
                     '|                                                                        |\n'
                     '|                                 FINISHED                               |\n'
                     '|                                                                        |\n'
                     '==========================================================================\n'
                     '\n'
                     'Date: {}.'.format(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')) + 
                     '\n\n'
                     '--------------------------------------------------------------------------\n')
    
    return 0