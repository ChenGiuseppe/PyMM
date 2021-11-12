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
    # Get the QC QM properties from text files.
    parser_pmm.add_argument('-g', '--ref-geom', action='store', type=str,
                        help='QC reference (QM) geometry filename')
    parser_pmm.add_argument('-gu', '--geom-units', choices=['angstrom', 'bohr',
                        'nm'], help='specify units used in the reference '
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
    parser_pmm.add_argument('-m', '--match', action='store_true',
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
                            default='eigvals.dat', help='Perturbed '
                            'eigenvalues trajectory (default: eigvals.dat)')
    parser_abs.add_argument('-ec', '--eigvecs', action='store', type=str,
                            default='eigvecs.npy', help='Perturbed '
                            'eigenvectors trajectory (default: eigvecs.npy)')
    parser_abs.add_argument('-sigma', action='store', type=float,
                            default=0.0003, help='Sigma value of the gaussian '
                            'broadening of the signal (expressed as frequency)')
    parser_abs.add_argument('-xtr', '--extra-range', action='store', type=float,
                            default=0.005, help='Additional value to add to the'
                            'frequency range considered for the spectrum '
                            'calculation')
    parser_abs.add_argument('-ot', '--output', action='store', type=str,
                            default='abs_spectrum',
                            help='Calculated absorption spectra names'
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
    cmdline = parser.parse_args()

    if cmdline.command == 'run_pmm':
        start = timer()
        logging.basicConfig(#filename='output.log', filemode='w',
                            format='%(message)s',
                            level=logging.INFO)
        logging.info('PyMM: A computational package for PMM-MD simulations in Python.\n\n' +
                     'User: {}\n'.format(getpass.getuser()) +
                     'Date: {}'.format(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')))
        logging.info('Command used:')
        logging.info('{}\n'.format(' '.join(sys.argv)))
        logging.info('=========================================================')
        logging.info('Launching PMM-MD calculation:\n')
        logging.info('=========================================================')
        pmm(cmdline)
        end = timer()
        logging.info('\nFINISHED\nDate: {}.'.format(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')) + 
                     '\nThe calculation took: {}'.format(end - start))
    elif cmdline.command == 'calc_abs':
        dip_matrix = read_raw_matrix(cmdline.dip_matrix)
        eigvals = np.loadtxt(cmdline.eigvals)
        eigvecs = np.load(cmdline.eigvecs)
        pert_matrix = calc_pert_matrix(dip_matrix, eigvecs)
        calc_abs(eigvals, pert_matrix, cmdline.output, cmdline.sigma, cmdline.extra_range)
    elif cmdline.command == 'free_en':
        col = 0
        en_in_in = np.loadtxt(cmdline.en_in_in)[:,col]
        en_fin_in = np.loadtxt(cmdline.en_fin_in)[:,col]
        en_in_fin = np.loadtxt(cmdline.en_in_fin)[:,col]
        en_fin_fin = np.loadtxt(cmdline.en_fin_fin)[:,col]
        dA = calc_dA_mean(cmdline.temperature, en_in_in, en_fin_in,
                     en_in_fin, en_fin_fin)
        with open(cmdline.output, 'w') as file_out:
            file_out.write('Free Energy\n')
            file_out.write(f'{dA} J/mol')
    
    return 0