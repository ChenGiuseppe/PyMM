from pymm.inputs import read_raw_matrix
from pymm.pmm import pmm
from pymm.spectra import calc_pert_matrix, calc_abs
import pymm.free_en import calc_dA
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
    parser_pmm.add_argument('-o', '--output', action='store', type=str,
                        default='eigvals.txt', help='perturbed QC energies '
                        '(default: eigenval.txt)')
    parser_pmm.add_argument('-oc', '--output_vecs', action='store', type=str,
                        default='eigvecs', help='perturbed eigenvectors '
                        '(default: eigvecs -> eigvecs.npy)')
    # Calculate absorption spectra.
    parser_abs = subparsers.add_parser('calc_abs',
                                       help='Calculate absorption spectra')
    parser_abs.add_argument('-dm', '--dip-matrix', action='store', type=str,
                            help='QC unperturbed electric dipole moment matrix '
                            'filename')
    parser_abs.add_argument('-el', '--eigvals', action='store', type=str,
                            default='eigvals.txt', help='Perturbed '
                            'eigenvalues trajectory (default: eigvals.npy)')
    parser_abs.add_argument('-ec', '--eigvecs', action='store', type=str,
                            default='eigvecs.npy', help='Perturbed '
                            'eigenvectors trajectory (default: eigvecs.npy)')
    parser_abs.add_argument('-sigma', action='store', type=float,
                            default=0.0003, help='Sigma value of the gaussian '
                            'broadening of the signal')
    parser_abs.add_argument('-ot', '--output', action='store', type=str,
                            default='abs_spectrum',
                            help='Calculated absorption spectra names'
                            '(default: abs_spectrum)')
    # Calculate free energy
    parser_dA = subparser.add_parser('calc_dA',
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
    parser_dA.add_argument('-eif', '--en_in_in', action='store', type=str,
                           help='Filename of the MD-PMM energies trajectory '
                           'for the initial state in the final ensemble')
    parser_dA.add_argument('-eff', '--en_fin_in', action='store', type=str,
                           help='Filename of the MD-PMM energies trajectory '
                           'for the final state in the final ensemble')   
    parser_dA.add_argument('-o', '--output', action='store', type=str,
                            default='dA_mean.txt',
                            help='Calculated free energy differences'
                            '(default: dA_mean.txt)')                    
    cmdline = parser.parse_args()

    start = timer()
    if cmdline.command == 'run_pmm':
        pmm(cmdline)
    elif cmdline.command == 'calc_abs':
        dip_matrix = read_raw_matrix(cmdline.dip_matrix)
        eigvals = np.loadtxt(cmdline.eigvals)
        eigvecs = np.load(cmdline.eigvecs)
        pert_matrix = calc_pert_matrix(dip_matrix, eigvecs)
        calc_abs(eigvals, pert_matrix, cmdline.output, cmdline.sigma)
    elif cmdline.command == 'calc_dA':
        dA = calc_dA(cmdline.T, cmdline.en_in_in, cmdline.en_fin_in,
                     cmdline.en_in_fin, cmdline.en_fin_fin)
        np.loadtxt(cmdline.output, dA)
    end = timer()
    print('The calculation took: ', end - start)
    return 0