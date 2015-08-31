#! /usr/bin/env python


#============================================================================
#
#  PROJECT: vasp_unfold
#  FILE:    __main__.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.3
#  DATE:    Aug 31st 2015
#
#
# vasp_unfold is a code whose purpose is to perform the unfolding
# of band structures calculated with VASP. It is based on the method
# described in (Phys. Rev. B 90, 195121). If this code have been used to 
# obtain the data in your publication please cite the given reference.
#
# This code is provided as is. There is no guarantee that it will work.
# However, if you encounter errors or unexpected behavor, please report
# it to tomic@itp.uni-frankfurt.de and I will try to improve it.
#
#============================================================================


import numpy as np
import argparse
import sys
from utils import post_error, translation
from unfolding import build_translations, build_operators, build_projectors
from parse import parse_poscar, parse_procar
from write import write_procar
import errors

def main():
    desc_str = 'Unfold bands calculated by VASP. For this, phase '\
               'information needs to be present in the PROCAR file '\
               'which means that bands need to be calculated with '\
               'the option LORBIT=12 specified in the INCAR file. '\
               'The required input consists of POSCAR file, PROCAR '\
               'file and a set of up to three fractional translations. '\
               'There is no need to specify all translations, just the '\
               'generators. The programs will generate all distinct '\
               'translations from the generators and generate all irreps. '\
               'NOTE: In cas of non-collinear spin-polarized calculations ' \
               'mx, my and mz components of orbital weights are not ' \
               'unfolded, only spin up and spin down totals.'
               
               
    parser = argparse.ArgumentParser(prog='vasp_unfold', description = desc_str)

    parser.add_argument('poscar', type=str, help='POSCAR file')     
    parser.add_argument('procar', type=str, help='PROCAR file')

    parser.add_argument('--tgen', type=translation, action='append',
                        metavar='SX,SY,SZ', help='Fractional translation '
                        'generator. No whitespaces are allowed between the '
                        'components! SX, SY and SZ can be either 0 or 1/n, '
                        'where n is an integer. Up to three linearly '
                        'independant generators can be specified.')

    parser.add_argument('--out', type=str, help='Output filename. If left '
                        'unspecified  output is writen to PROCAR.irrep.n '
                        'where PROCAR is location of the input PROCAR file. '
                        'If specified, output is written to OUT.irrep.n.')
    
    parser.add_argument('--eps', type=float, default=1e-6, help='Numerical '
                        'precision. When building permutation representation '
                        'of the fractional translations this parameter is used '
                        'to determine whether two fractional positions are '
                        'identical. For irregular structures, it may need to '
                        'be increased. If this does not help, try to tweak '
                        'the atomic positions in POSCAR to make the structure '
                        'more regular. Default is 1e-6.')
    
    parser.add_argument('--all-irreps', default=False, action='store_true',
                        help='Flag specifying that all irreps from the unfolding '
                        'will be written to the output. By default, only irrep 0 '
                        '(unit irrep) is written out.')
                        
    parser.add_argument('--check-mapping', action='store_true', default=False,
                        help='Specifies to check if supplied translation '
                        'generators generate fractional translations which are '
                        'one-to-one, ie. map every atom on exactly one other atom '
                        'in the unit cell. This MUST not be enabled for the cases '
                        'where vacancies or excess atoms are present.')
                                                                     
    args = parser.parse_args()
    
    tgens = args.tgen
    
    trans, irreps = build_translations(tgens)
    
    try:
        cell, spos, symbols = parse_poscar(args.poscar)
    except:
        post_error('Unable to parse the input POSCAR file. Please '
            'check if the file exists and is formatted properly')
    
    ops = build_operators(spos, trans, args.check_mapping, args.eps)
    
    try:
        data = parse_procar(args.procar)
    except:
        post_error(errors.poscar_parse_error)
    
    if data[-1] is None:
        post_error('Phase information has to be present in the PROCAR '
            'file. Please repeat the calculation with LORBIT=12.')

    phases = np.copy(data[-1])
    
    norbs = phases.shape[1]/len(spos)
    
    projs = build_projectors(irreps, ops, norbs)
    
    if args.out is None:
        output = args.procar
    else:
        output = args.out
    
    if args.all_irreps:
        nirrep = len(projs)
    else:
        nirrep = 1
        
    for i, p in enumerate(projs[:nirrep]):
        for s in xrange(data[-1].shape[-1]):
            try:
                data[-1][:,:,:,s] = np.dot(p, phases[:,:,:,s]).swapaxes(0, 1)
            except:
                post_error('Unable to apply projectors. Are you sure '
                    'that specified POSCAR and PROCAR file belong to '
                    'the same crystal structure?')
        
        # We update total absolute weights to correspond to
        # unfolded phases (by multiplying them by the magnitude
        # ratio of unfolded and folded phases 
        phase_ratio = np.abs(data[-1])/(np.abs(phases)+1e-4)
        
        for idim in xrange(data[-2].shape[3]):
            data[-2][:,:,:,idim,:] *= phase_ratio

        write_procar('{0}.irrep.{1}'.format(output, i), *data)
       
       
if __name__ == '__main__':
    main()
    

