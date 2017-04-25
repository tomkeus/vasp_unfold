#===========================================================
#
#  PROJECT: vasp_unfold
#  FILE:    unfolding.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.32
#  DATE:    Apr 25th 2017
#
#===========================================================

import numpy as np
from utils import post_error, frac_translation_order


def build_operators(spos, trans, check_mapping=False, eps=1e-6):
    '''Given a list of fractional positions, and a list of
    fractional translations, produce a set of matrices, 
    describing how fractional translations permute atomic
    positions within the unit cell. Two fractional positions
    si and sj are considered to be identical when |si-sj|<eps.
    '''
    ntrans = len(trans)
    natoms = len(spos)
    
    ops = np.zeros((ntrans, natoms, natoms), int)
    
    for i, ti in enumerate(trans):
        for j, sj in enumerate(spos):
            for k, sk in enumerate(spos):
                # If displacement between two atomic positions
                # differs from the fractional translations by
                # by a lattice translation+-eps we consider the
                # that two atomic positions to be map onto each
                # other by the fractional translation
                disp = sk-sj-ti
                ops[i, j, k] = np.linalg.norm(disp-np.rint(disp)) < eps

    if check_mapping:
        # Every row and every column of every operator must
        # contain exactly one unity, otherwise translations
        # are not maping atoms one-to-one.
        if np.any(np.sum(ops, axis=1) != 1) or np.any(np.sum(ops, axis=2) != 1):
            post_error('Translations are not one-to-one. '
                'Try changing the matching tolerance, or try using '
                'the POSCAR file with more regular positions.')

    return ops
    
    
def build_translations(tgens):
    '''Build a list of translations and irreps from at most 
    three linearly independent generators specified as lists
    of three fractions.
    '''
    if len(tgens) > 3:
        post_error('There can be at most three generators '
            'of fractional translations.')
    
    # Check if generators are linearly independent
    if len(tgens) == 2 and np.all(np.cross(tgens[0], tgens[1]) == 0):
        post_error('Generators are not linearly independant.')
    elif len(tgens) == 3 and np.linalg.det(tgens) == 0:
        post_error('Generators are not linearly independant.')
        
    # Expand the generator list to be a 3x3 matrix
    tgens = np.append(tgens, np.ones((3-len(tgens), 3)), axis=0)
    
    # Get the order of every generator
    order = np.array([frac_translation_order(g) for g in tgens], int)
    
    # Get the corresponding roots of unity which will be
    # used to construct irreps
    deltag = np.exp(-2*np.pi*1j/order)
    
    # Fold the translation vectors into the unit cell
    tgens = tgens % 1
    
    # Total number of translations
    ntrans = np.prod(order)
    
    # Storage for translations
    trans = np.zeros((ntrans, 3), float)
    
    # Calculate irreps of every generator
    irrepg = np.zeros((ntrans, 3), complex)
    
    # Loop over all possible products of powers of generators
    # and store the irreps
    for i in xrange(order[0]):
        for j in xrange(order[1]):
            for k in xrange(order[2]):
                ii = i*order[1]*order[2]+j*order[2]+k
                
                irrepg[ii] = deltag**[i, j, k]
    
    # Irrep table for all translations          
    irreps = np.zeros((ntrans, ntrans), complex)
    
    # Loop over all posible products of powers of generators
    for i in xrange(order[0]):
        for j in xrange(order[1]):
            for k in xrange(order[2]):
                ti = i*order[1]*order[2]+j*order[2]+k
                
                # Get the translation vector
                trans[ti] = np.dot([i, j, k], tgens)%1
                
                # Get all irreps of that translations
                for l in xrange(ntrans):
                    irreps[l, ti] = np.prod(irrepg[l]**[i, j, k])
                    
    return trans, irreps


def build_projectors(irreps, ops, intdim=1):
    '''Builds irrep projection operators, given irreps
    and corresponding translation operator matrices.
    Intdim specifies how many orbitals per atomic site 
    there are.
    '''
    projs = np.zeros((len(irreps), ops.shape[1], ops.shape[2]), complex)
    
    # Loop over irreps
    for i, ii in enumerate(irreps):
        # Sum operators multiplied by their irrep
        for j, oj in enumerate(ops):
            projs[i] += irreps[i, j]*oj
    
    # Expand onto the orbital space and normalize
    return np.kron(projs, np.eye(intdim))/len(ops)
