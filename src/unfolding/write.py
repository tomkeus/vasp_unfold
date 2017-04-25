#===========================================================
#
#  PROJECT: vasp_unfold
#  FILE:    write.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.32
#  DATE:    Apr 25th 2017
#
#===========================================================

import numpy as np
from utils import post_error


def write_procar(fname, orbitals, kpoints, kweights, bands, 
                 occupations, weights, phases):
    '''Write PROCAR file based on supplied data
    '''
    # Labels for orbitals
    orblabels = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                 'f-3', 'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3']
    
    norb = len(orbitals)
    npoints = len(kpoints)
    nbands = bands.shape[1]
    nions = weights.shape[1]/norb
    nspin = weights.shape[-1]
    ndim = weights.shape[-2]
    
    try:
        out = open(fname, 'w')
    except:
        post_error('Unable to open "{0}" for writing'.format(fname))
    
    # Write the first line of the PROCAR file
    if phases is not None:
        out.write('PROCAR lm decomposed + phase\n')
    else:
        out.write('PROCAR lm decomposed\n')

    # Format out the column title line for orbital weights
    orb_ttl_1 = 'ion '
    orb_ttl_2 = 'ion '
    
    for i in xrange(norb):
        orb_ttl_1 += '{0: >6} '.format(orblabels[i])
        orb_ttl_2 += '{0: >6} '.format(orblabels[i])
    
    orb_ttl_1 += '{0: >6}\n'.format('tot')
    orb_ttl_2 += '\n'
    
    for s in xrange(nspin):
        # Write the second line containing the sizes
        out.write('# of k-points:  {0}         # of bands:  {1}'
              '         # of ions:   {2}\n\n'.format(npoints, nbands, nions))
              
        for i, k in enumerate(kpoints):
            # Write the k-point info
            out.write(' k-point {0: >4} :    '.format(i+1))
            out.write('{0:.8f} {1:.8f} {2:.8f}     '.format(*k))
            out.write('weight = {0:.8f}\n\n'.format(kweights[i]))
            
            for j, b in enumerate(bands[i,:,s]):
                # Write the band info
                out.write('band {0: >4} # '.format(j+1))
                out.write('energy {0: >13.8f} # '.format(b))
                out.write('occ. {0: >11.8f}\n\n'.format(occupations[i, j, s]))
                
                # Write absolute weight blocks. In case of 
                # non-collinear calculation, there is four
                # such blocks
                for d in xrange(ndim):
                    if d == 0:
                        # Write names of orbitals (s, px, py etc...)
                        out.write(orb_ttl_1)
                    
                    # Allocate storage for accumulation of orbital totals
                    tot_orb = np.zeros(norb, float)
                    
                    # Loop over individual atoms
                    for k in xrange(nions):
                        # Write atom's index
                        out.write('{0: >3} '.format(k+1))
                        
                        # Extract corresponding weights
                        w = weights[i, k*norb:(k+1)*norb, j,d,s]
                        
                        # Add to totals
                        tot_orb += w
                        
                        # Write out row of weights
                        for l in xrange(norb):
                            out.write('{0: >6.3f} '.format(w[l]))
                        
                        # Write atom total
                        out.write('{0: >6.3f}\n'.format(np.sum(w)))
                    
                    # We will now write line with orbital totals
                    out.write('tot ')
                    
                    for l in xrange(norb):
                        out.write('{0: >6.3f} '.format(tot_orb[l]))
                    
                    # Finally, atom+orbital total
                    out.write('{0: >6.3f}\n'.format(np.sum(tot_orb)))
                
                # If we have phases we write them now
                if phases is not None:
                    # Write again names of orbitals
                    out.write(orb_ttl_2)
                    
                    # Loop over atoms
                    for k in xrange(nions):
                        # Write atom index
                        out.write('{0: >3} '.format(k+1))
                        
                        # Extract corresponding phase
                        phs = phases[i, k*norb:(k+1)*norb, j, s]
                        
                        # Write real part
                        for l in xrange(norb):
                            out.write('{0: >6.3f} '.format(phs[l].real))
                        
                        # Write atom index
                        out.write('\n{0: >3} '.format(k+1))
                        
                        # Write imaginary part
                        for l in xrange(norb):
                            out.write('{0: >6.3f} '.format(phs[l].imag))
                        
                        out.write('\n')
                    
                    out.write('\n')
                
            out.write('\n')
                
    out.close()

