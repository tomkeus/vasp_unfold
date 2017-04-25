#===========================================================
#
#  PROJECT: vasp_unfold
#  FILE:    parse.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.32
#  DATE:    Apr 25th 2017
#
#===========================================================

import numpy as np
from utils import Getlines, post_error

def parse_poscar(filename):
    '''Parses POSCAR file. Returns 3x3 float array
    containing unit cell vectors as rows, Nx3
    array containing fractional positions of atoms,
    N being number of atoms and a list of chemical
    symbols.
    '''
    try:
        gl = Getlines(filename)
    except:
        post_error('Unable to open "{0}" for reading.'.format(filename))
    
    # Skip the comment line
    gl.readline()
    
    # Read the scaling factor
    f = float(gl.readline())
    
    # Read the unit cell vectors
    cell = np.zeros((3, 3), float)
    
    cell[0] = f*np.array(gl.readline().split(), float)
    cell[1] = f*np.array(gl.readline().split(), float)
    cell[2] = f*np.array(gl.readline().split(), float)
    
    # Read the chemical symbols
    syms = gl.readline().split()
    
    # Read the atom counts
    counts = np.array(gl.readline().split(), int)
    
    # Get the number of atoms
    natoms = np.sum(counts)
    
    # Rearrange into the long list of chemical symbols
    symbols = []
        
    for s, c in zip(syms, counts):
        symbols += [s]*c
    
    # Cartesian or fractional coordinates?
    ctype = gl.readline()[0].lower()
    
    if ctype == 'c':
        mult = np.linalg.inv(cell)
    elif ctype == 'd':
        mult = np.eye(3)
    else:
        post_error('"{0}" is unknown POSCAR option'.format(plines[7].strip()))
    
    # Allocate storage for positions
    spos = np.zeros((len(symbols), 3))
    
    # Read the positions
    for i in xrange(len(symbols)):
        spos[i] = np.array(gl.readline().split()[:3], float)
    
    # If necessary, this will convert from
    # Cartesian to fractional coordinates
    spos = np.dot(spos, mult)
    
    return cell, spos, symbols
    
    
def parse_procar(filename):
    '''This function parses a PROCAR file. It returns a tuple
    consisting of following elements:
    
    orbitals    - (norbs) string array of orbital labels (s, px, py etc...)
    kpoints     - (npoints,3) float array of k-point coordinates
    kweights    - (npoints) float array of k-point weights
    bands       - (npoints,nbands,nspin) float array of band energies
    occupancies - (npoints,nbands,nspin) float array of band occupancies
    weights     - (npoints,nions*norbs,nbands,ndim,nspin) float array
                  of orbital weights
    phases      - (npoints,nions*norbs,nbands,nspin) complex array of
                  phases of orbital weights if LORBIT=12, otherwise None
                  
    Where:
    
    norbs   - number of orbitals (It can be 9 or 16 with f orbitals)
    npoints - number of k-points
    nbands  - number of bands
    nspin   - number of spins (1 for non spin-polarized, 2 otherwise)
    nions   - number of atoms
    ndim    - orbital weight dimensionality (1 for collinear, 4 otherwise)
    '''

    try:
        gl = Getlines(filename)
    except:
        post_error('Unable to open "{0}" for reading.'.format(filename))

    header_1 = gl.readline()
    header_2 = gl.readline().split()
    
    npoints = int(header_2[3])
    nbands = int(header_2[7])
    nions = int(header_2[-1])
    
    # Remember the position in file
    start = gl.tell()
    
    # Skip two lines containing first k-point and band
    gl.readline()
    gl.readline()
    
    # Determine the number of orbitals
    orbitals = gl.readline().split()[1:-1]
    
    norbs = len(orbitals)
    
    # Allocate maximal storage. In case calculation was
    # non spin-polarized or collinear we can just trim
    # the excess components at the end
    kpoints = np.zeros((npoints, 3), float)
    kweights = np.zeros(npoints, float)
    bands = np.zeros((npoints, nbands, 2), float)
    occupancies = np.zeros((npoints, nbands, 2), float)
    weights = np.zeros((npoints, nions*norbs, nbands, 4, 2), float)
    
    dim = 0
    
    # Determine if the calculation was non-collinear
    # by counting how many lines in the first band
    # block begin with tot. That number will be equal
    # to the number of sub-blocks for orbital weights
    # (1 in case of collinear and 4 otherwise)
    while True:
        line = gl.readline()
        
        if line.startswith('tot'):
            dim += 1
        elif line.startswith('band'):
            break
    
    # This function will read block of absolute weights
    # for i-th k-point, j-th band and s-th spin component
    def get_absweights(i, j, s):
        # Skip line with orbital names        
        gl.readline()
        
        for k in xrange(dim):
            # Fetch entire orbital weight block
            data = np.fromfile(gl, sep=" ", count=nions*(norbs+2))
            # Cast it into tabular shape
            data = data.reshape((nions, norbs+2))
            
            # Discard first and last columns and store weights
            weights[i,:,j,k,s] = data[:,1:-1].flatten()
            
            # Skip line with the totals
            gl.readline()
            
    # Check whether phase information is included
    if '+ phase' in header_1:
        # Allocate storage for phases
        phases = np.zeros((npoints, nions*norbs, nbands, 2), complex)

        # Declare nested function that handles 
        # parsing of complex weights
        def get_weights(i, j, s):
            # Read abs values of weights
            get_absweights(i, j, s)
            
            # Skip line with orbital names
            gl.readline()
            
            # Fetch entire phase block
            data = np.fromfile(gl, sep=" ", count=2*nions*(norbs+1))
            # Cast it into tabular shape
            data = data.reshape((2*nions, norbs+1))

            # Discard first column and store real and imaginary
            # parts respectively
            phases[i,:,j,s] =  data[::2,1:].flatten()
            phases[i,:,j,s] += 1j*data[1::2,1:].flatten()
    else:
        # Phases are None in this case
        phases = None
        
        # In this case we just have absolutes of weights
        # No need for a new function
        get_weights = get_absweights
    
    # Go back to the beginning of the first k-point
    gl.seek(start)
    
    for i in xrange(npoints):
        # Parse k-point coordinates
        k_line = gl.readline().split()
        
        kpoints[i] = [float(k_line[c]) for c in [3, 4, 5]]
        kweights[i] = float(k_line[-1])
        
        for j in xrange(nbands):
            # Parse band energy
            band_line = gl.readline().split()

            bands[i, j, 0] = float(band_line[4])
            occupancies[i, j, 0] = float(band_line[-1])
            
            # Parse orbital weights
            get_weights(i, j, 0)
    
    # Seek now for the second spin component
    res = gl.readline(False)
    
    # If there is no second spin component, finish by
    # returning just the first component
    if res is None:
        # Trim the second spin component
        bands = bands[:,:,:1]
        occupancies = occupancies[:,:,:1]
        weights = weights[:,:,:,:dim,:1]
        
        if phases is not None:
            phases = phases[:,:,:,:1]
        
        return [orbitals, kpoints, kweights, bands, occupancies, \
            weights, phases]
    
    # Otherwise, read bands and weights for the second component
    for i in xrange(npoints):
        # Skip k-point coordinates
        gl.readline()
        
        for j in xrange(nbands):
            # Parse band energy
            band_line = gl.readline().split()

            bands[i, j, 1] = float(band_line[4])
            occupancies[i, j, 1] = float(band_line[-1])
            
            # Parse orbital weights
            get_weights(i, j, 1)
    
    return [orbitals, kpoints, kweights, bands, occupancies, \
        weights[:,:,:,:dim,:], phases]


