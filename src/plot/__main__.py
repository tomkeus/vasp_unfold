#===========================================================
#
#  PROJECT: vasp_unfold/fatplot
#  FILE:    __main__.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.32
#  DATE:    Apr 25th 2017
#
#===========================================================


import os
import numpy as np
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot

# Handle color arguments
def color(string):
    try:
        # Try to convert to tuple
        return eval(string)
    except:
        # If cannot convert to tuple just
        # return the original string
        return string

desc_str = '''Simple program used to quickly plot the orbital weights
of the band structure contained in the specified PROCAR file.
'''

parser = argparse.ArgumentParser(prog='fatplot', description = desc_str)

parser.add_argument('procar', type=str, help='POSCAR file')
parser.add_argument('output', type=str, help='Filename for the plot. It accepts '
                    'all image formats supported by matplotlib.')

parser.add_argument('--figsize', type=eval, default=(6, 4),
                    help='Figure size in inches formatted as width,height '
                    '(no spaces allowed). Default is 6,4')
parser.add_argument('--dpi', type=float, default=300,
                    help='Plot resolution in dots per inch. Default is 300.')
parser.add_argument('--marker', type=str, default='o', 
                    help='Marker for the fatband plot. Default is o.')
parser.add_argument('--markersize', type=float, default=20,
                    help='Marker size. Default is 20.')
parser.add_argument('--color', type=color, default='"steelblue"',
                    help='Color for the marker. It accepts any color specification '
                    'accepted by matplotlib. Color specified as r,g,b tuple should ' 
                    'not have any spaces.')
parser.add_argument('--elim', type=eval, default=(1,-1),
                    help='Energy range for the plot, specified as emin,emax '
                    '(no spaces allowed). Default is entire band range.')
parser.add_argument('--efermi', type=float, default=0,
                    help='Fermi energy. Default is 0.')
parser.add_argument('--pow', type=float, default=1,
                    help='Raise orbital weights to the specified integer power. '
                    'Powers larger than 1 help to filter out the ghost bands '
                    'in the unfolded band structures.')
                    
args = parser.parse_args()


# Commands to extract the needed information from the PROCAR file
grep_npoints = 'grep -m 1 "of k\-points" {0} | tr -s " " | cut -d" " -f4'
grep_kpoints = 'grep -E "^ k\-point " {0} | tr -s " " | cut -d" " -f5-7'
grep_bands = 'grep -E "^band" {0} | tr -s " " |  cut -d" " -f5'
grep_weights = 'grep -E "^tot" {0} | tr -s " " |  cut -d" " -f11'

# Extract the number of k-points
npoints = int(os.popen(grep_npoints.format(args.procar)).read())

# Extract the k-points
kpoints = os.popen(grep_kpoints.format(args.procar))
kpoints = np.fromfile(kpoints, count=3*npoints, dtype=float, sep=' ')

# Extract the band energies
bands = os.popen(grep_bands.format(args.procar))
bands = np.fromfile(bands, count=-1, dtype=float, sep=' ')

# Extract the total orbital weights for each band
weights = os.popen(grep_weights.format(args.procar))
weights = np.fromfile(weights, count=-1, dtype=float, sep=' ')

# Figure out the number of bands
nbands = len(bands)/npoints

# Figure out the number of orbital blocks per band
# 1 for collinear calculation, 4 for non-collinear
# In non-collinear case we just need the first one
wdim = len(weights)/(npoints*nbands)

# Reshape the arrays into their proper shapes
kpoints = kpoints.reshape((npoints, 3))
bands = bands.reshape((npoints, nbands))-args.efermi
weights = weights[::wdim].reshape((npoints, nbands))

# Raise the weights to the specified power
if args.pow != 1:
    weights = np.power(weights, args.pow)


# Try to guess where the high symmetry points 
# are by looking for the repeated k-points
dk = np.zeros((npoints, 3), float)

# Displacement between i+1st and ith k-point
dk[1:] = kpoints[1:]-kpoints[:-1]

# Get the magnitude of displacements
dk = np.sqrt(np.sum(dk*dk, axis=1))

# Locate the high-symmetry points
ipoint = np.where(dk < 1e-6)[0]

# Generate plot's x-axis from the adjacent k-point
# displacement magnitudes
x = np.cumsum(dk)

# Get the plot's x-coordinates for high-symmetry k-points
xsym = x[ipoint]

# Figure out the maximum possible energy boundaries 
# based on the energy extend of all bands
e_low = np.min(bands)-1.0
e_high = np.max(bands)+1.0

plot.figure(figsize=args.figsize)

# Plot horizontal line for the Fermi level
plot.plot(x[[0,-1]], [0, 0], color='gray', zorder=-1)

# Plot vertical lines for the high-symmetry points
for xi in xsym:
    plot.plot([xi, xi], [e_low, e_high], color='gray', zorder=-1)

# Plot the weights
for bi, wi in zip(bands.T, weights.T):
    plot.scatter(x, bi, s=args.markersize*wi, marker=args.marker, 
                 color=args.color, lw=0)

# Fix x and y-axis boundaries
plot.xlim(x[0], x[-1])

if args.elim[0] < args.elim[1]:
    plot.ylim(*args.elim)
else:
    plot.ylim(e_low, e_high)

# Set x-axis ticks
plot.xticks(xsym, ['']*len(xsym))

# Save the figure
plot.savefig(args.output, dpi=args.dpi, bbox_inches='tight')


