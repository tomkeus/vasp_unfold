#===========================================================
#
#  PROJECT: vasp_unfold
#  FILE:    __errors__.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.4
#  DATE:    December 12th 2017
#
#===========================================================

poscar_parse_error = """Unable to parse the input PROCAR file. Please check if the 
PROCAR file is properly formatted.

PLEASE READ CAREFULLY

Beware that output formatting in VASP is broken. Output field 
widths are not set to a large enough value in the source code. 
As a result, two common types of errors occur:

 1. When a k-point has negative component, k-point output merges
    into a single string, for example

  k-point   75 :    0.25510204-0.74489796-0.00000000     weight = 0.00500000 
                              ^          ^
                            error       error
                          
    HOW TO FIX:
    
    Use sed to insert spaces between the components
    
  sed "s/\([0-9]\)\-\([0-9]\)/\\1 -\\2/g" PROCAR.broken > PROCAR.fixed
  

 2. Orbital weight/phase precision requires field width larger than 
    the field width allows. As a result, instead of orbital weight/
    phase, asterisks are written. For example
    
  ion      s     py     pz     px    dxy    dyz    dz2 dxz    dx2
    1  0.000 -0.292  0.000 -0.503  0.000  0.001  0.000 -0.001  0.000
    1  0.000 -1.286  0.000  0.049  0.000  0.003  0.000  0.000  0.000
    2  0.000  0.292  0.000  0.503  0.000  0.001  0.000 -0.001  0.000
    2  0.000  1.286  0.000 -0.049  0.000  0.003  0.000  0.000  0.000
    3  0.000 -0.003  0.000  0.006  0.000  0.002  0.000 -0.003  0.000
    3  0.000 -0.015  0.000 -0.001  0.000  0.008  0.000  0.000  0.000
    4  *****  0.003  0.000 -0.006  0.000  0.002  0.000 -0.003  0.000
    4  *****  0.015  0.000  0.001  0.000  0.008  0.000  0.000  0.000 
        ^
       errors
       
    HOW TO FIX:
   
    Use sed to replace every asterisk with a zero
   
  sed "s/\*/0/g" PROCAR > PROCAR.edited
  
  
So far no other sources of errors have been encountered. In order 
to solve them for good, either edit VASPs source code or bug VASPs
developer until they fix it.
"""
