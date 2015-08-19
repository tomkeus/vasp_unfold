# Bandstructure unfolding for VASP

Here you can download vasp_unfold, a Python script which you can use to unfold the bandstructures obtained with VASP. This code is based on the method discussed in [Ref. 1](#ref_1).

You can use the code in whatever way you see fit, but if it is used to produce the data for the publication, please cite [Ref. 1](#ref_1).

The code is given as is, ie. there is no guarantee that it will work. However, if you have some problems, or the code does not behave in a way you expect it to, I encourage you to report the problem to my [e-mail](mailto:tomic@itp.uni-frankfurt.de) and I will try to fix it as soon as possible.

## Installation

Dependencies are

* [Python](http://www.python.org/downloads/) ofcourse. The script requires Python2.7 or greater
* [NumPy](http://www.scipy.org/scipylib/download.html) library for Python

No special installation is required. Just place it wherever it suits you and run it. 

## Usage

vasp_unfold works as a postprocessing tool. The unfolding is achieved through manipulation of orbital weights. Basically, in the unfolded bandstructure, orbital weights of some bands are set to zero. This means that **the effect of unfolding can only be seen with the fatband plots**. In order for vasp_unfold to be able to work, the phase information must be present in the PROCAR file which means that **LORBIT=12 must be set in the INCAR file**.

The workflow is as follows: once the bandstructure is calculated with VASP, the vasp_unfold is operated on the PROCAR file obtained in the bandstructure calculation, whereby a number of PROCAR files is produced, containing the unfolded bandstructure. Every PROCAR file generated in the output contains projection to one of the irreps as explained in [Ref. 1](#ref_1).

The command line usage of vasp_unfold is following

```
usage: vasp_unfold [-h] [--tgen SX,SY,SZ] [--out OUT] [--eps EPS]
                   [--all-irreps] [--check-mapping]
                   poscar procar
 ```

where the parameters are

```
--tgen           Fractional translation generator
--out            Output filename prefix
--eps            Numerical tolerance for position discrimination
--all-irreps     Write all irreps from the unfolding
--check-mapping  Verify if fractional translations map atoms one-to-one
poscar           Location of POSCAR file
procar           Location of PROCAR file
```

Let us say that we have a 2x3x1 supercell. In total we have six unit cells contained within the supercell. To every unit cell in the supercell, corresponds one fractional translation (relative to the supercell), so that we have six fractional translations in total. These six fractional translations can be generated from two fractional translations given by vectors (1/2,0,0) and (0,1/3,0). We can then unfold the supercell bandstructure with 

```
vasp_unfold --tgen 1/2,0,0 --tgen 0,1/3,0 POSCAR PROCAR 
```

The unfolded bandstructures will be located in PROCAR.irrep.0 file. In case --all-irreps flag was specified, the unfolded bandstructure will be located in PROCAR.irrep.0 through PROCAR.irrep.5 files. 

**NOTE 1**: no whitespace is allowed in the fractional translation generator specification. Also, the components can be either 0, or 1/N, where N is an integer. Floating point values are not allowed. 

**NOTE 2**: in case of spin polarized non-collinear calculation only spin-up and spin-down orbital weight totals will be unfolded. Orbital weight components for x,y and z component of spin won't be unfolded.

**NOTE 3**: do not enable --check-mapping flag if your structure has vacancies or excess atoms, since in this case fractional translations do not map every atom onto some other atom.
 
## Resolving the issues with the code

Here is a little advice pertaining to the "Translations are not one-to-one" error when --chek-mapping flag is enabled. This problem arises because the vasp_unfold script tries to figure out which atoms are mapped onto which atoms under the action of the fractional translations. If the supercell would be perfectly symmetrical under the fractional translations, this issue would not occur. However, in real life, the supercell will usually break this translational symmetry which means that atoms wont be mapped exactly onto each other by the fractional translations.

The resolution to this problem is quite simple. If increase of the numerical tolerance parameter does not give any result, simply use the POSCAR file where the atomic positions are restored to their symmetryc sites. Practically, if your supercell is a result of the structural optimization, use the initial POSCAR file instead of the final one. 

## References

1. <a name="ref_1"></a> M. Tomić, H. O. Jeschke and R. Valentí - Unfolding of the electronic structure through the induced representations of space groups: Application to Fe-based superconductors, [Phys. Rev. **B** 90, 195121  (2014)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.195121) ([arXiv preprint](http://arxiv.org/abs/1408.2258))
