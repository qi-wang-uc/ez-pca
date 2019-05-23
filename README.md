# EZ-PCA
Principal Component Analysis of Molecular Dynamics Trajectories.

## Example input
```
# Input file for EZ-PCA

job_name  test-pca      # Output file name will be $job_name_PCA.dat
psf_name  test-mol.psf  # protein structure file
dcd_name  test-mol.dcd  # trajectory file of Cartesian coordinates
num_of_pc 10            # Number of princpal components to write.
```

## Notes
A short demo of before and after applying the alignment:
- Before alignment (atomistic representation), protein is drifting in water

<img src="demo/demo-allatom.gif" width="40%" height="40%" />

- After alignment (coarse-grained representation using CA backbone atoms)

<img src="demo/demo-ca.gif" width="40%" height="40%" />

The top 3 principal components are:

<img src="demo/top-pcs.png" width="80%" height="80%" />
