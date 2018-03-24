# ez-pca
Principal Component Analysis of Molecular Dynamics Trajectories

Program under development...

### TODO list:
- [X] Re-organize the program for easier maintainance
- [X] Read protein structure file
- [X] Read binary trajectory file
- [X] Apply proper algorithm to fit each frame to a reference (or average) structure 
- [ ] Build covariance matrix of Cartesian coordinates
- [ ] Matrix diagonalization
- [ ] Change coordinate data structure to avoid format conversion back and forth.

A short demo of before and after applying the alignment:
- Before alignment (atomistic representation)
<img src="https://github.com/wangqi1990uc/ez-pca/blob/master/demo-allatom.gif" width="40%" height="40%" />

- After alignment (coarse-grained representation using CA backbone atoms)
<img src="https://github.com/wangqi1990uc/ez-pca/blob/master/demo-ca.gif" width="40%" height="40%" />

