#ifndef PCA_H
#define PCA_H

#include <vector>
#include "main.h"
#include "config.h"
#include "dcd.h"
#include "psf.h"

struct PCA {
    Coor_Sets coor_CA_avrg;
    Coor_Sets coor_CA_frame;
    std::vector<float> C;                //[3Nx3N] array of correlation matrix
    std::vector<float> E;                //[3Nx3N] array of eigenvector matrix

    bool build_corr(const Config& config, const PSF& psf);
    void diag_corr();
    void write_pca();
};

#endif
