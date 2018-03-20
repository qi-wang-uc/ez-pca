#ifndef PCA_H
#define PCA_H

#include <vector>
// #include <string>
#include "main.h"
#include "config.h"
#include "psf.h"

struct PCA {
    std::vector<float> coor_buff;       //[FrameSize] array to store data of each frame
    std::vector<real> coor_CA_avrg;     //[Nx3] array of average coor of CAs
    std::vector<real> coor_CA_frame;    //[Nx3] array of frame coor of CAs
    std::vector<real> C;                //[3Nx3N] array of correlation matrix
    std::vector<real> E;                //[3Nx3N] array of eigenvector matrix

    bool build_corr(const Config& config, const PSF& psf);
    void diag_corr();
    void write_pca();
};

#endif
