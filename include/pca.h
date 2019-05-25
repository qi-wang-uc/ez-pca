#ifndef PCA_H
#define PCA_H

#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
// #include <Eigen/Eigenvalues>
#include "main.h"
#include "config.h"
#include "dcd.h"
#include "psf.h"

struct PCA {
    Coor_Sets coor_CA_avrg;
    Coor_Sets coor_CA_frame;
    std::string aligned_dcd_name;        // file name of aligned trajectory file
    // std::vector<float> C;                //[3Nx3N] array of correlation matrix
    // std::vector<float> E;                //[3Nx3N] array of eigenvector matrix
    Eigen::MatrixXf C;  //[3Nx3N] array of correlation matrix
    bool align_coor(const Config& config, const PSF& psf);
    void build_corr();
    void diag_corr(const std::string& job_name, const size_t& num_of_pc);
};

#endif