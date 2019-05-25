#include <iostream>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "../include/main.h"
#include "../include/align.h"

void align_translate(const Coor_Sets& ref_coor_sets, const size_t& Natom, Coor_Sets& inp_coor_sets) {
    float dx = (ref_coor_sets.xcoor.sum() - inp_coor_sets.xcoor.sum())/Natom;
    float dy = (ref_coor_sets.ycoor.sum() - inp_coor_sets.ycoor.sum())/Natom;
    float dz = (ref_coor_sets.zcoor.sum() - inp_coor_sets.zcoor.sum())/Natom;
    inp_coor_sets.xcoor += dx;
    inp_coor_sets.ycoor += dy;
    inp_coor_sets.zcoor += dz;
}

void align_rotate(const Coor_Sets& ref_coor_sets, const size_t& Natom, Coor_Sets& inp_coor_sets) {
    // Build cross-covariance matrix
    float Sxx = (ref_coor_sets.xcoor*inp_coor_sets.xcoor).sum();
    float Sxy = (ref_coor_sets.xcoor*inp_coor_sets.ycoor).sum();
    float Sxz = (ref_coor_sets.xcoor*inp_coor_sets.zcoor).sum();
    float Syx = (ref_coor_sets.ycoor*inp_coor_sets.xcoor).sum();
    float Syy = (ref_coor_sets.ycoor*inp_coor_sets.ycoor).sum();
    float Syz = (ref_coor_sets.ycoor*inp_coor_sets.zcoor).sum();
    float Szx = (ref_coor_sets.zcoor*inp_coor_sets.xcoor).sum();
    float Szy = (ref_coor_sets.zcoor*inp_coor_sets.ycoor).sum();
    float Szz = (ref_coor_sets.zcoor*inp_coor_sets.zcoor).sum();

    // SVD: A = USV', W=UV' (https://en.wikipedia.org/wiki/Kabsch_algorithm)
    Eigen::Matrix3f A; Eigen::Matrix3f I; I.setIdentity();
    A << Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    Eigen::JacobiSVD<Eigen::MatrixXf> svd( A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXf W = svd.matrixU() * svd.matrixV().transpose();
    I(2,2) = W.determinant();
    Eigen::MatrixXf R = svd.matrixU() * I * svd.matrixV().transpose();

    // Apply rotation matrix
    for (size_t i = 0; i < Natom; i++) {
        float x = R(0,0)*inp_coor_sets.xcoor[i] + R(0,1)*inp_coor_sets.ycoor[i] + R(0,2)*inp_coor_sets.zcoor[i];
        float y = R(1,0)*inp_coor_sets.xcoor[i] + R(1,1)*inp_coor_sets.ycoor[i] + R(1,2)*inp_coor_sets.zcoor[i];
        float z = R(2,0)*inp_coor_sets.xcoor[i] + R(2,1)*inp_coor_sets.ycoor[i] + R(2,2)*inp_coor_sets.zcoor[i];
        inp_coor_sets.xcoor[i] = x;
        inp_coor_sets.ycoor[i] = y;
        inp_coor_sets.zcoor[i] = z;
    }
}