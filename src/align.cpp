#include <iostream>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "../include/main.h"
#include "../include/align.h"

void align_translate(const std::vector<real>& ref_coor, std::vector<real>& inp_coor, const integer& Natom) {
	// 1. Remove displacement of center (of mass)
	real xref = 0.0; real yref = 0.0; real zref = 0.0;
	real xinp = 0.0; real yinp = 0.0; real zinp = 0.0;
	for (integer i=0; i<Natom; i++) {
		xref += ref_coor[i*3+0];
		yref += ref_coor[i*3+1];
		zref += ref_coor[i*3+2];
		xinp += inp_coor[i*3+0];
		yinp += inp_coor[i*3+1];
		zinp += inp_coor[i*3+2];
	}
	xref /= Natom; yref /= Natom; zref /= Natom;
	xinp /= Natom; yinp /= Natom; zinp /= Natom;
	real dx = xref - xinp;  real dy = yref - yinp;  real dz = zref - zinp;
	for (integer i=0; i<Natom; i++) {
		inp_coor[i*3+0] += dx;
		inp_coor[i*3+1] += dy;
		inp_coor[i*3+2] += dz;
	}
}

void align_rotate(const std::vector<real>& ref_coor, std::vector<real>& inp_coor, const integer& Natom) {
	// 1. Retrieve elements of matrix M to build quaternion-based matrix K, implicitly.
	real Sxx = 0.0, Sxy = 0.0, Sxz = 0.0;
	real Syx = 0.0, Syy = 0.0, Syz = 0.0;
	real Szx = 0.0, Szy = 0.0, Szz = 0.0;

	for (integer i=0; i<Natom; i++) {
		real x1 = ref_coor[i*3+0]; real y1 = ref_coor[i*3+1]; real z1 = ref_coor[i*3+2];
		real x2 = inp_coor[i*3+0]; real y2 = inp_coor[i*3+1]; real z2 = inp_coor[i*3+2];
		Sxx +=  (x1 * x2); Sxy +=  (x1 * y2); Sxz +=  (x1 * z2);
		Syx +=  (y1 * x2); Syy +=  (y1 * y2); Syz +=  (y1 * z2);
		Szx +=  (z1 * x2); Szy +=  (z1 * y2); Szz +=  (z1 * z2);
	}

	// A = USV', W=UV'
	Eigen::Matrix3f A; Eigen::Matrix3f I; I.setIdentity();

	A << Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
	Eigen::JacobiSVD<Eigen::MatrixXf> svd( A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXf W = svd.matrixU() * svd.matrixV().transpose();
	I(2,2) = W.determinant();
	Eigen::MatrixXf R = svd.matrixU() * I * svd.matrixV().transpose();

	// 4. Apply rotation matrix
    for (integer i = 0; i < Natom; i++) {
        real x = R(0,0)*inp_coor[i*3+0] + R(0,1)*inp_coor[i*3+1] + R(0,2)*inp_coor[i*3+2];
        real y = R(1,0)*inp_coor[i*3+0] + R(1,1)*inp_coor[i*3+1] + R(1,2)*inp_coor[i*3+2];
        real z = R(2,0)*inp_coor[i*3+0] + R(2,1)*inp_coor[i*3+1] + R(2,2)*inp_coor[i*3+2];

        inp_coor[i*3+0] = x;
        inp_coor[i*3+1] = y;
        inp_coor[i*3+2] = z;
    }
}
