#include <iostream>
#include "../include/main.h"
#include "../include/align.h"

void align(const std::vector<real> ref_coor, std::vector<real> inp_coor, const integer dim3) {
	// 1. Remove displacement of center (of mass)
	integer natom = dim3/3;
	real xref = 0.0; real yref = 0.0; real zref = 0.0;
	real xinp = 0.0; real yinp = 0.0; real zinp = 0.0;
	for (integer i=0; i<natom; i++) {
		xref += ref_coor[i*3+0];
		yref += ref_coor[i*3+1];
		zref += ref_coor[i*3+2];
		xinp += inp_coor[i*3+0];
		yinp += inp_coor[i*3+1];
		zinp += inp_coor[i*3+2];
	}
	xref = xref/natom; yref = yref/natom; zref = zref/natom;
	xinp = xinp/natom; yinp = yinp/natom; zinp = zinp/natom;
	real dx = xref - xinp;  real dy = yref - yinp;  real dz = zref - zinp;
	for (integer i=0; i<natom; i++) {
		inp_coor[i*3+0] += dx;
		inp_coor[i*3+1] += dy;
		inp_coor[i*3+2] += dz;
	}
	// 2. Perform best rotation such that RMS is minimized.
	
}
