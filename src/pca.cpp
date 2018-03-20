#include <iostream>
#include <fstream>
#include "../include/main.h"
#include "../include/pca.h"
#include "../include/dcd.h"

bool PCA::build_corr(const Config& config, const PSF& psf) {
	std::cout << "BuildCorr> Building correlation matrix." << std::endl;
	DCD inp_dcd; DCD_Info dcd_info;
	// 0-th read: Get header info
	if(!inp_dcd.read_dcdheader(config.dcd_name, dcd_info)) return false;
	if(dcd_info.n_atom != psf.n_atom) {
		std::cout << "ERROR> Numbers of atoms in PSF and DCD files do not match." << std::endl;
		return false;
	}
	// Till now the dcd file is guaranteed to exist with compatible format.
	// 1st read: Calculate average structure.
	const integer N_CAs  = psf.n_CAs;	// For resizing coor arrays
	this->coor_CA_avrg.resize(3*N_CAs);
	this->coor_CA_frame.resize(3*N_CAs);
	this->coor_buff.resize(dcd_info.sz_frame/sizeof(float));
	std::ifstream inp_file(config.dcd_name, std::ios::binary);
	const integer NFRAME = dcd_info.n_frame;
	for (integer iframe=0; iframe < NFRAME; iframe++) {
		inp_dcd.read_dcdframe(inp_file, coor_buff, iframe, dcd_info,
			psf.index_CAs, this->coor_CA_frame);
	}
	inp_file.close();
	// 2nd read: Align each frame wrt average structure, write aligned DCD and calculate correlation matrix.

	return true;
}

void PCA::diag_corr() {

}

void PCA::write_pca() {

}