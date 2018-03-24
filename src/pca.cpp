#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include "../include/main.h"
#include "../include/pca.h"
#include "../include/dcd.h"
#include "../include/align.h"

bool PCA::build_corr(const Config& config, const PSF& psf) {
	std::cout << "BuildCorr> Building correlation matrix." << std::endl;
	DCD inp_dcd, out_dcd;	// inp_dcd and out_dcd share almost same header but different coordinates.
	DCD_Info dcd_info; 		// info shared by inp_dcd and out_dcd.
	DCD_Pads dcd_pads;		// pads shared by inp_dcd and out_dcd, but only used by out_dcd.

	// 0-th read: Get header info and write corresponding header of CA-based trajectory.
	if(!inp_dcd.read_dcdheader(config.dcd_name, dcd_info)) return false;
	if(dcd_info.n_atom != psf.n_atom) {
		std::cout << "ERROR> Numbers of atoms in PSF and DCD files do not match." << std::endl;
		return false;
	}

	// Resize coor arrays
	const size_t N_CAs  = psf.n_CAs;
	const size_t NFRAME = dcd_info.n_frame;
	std::string out_dcdname = config.job_name + "_CA_aligned.dcd";
	std::ofstream out_file(out_dcdname, std::ios::binary);
	out_dcd.write_dcdheader(N_CAs, dcd_info, dcd_pads, out_file);

	// 1st read: Calculate average structure.
	this->coor_CA_avrg.resize(N_CAs);
	this->coor_CA_frame.resize(N_CAs);

	std::ifstream inp_file(config.dcd_name, std::ios::binary);
	for (size_t iframe=0; iframe < NFRAME; iframe++) {
		inp_dcd.read_dcdframe(inp_file, iframe, dcd_info, psf.index_CAs, this->coor_CA_frame);
		this->coor_CA_avrg += this->coor_CA_frame;		
		break;	// Currently only align wrt the 1st frame
	}
	// this->coor_CA_avrg /= NFRAME;
	inp_file.close();
	
	// 2nd read: Align each frame wrt average structure, write aligned DCD and calculate correlation matrix.
	inp_file.open(config.dcd_name, std::ios::binary);
	for (size_t iframe=0; iframe < NFRAME; iframe++) {
		inp_dcd.read_dcdframe(inp_file, iframe, dcd_info, psf.index_CAs, this->coor_CA_frame);
		align_translate(this->coor_CA_avrg,  N_CAs, this->coor_CA_frame);
		align_rotate(this->coor_CA_avrg, N_CAs, this->coor_CA_frame);
		out_dcd.write_dcdframe(this->coor_CA_frame, N_CAs, dcd_pads, out_file);
	}
	inp_file.close();
	out_file.close();

	return true;
}

void PCA::diag_corr() {

}

void PCA::write_pca() {

}