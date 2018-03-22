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

	std::string out_dcdname = config.job_name + "_CA_aligned.dcd";
	std::ofstream out_file(out_dcdname, std::ios::binary);
	out_dcd.write_dcdheader(out_file, psf.n_CAs, dcd_info, dcd_pads);
	out_dcd.xcoor.resize(psf.n_CAs);
	out_dcd.ycoor.resize(psf.n_CAs);
	out_dcd.zcoor.resize(psf.n_CAs);

	// Till now the dcd file is guaranteed to exist with compatible format.
	// 1st read: Calculate average structure.
	const integer N_CAs  = psf.n_CAs;	// For resizing coor arrays
	this->coor_CA_avrg.resize(3*N_CAs);
	this->coor_CA_frame.resize(3*N_CAs);
	this->coor_buff.resize(dcd_info.sz_frame/sizeof(float));

	std::ifstream inp_file(config.dcd_name, std::ios::binary);
	const integer NFRAME = dcd_info.n_frame;
	for (integer iframe=0; iframe < NFRAME; iframe++) {
		// Get coordinates into coor_CA_frame;
		inp_dcd.read_dcdframe(inp_file, coor_buff, iframe, dcd_info, psf.index_CAs, this->coor_CA_frame);
		// Add to coor_CA_avrg.
		std::transform(this->coor_CA_avrg.begin(), this->coor_CA_avrg.end(), 
			this->coor_CA_frame.begin(), this->coor_CA_avrg.begin(), std::plus<real>());
	}
	// Devided by NFRAME to get average coordinates.
	std::transform(this->coor_CA_avrg.begin(), this->coor_CA_avrg.end(), 
		this->coor_CA_avrg.begin(), std::bind1st(std::multiplies<real>(),real(1)/NFRAME));
	inp_file.close();
	
	// 2nd read: Align each frame wrt average structure, write aligned DCD and calculate correlation matrix.
	inp_file.open(config.dcd_name, std::ios::binary);
	for (integer iframe=0; iframe < NFRAME; iframe++) {
		// Get coordinates into coor_CA_frame;
		inp_dcd.read_dcdframe(inp_file, coor_buff, iframe, dcd_info, psf.index_CAs, this->coor_CA_frame);
		align_translate(this->coor_CA_avrg, this->coor_CA_frame, N_CAs);
		align_rotate(this->coor_CA_avrg, this->coor_CA_frame, N_CAs);
		// Write aligned coordinates of CAs.
		out_dcd.convert2dcdcoor(this->coor_CA_frame, N_CAs);
		out_dcd.write_dcdframe(out_file, N_CAs, dcd_pads);
	}
	inp_file.close();
	out_file.close();

	return true;
}

void PCA::diag_corr() {

}

void PCA::write_pca() {

}