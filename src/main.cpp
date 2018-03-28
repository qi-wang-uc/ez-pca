#include <iostream>
#include <string>
#include "../include/config.h"
#include "../include/psf.h"
#include "../include/pca.h"

int main(int argc, char* argv[]) {
	std::cout << "EZPCA> Principal Component Analysis of Molecular Dynamics Trajectories" << std::endl;
	if (argc < 2) {
		std::cout << "ERROR> Missing input file." << std::endl;
		return EXIT_FAILURE;
	}
	std::string inp_name = argv[1];
	Config config;
	if(!config.read_config(inp_name)) return EXIT_FAILURE;
	PSF psf;
	if(!psf.read_psf(config.psf_name)) return EXIT_FAILURE;
	if(!psf.write_psf(config.job_name))return EXIT_FAILURE;
	PCA pca;
	pca.align_coor(config, psf);
	pca.build_corr();
	pca.diag_corr(config.job_name, config.num_of_pc);
	return EXIT_SUCCESS;
}
