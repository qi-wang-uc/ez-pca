#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <numeric>
#include <iomanip>
#include "../include/main.h"
#include "../include/pca.h"
#include "../include/dcd.h"
#include "../include/align.h"

bool PCA::align_coor(const Config& config, const PSF& psf) {
    
    DCD inp_dcd, out_dcd; // inp_dcd and out_dcd share almost same header but different coordinates.
    DCD_Info dcd_info;    // info shared by inp_dcd and out_dcd.
    DCD_Pads dcd_pads;    // pads shared by inp_dcd and out_dcd, but only used by out_dcd.

    // 0-th read: Get dcd info and write corresponding header of CA-based trajectory.
    if(!inp_dcd.read_dcdheader(config.dcd_name, dcd_info)) return false;
    if(dcd_info.n_atom != psf.n_atom) {
        std::cout << "ERROR> Numbers of atoms in PSF and DCD files do not match." << std::endl;
        return false;
    }

    // Initialize coor arrays
    const size_t N_CAs  = psf.n_CAs;
    const size_t NFRAME = dcd_info.n_frame;
    this->aligned_dcd_name = config.job_name + "_CA_aligned.dcd";
    std::cout << "AlignCoor> Alignning coordinates and writing to trajectory [" 
              << this->aligned_dcd_name << "]" << std::endl;
    std::ofstream out_file(this->aligned_dcd_name, std::ios::binary);
    out_dcd.write_dcdheader(N_CAs, dcd_info, dcd_pads, out_file);

    // 1st read: Read the 1st frame into average as reference.
    this->coor_CA_avrg.initialize(N_CAs, 0.0);
    this->coor_CA_frame.initialize(N_CAs, 0.0);
    std::ifstream inp_file(config.dcd_name, std::ios::binary);
    inp_dcd.read_dcdframe(inp_file, 0, dcd_info, psf.index_CAs, this->coor_CA_avrg);
    inp_file.close();
    
    // 2nd read: Align each frame wrt the 1st structure, write aligned DCD.
    inp_file.open(config.dcd_name, std::ios::binary);
    for (size_t iframe=0; iframe < NFRAME; iframe++) {
        inp_dcd.read_dcdframe(inp_file, iframe, dcd_info, psf.index_CAs, this->coor_CA_frame);
        align_translate(this->coor_CA_avrg,  N_CAs, this->coor_CA_frame);
        align_rotate(this->coor_CA_avrg, N_CAs, this->coor_CA_frame);
        out_dcd.write_dcdframe(this->coor_CA_frame, N_CAs, dcd_pads, out_file);
    }
    inp_file.close();
    out_file.close();
    std::cout << "AlignCoor> Done." << std::endl << std::endl;
    return true;
}

void PCA::build_corr() {
    std::cout << "BuildCorr> Building correlation matrix ..." << std::endl;
    DCD inp_CA_dcd;
    DCD_Info dcd_info;
    inp_CA_dcd.read_dcdheader(this->aligned_dcd_name, dcd_info);
    const size_t NFRAME = dcd_info.n_frame;
    const size_t NATOM  = dcd_info.n_atom;
    // 1st read: Calculate average structure.
    this->C.resize(3*NATOM, 3*NATOM);
    this->coor_CA_avrg.initialize(NATOM, 0.0);
    std::vector<size_t> dummy_index(NATOM);
    std::iota(dummy_index.begin(), dummy_index.end(), 1);
    std::ifstream inp_file(this->aligned_dcd_name, std::ios::binary);
    for (size_t iframe=0; iframe < NFRAME; iframe++) {
        inp_CA_dcd.read_dcdframe(inp_file, iframe, dcd_info, dummy_index, this->coor_CA_frame);
        this->coor_CA_avrg += this->coor_CA_frame;
    }
    this->coor_CA_avrg /= NFRAME;
    for(size_t iframe=0; iframe < NFRAME; iframe++) {
        inp_CA_dcd.read_dcdframe(inp_file, iframe, dcd_info, dummy_index, this->coor_CA_frame);
        for(size_t i=0; i < NATOM; i++) {
            float dxi = this->coor_CA_frame.xcoor[i] - this->coor_CA_avrg.xcoor[i];
            float dyi = this->coor_CA_frame.ycoor[i] - this->coor_CA_avrg.ycoor[i];
            float dzi = this->coor_CA_frame.zcoor[i] - this->coor_CA_avrg.zcoor[i];
            for(size_t j=0; j < NATOM; j++) {
                float dxj = this->coor_CA_frame.xcoor[j] - this->coor_CA_avrg.xcoor[j];
                float dyj = this->coor_CA_frame.ycoor[j] - this->coor_CA_avrg.ycoor[j];
                float dzj = this->coor_CA_frame.zcoor[j] - this->coor_CA_avrg.zcoor[j];
                this->C(i*3+0, j*3+0) += dxi*dxj;
                this->C(i*3+0, j*3+1) += dxi*dyj;
                this->C(i*3+0, j*3+2) += dxi*dzj;
                this->C(i*3+1, j*3+0) += dyi*dxj;
                this->C(i*3+1, j*3+1) += dyi*dyj;
                this->C(i*3+1, j*3+2) += dyi*dzj;
                this->C(i*3+2, j*3+0) += dzi*dxj;
                this->C(i*3+2, j*3+1) += dzi*dyj;
                this->C(i*3+2, j*3+2) += dzi*dzj;
            }
        }
    }
    this->C /= NFRAME;
    inp_file.close();
    std::cout << "BuildCorr> Done." << std::endl << std::endl;
}

void PCA::diag_corr(const std::string& job_name, const size_t& num_of_pc) {
    std::cout << "DiagCorr> Diagonalizing correlation matrix ..." << std::endl;
    Eigen::JacobiSVD<Eigen::MatrixXf> svd( this->C, Eigen::ComputeFullU | Eigen::ComputeFullV);
    std::string out_name = job_name + "_PCA_data.dat";
    std::cout << "DiagCorr> Writing eigenvalues and eigenvectors to file ["
              << out_name << "]" << std::endl;
    std::ofstream out_file(out_name);
    for (size_t i=0; i<num_of_pc; i++) {
        out_file << std::setw(10) << std::fixed << std::setprecision(6) 
                 << svd.singularValues()[i] << " : " << svd.matrixV().transpose().row(i)
                 << std::endl;
    }
    out_file.close();
    std::cout << "DiagCorr> Done." << std::endl;
}