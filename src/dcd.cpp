#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include "../include/dcd.h"

void Coor_Sets::initialize(size_t dim, float val) {
    this->xcoor.resize(dim, val);
    this->ycoor.resize(dim, val);
    this->zcoor.resize(dim, val);
}

bool DCD::read_dcdheader(const std::string& inp_name, DCD_Info& dcd_info) {
    std::cout << "ReadDCD> Reading dcd header info from [" << inp_name << "]" << std::endl;
    std::ifstream inp_file(inp_name, std::ios::binary);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> Cannot open DCD file." << std::endl;
        return false;
    }
    char dcd_head1[100];
    char dcd_head2[16];
    char title[80];
    inp_file.read(dcd_head1, 100); dcd_info.dcd_header1.assign(dcd_head1, 100);
    if(strncmp(&dcd_head1[4], "CORD", 4)!=0) {
        std::cout << "ERROR> Wrong DCD format" << std::endl;
        return false;
    }
    int n_file = *(int*)(&dcd_head1[8]);   dcd_info.n_frame = n_file;
	int n_priv = *(int*)(&dcd_head1[12]);
	int n_savc = *(int*)(&dcd_head1[16]);
	int n_step = *(int*)(&dcd_head1[20]);
	float delta = *(float*)(&dcd_head1[44]);
	int q_cell = *(int*)(&dcd_head1[48]); dcd_info.q_cell = q_cell;
	int c24tag = *(int*)(&dcd_head1[84]);
	if(c24tag!=24) std::cout << "ReadDCD> NOT NAMD trajectory format" << std::endl;
    
    std::cout << "ReadDCD> After reading, the following information were found:" << std::endl
              << "ReadDCD> NFILE=" << n_file << " NPRIV=" << n_priv
              <<" NSAVC=" << n_savc << " NSTEP=" << n_step << std::endl
              << "ReadDCD> DELTA=" << delta
              << (q_cell==1?"  PBC cells detected.":"  PBC cells NOT detected.") << std::endl;
    
    inp_file.read(title, 80); dcd_info.dcd_remark1.assign(title, 80);
	std::cout << "ReadDCD> " << title << std::endl;
	inp_file.read(title, 80); dcd_info.dcd_remark2.assign(title, 80);
	std::cout << "ReadDCD> " << title << std::endl;
	inp_file.read(dcd_head2, 16);
	int n_atom = *(int*)(&dcd_head2[8]); dcd_info.n_atom = n_atom;
	std::cout << "ReadDCD> (" << n_atom << ") atoms found in trajectory file." << std::endl;
    dcd_info.x_offset = (q_cell==1) ? 15           : 1;	
	dcd_info.y_offset = (q_cell==1) ? n_atom+17   : n_atom + 3;
	dcd_info.z_offset = (q_cell==1) ? 2*n_atom+19 : 2*n_atom+5;
	dcd_info.sz_frame = (q_cell==1) ? (3*(4*n_atom+8)+56) : (3*(4*n_atom+8));
    inp_file.close();
    return true;
}

void DCD::read_dcdframe(std::ifstream& inp_file,
                        const size_t& iframe, const DCD_Info& dcd_info, 
                        const std::vector<size_t>& index_CAs, Coor_Sets& coor_sets) {
	// Read coordinates of a single frame. If PBC cells detected, use (56) as an offset of each frame.
	// Each frame is organized as follows: 
	// cell_offset + coor_pad + xcoor + 2*coor_pad + y_coor + 2*coor_pad + z_coor + coor_pad.
	// the corresponding size configuration is:
	// (14 or 0)   +    (1)      + n_atom+      (2)      + n_atom+      (2)      + natom  +    (1).
	// Thus, the Cartesian coordinates of the (i-1)-th atom are: (x_offset+i, y_offset+i, z_offset+i).
    if(0==iframe) {
        // Rewind dcd file and skip header information.(100+80+80+16)
        inp_file.seekg(276, std::ios::beg);	
        this->frame_buff.resize(dcd_info.sz_frame/sizeof(float));
    }
    inp_file.read(reinterpret_cast<char*>(&this->frame_buff[0]), dcd_info.sz_frame);
	for(size_t i=0; i<index_CAs.size(); i++){
		coor_sets.xcoor[i] = this->frame_buff[dcd_info.x_offset-1 + index_CAs[i]];
		coor_sets.ycoor[i] = this->frame_buff[dcd_info.y_offset-1 + index_CAs[i]];
		coor_sets.zcoor[i] = this->frame_buff[dcd_info.z_offset-1 + index_CAs[i]];
	}
}

void DCD::write_dcdheader(const size_t& N_atom, const DCD_Info& dcd_info, const DCD_Pads& dcd_pads, std::ofstream& out_file) {
    out_file.write(dcd_info.dcd_header1.c_str(), 100);
    out_file.write(dcd_info.dcd_remark1.c_str(), 80);
    out_file.write(dcd_info.dcd_remark2.c_str(), 80);
    out_file.write(reinterpret_cast<const char*>(&dcd_pads.pad164), 4);
    out_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
    out_file.write(reinterpret_cast<const char*>(&N_atom), 4);
    out_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
}

void DCD::write_dcdframe(const Coor_Sets& coor_sets, const size_t& N_atom, const DCD_Pads& dcd_pads, std::ofstream& out_file) {
    const int pad4N = 4*N_atom;
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(&coor_sets.xcoor[0]), pad4N);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(&coor_sets.ycoor[0]), pad4N);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(&coor_sets.zcoor[0]), pad4N);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
}