#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include "../include/dcd.h"

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
    int n_file = *(unsigned int*)(&dcd_head1[8]);   dcd_info.n_frame = n_file;
	int n_priv = *(unsigned int*)(&dcd_head1[12]);
	int n_savc = *(unsigned int*)(&dcd_head1[16]);
	int n_step = *(unsigned int*)(&dcd_head1[20]);
	float delta = *(float*)(&dcd_head1[44]);
	int q_cell = *(unsigned int*)(&dcd_head1[48]); dcd_info.q_cell = q_cell;
	int c24tag = *(unsigned int*)(&dcd_head1[84]);
	if(c24tag!=24) std::cout << "WARNNING> NOT NAMD trajectory format" << std::endl;
    
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

void DCD::read_dcdframe(std::ifstream& inp_file, std::vector<real>& coor_buff,
                        const integer& iframe, const DCD_Info& dcd_info, 
                        const std::vector<integer>& index_CAs, std::vector<real>& coor_frame) {
	// Read coordinates of a single frame. If PBC cells detected, use (56) as an offset of each frame.
	// Each frame is organized as follows: 
	// cell_offset + coor_pad + xcoor + 2*coor_pad + y_coor + 2*coor_pad + z_coor + coor_pad.
	// the corresponding size configuration is:
	// (14 or 0)   +    (1)      + n_atom+      (2)      + n_atom+      (2)      + natom  +    (1).
	// Thus, the Cartesian coordinates of the (i-1)-th atom are: (x_offset+i, y_offset+i, z_offset+i).
    if(0==iframe) inp_file.seekg(276, std::ios::beg);	// Rewind dcd file and skip header information.(100+80+80+16)
    // std::cout << "ReadDCD> Processing frame " << std::setw(8) << iframe+1 << "..." << std::endl;
    inp_file.read(reinterpret_cast<char*>(coor_buff.data()), dcd_info.sz_frame);
	for(integer i=0; i<index_CAs.size(); i++){
		coor_frame[i*3+0] = coor_buff[dcd_info.x_offset + index_CAs[i]];
		coor_frame[i*3+1] = coor_buff[dcd_info.y_offset + index_CAs[i]];
		coor_frame[i*3+2] = coor_buff[dcd_info.z_offset + index_CAs[i]];
	}
}

void DCD::write_dcdheader(std::ofstream& out_file, const integer& N_atom, const DCD_Info& dcd_info, const DCD_Pads& dcd_pads) {
    out_file.write(dcd_info.dcd_header1.c_str(), 100);
    out_file.write(dcd_info.dcd_remark1.c_str(), 80);
    out_file.write(dcd_info.dcd_remark2.c_str(), 80);
    out_file.write(reinterpret_cast<const char*>(&dcd_pads.pad164), 4);
    out_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
    out_file.write(reinterpret_cast<const char*>(&N_atom), 4);
    out_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
}

void DCD::write_dcdframe(std::ofstream& out_file, const integer& N_atom, const DCD_Pads& dcd_pads) {
    const int pad4N = 4*N_atom;
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(this->xcoor.data()), pad4N);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(this->ycoor.data()), pad4N);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    out_file.write(reinterpret_cast<const char*>(this->zcoor.data()), pad4N);
    out_file.write(reinterpret_cast<const char*>(&pad4N), 4);
}

void DCD::convert2dcdcoor(const std::vector<real>& inp_coor, const integer& N_atom) {
    for (integer i=0; i<N_atom; i++) {
        this->xcoor[i] = inp_coor[i*3+0];
        this->ycoor[i] = inp_coor[i*3+1];
        this->zcoor[i] = inp_coor[i*3+2];
    }
}