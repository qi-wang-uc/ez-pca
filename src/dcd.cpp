#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include "../include/dcd.h"

bool DCD::read_dcdheader(const std::string& inp_name, DCD_Info& out_dcdinfo) {
    std::cout << "ReadDCD> Reading dcd header info from [" << inp_name << "]" << std::endl;
    std::ifstream inp_file(inp_name, std::ios::binary);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> Cannot open DCD file." << std::endl;
        return false;
    }
    char dcd_head1[100];
    char dcd_head2[16];
    char title[80];
    // char test[56];  // ??? wtf is this ???
    inp_file.read(dcd_head1, 100);
    if(strncmp(&dcd_head1[4], "CORD", 4)!=0) {
        std::cout << "ERROR> Wrong DCD format" << std::endl;
        return false;
    }
    int n_file = *(unsigned int*)(&dcd_head1[8]);
	int n_priv = *(unsigned int*)(&dcd_head1[12]);
	int n_savc = *(unsigned int*)(&dcd_head1[16]);
	int n_step = *(unsigned int*)(&dcd_head1[20]);
	float delta = *(float*)(&dcd_head1[44]);
	int q_cell = *(unsigned int*)(&dcd_head1[48]);
	int c24tag = *(unsigned int*)(&dcd_head1[84]);
	if(c24tag!=24) std::cout << "WARNNING> NOT NAMD trajectory format" << std::endl;
    
    std::cout << "ReadDCD> After reading, the following information were found:" << std::endl
              << "ReadDCD> NFILE=" << n_file << " NPRIV=" << n_priv
              <<" NSAVC=" << n_savc << " NSTEP=" << n_step << std::endl
              << "ReadDCD> DELTA=" << delta
              << (q_cell==1?"  PBC cells detected.":"  PBC cells NOT detected.") << std::endl;
    
    inp_file.read(title, 80);
	std::cout << "ReadDCD> " << title << std::endl;
	inp_file.read(title, 80);
	std::cout << "ReadDCD> " << title << std::endl;
	inp_file.read(dcd_head2, 16);
	int n_atom = *(int*)(&dcd_head2[8]);
	std::cout << "ReadDCD> (" << n_atom << ") atoms found in trajectory file." << std::endl;
    integer x_offset = (q_cell==1) ? 15           : 1;	
	integer y_offset = (q_cell==1) ? n_atom+17   : n_atom + 3;
	integer z_offset = (q_cell==1) ? 2*n_atom+19 : 2*n_atom+5;
	integer sz_frame = (q_cell==1) ? (3*(4*n_atom+8)+56) : (3*(4*n_atom+8));
    out_dcdinfo = DCD_Info(n_atom, n_file, q_cell, x_offset, y_offset, z_offset, sz_frame);
    inp_file.close();
    return true;
}

void DCD::read_dcdframe(std::ifstream& inp_file, std::vector<float>& coor_buff,
                        const integer& iframe, const DCD_Info& dcd_info, 
                        const std::vector<integer>& index_CAs, std::vector<float>& coor_frame) {
	// Read coordinates of a single frame. If PBC cells detected, use (56) as an offset of each frame.
	// Each frame is organized as follows: 
	// cell_offset + coor_pad + x_coor + 2*coor_pad + y_coor + 2*coor_pad + z_coor + coor_pad.
	// the corresponding size configuration is:
	// (14 or 0)   +    (1)      + n_atom+      (2)      + n_atom+      (2)      + natom  +    (1).
	// Thus, the Cartesian coordinates of the (i-1)-th atom are: (x_offset+i, y_offset+i, z_offset+i).
    if(0==iframe) inp_file.seekg(276, std::ios::beg);	// Rewind dcd file and skip header information.
    std::cout << "ReadDCD> Processing frame " << std::setw(8) << iframe+1 << "..." << std::endl;
    inp_file.read(reinterpret_cast<char*>(coor_buff.data()), dcd_info.sz_frame);
	for(integer i=0; i<index_CAs.size(); i++){
		coor_frame[i*3+0] = coor_buff[dcd_info.x_offset + index_CAs[i]];
		coor_frame[i*3+1] = coor_buff[dcd_info.y_offset + index_CAs[i]];
		coor_frame[i*3+2] = coor_buff[dcd_info.z_offset + index_CAs[i]];
	}
}