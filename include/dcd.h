#ifndef DCD_H
#define DCD_H

#include <vector>
#include <string>
#include "main.h"

// For writing DCD files.
struct DCD_Pads {
    const int pad0   = 0;
    const int pad2   = 2;
    const int pad4   = 4;
    const int pad24  = 24;
    const int pad84  = 84;
    const int pad164 = 164;  
};

struct DCD_Info {
    // Header info, retrieved when reading header (except dcd_header2).
    std::string dcd_header1;    // [100] "CORD" and info of dcd stats.
    std::string dcd_remark1;    // [80] "REMARK CREATED BY ...
    std::string dcd_remark2;    // [80] "REMARK" + $DATE + $USER 
    std::string dcd_header2;    // [16] NATOM padded by dummy integers on both sides

    // Obtained from header.
    integer n_atom  = 0;
    integer n_frame = 0;
    integer q_cell  = 0;
    // Caclulated based on above.
    integer x_offset = 0;
	integer y_offset = 0;
	integer z_offset = 0;
	integer sz_frame = 0;
    // DCD_Info() {}
    // DCD_Info(integer na, integer nf, integer qc, integer x, integer y, integer z, integer sz)
    //     : n_atom(na), n_frame(nf), q_cell(qc), x_offset(x), y_offset(y), z_offset(z), sz_frame(sz)
    //     {}
};

/* The logic of this part is a little confusing.
- For trajectory reading, we need to read input_traj 3 times: 
    (1a) - 1st time to get header info to make sure it's compatible with PSF and retrieve NFILE as number of frames.
    (1b) - 2nd read to loop through each frame in order to calculate average coordinate.
    (1c) - 3rd read to visit each frame again in order to build correlation matrix.
- For trajectory writing, we need to write output traj 2 times:
    (2a) - 1st write to copy the header info of input traj, swapping NATOM with N_CAs.
    (2b) - 2nd write to append the coordinates of each frame.
- Thus, (1a) and (2a) can be done subsequently. (1b) should be done before (1c) and (2b) since at this time 
    we don't have the average coordinate as reference to calclualte correlation or align to. After (1b), (1c) and (2b)
    can be done subsequently.
*/

// Actually DCD_IO
struct DCD {    // TODO: Make DCD a parent class and specify input or output as children class.
    // Temporary coordinate arrays for input or output.
    std::vector<real> xcoor;
    std::vector<real> ycoor;
    std::vector<real> zcoor;

    // API for input dcd files
    bool read_dcdheader(const std::string& inp_name, DCD_Info& dcd_info);
    void read_dcdframe(std::ifstream& inp_file, std::vector<real>& coor_buff,
                        const integer& iframe, const DCD_Info& dcd_info, 
                        const std::vector<integer>& index_CAs, std::vector<real>& coor_frame);
    // void convert2regcoor()
    // API for out dcd files                        
    void write_dcdheader(std::ofstream& out_file, const integer& N_atom, const DCD_Info& dcd_info, const DCD_Pads& dcd_pads);
    void write_dcdframe(std::ofstream& out_file, const integer& N_atom, const DCD_Pads& dcd_pads);

    // Covert coordinates from regular format: 
    //   (X1, Y1, Z1), (X2, Y2, Z2), ..., (Xn, Yn, Zn) 
    // to dcd format:
    //   (X1, X2, ..., Xn), (Y1, Y2, ..., Yn), (Z1, Z2, ..., Zn)
    void convert2dcdcoor(const std::vector<real>& inp_coor, const integer& N_atom);
};

#endif
