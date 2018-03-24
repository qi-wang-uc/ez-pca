#ifndef DCD_H
#define DCD_H

#include <vector>
#include <string>
#include <algorithm>
#include <valarray>
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
    std::string dcd_header2;    // [16] NATOM padded by dummy size_ts on both sides

    // Obtained from header.
    size_t n_atom  = 0;
    size_t n_frame = 0;
    size_t q_cell  = 0;
    
    // Caclulated based on above.
    size_t x_offset = 0;
	size_t y_offset = 0;
	size_t z_offset = 0;
	size_t sz_frame = 0;
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
struct Coor_Sets {
    std::valarray<float> xcoor;
    std::valarray<float> ycoor;
    std::valarray<float> zcoor;
    Coor_Sets () {}
    Coor_Sets (const std::valarray<float>& x, const std::valarray<float>& y, const std::valarray<float>& z):
        xcoor(x), ycoor(y), zcoor(z) {}
    Coor_Sets& operator += (const Coor_Sets& rhs) {
        this->xcoor += rhs.xcoor;
        this->ycoor += rhs.ycoor;
        this->zcoor += rhs.zcoor;
        return *this;
    }
    Coor_Sets& operator /= (const size_t& rhs) {
        this->xcoor /= rhs;
        this->ycoor /= rhs;
        this->zcoor /= rhs;
        return *this;
    }
    void resize(size_t dim);
};

struct DCD {
    // [FRAME_SIZE] buffer array
    std::vector<float> frame_buff;  

    // API for input dcd files
    bool read_dcdheader(const std::string& inp_name, DCD_Info& dcd_info);
    void read_dcdframe(std::ifstream& inp_file, const size_t& iframe, 
                        const DCD_Info& dcd_info, const std::vector<size_t>& index_CAs,
                        Coor_Sets& coor_sets);

    // API for out dcd files                        
    void write_dcdheader(const size_t& N_atom, const DCD_Info& dcd_info, 
                         const DCD_Pads& dcd_pads,std::ofstream& out_file);
    void write_dcdframe(const Coor_Sets& coor_sets, const size_t& N_atom, 
                        const DCD_Pads& dcd_pads, std::ofstream& out_file);
};

#endif