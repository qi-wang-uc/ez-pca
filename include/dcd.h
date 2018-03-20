#ifndef DCD_H
#define DCD_H

#include <vector>
#include <string>
#include "main.h"

struct DCD_Info {
    integer n_atom  = 0;
    integer n_frame = 0;
    integer q_cell  = 0;
    // The following info are caclulated based on above.
    integer x_offset = 0;
	integer y_offset = 0;
	integer z_offset = 0;
	integer sz_frame = 0;
    DCD_Info() {}
    DCD_Info(integer na, integer nf, integer qc, integer x, integer y, integer z, integer sz)
        : n_atom(na), n_frame(nf), q_cell(qc), x_offset(x), y_offset(y), z_offset(z), sz_frame(sz)
        {}
};

struct DCD {
    std::string dcd_header1;    // [100] "CORD" and info of dcd stats.
    std::string dcd_remark1;    // [80] "REMARK CREATED BY ...
    std::string dcd_remark2;    // [80] "REMARK" + $DATE + $USER 
    std::string dcd_header2;    // [16] NATOM padded by dummy integers on both sides

    int n_file = 0;
    int n_priv = 0;
    int n_savc = 0;
    int n_step = 0;
    int q_cell = 0;
    int c24tag = 0;
    int n_atom = 0;
    float delta = 0.0;

    std::vector<float> xcoor;
    std::vector<float> ycoor;
    std::vector<float> zcoor;

    bool read_dcdheader(const std::string& inp_name, DCD_Info& out_dcdinfo);
    void read_dcdframe(std::ifstream& inp_file, std::vector<float>& coor_buff,
                        const integer& iframe, const DCD_Info& dcd_info, 
                        const std::vector<integer>& index_CAs, std::vector<float>& coor_frame);
    void write_dcdheader();
    void write_dcdframe();
};

#endif
