#ifndef ALIGN_H
#define ALIGN_H

#include <vector>
#include "dcd.h"

void align_translate(const Coor_Sets& ref_coor_sets, const size_t& Natom, Coor_Sets& inp_coor_sets);

void align_rotate(const Coor_Sets& ref_coor_sets, const size_t& Natom, Coor_Sets& inp_coor_sets);

#endif