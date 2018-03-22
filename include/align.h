#ifndef ALIGN_H
#define ALIGN_H

#include <vector>
#include "main.h"

void align_translate(const std::vector<real>& ref_coor, std::vector<real>& inp_coor, const integer& dim);

void align_rotate(const std::vector<real>& ref_coor, std::vector<real>& inp_coor, const integer& dim);

#endif
