#ifndef PSF_H
#define PSF_H

#include <vector>
#include <string>
#include "main.h"

struct PSFatom {
	size_t atomid;			// (1)-th column in PSF file.
	size_t resid;			// (3)-th column in PSF file.
	size_t unused;			// (9)-th column in PSF file.
	float charge;			// (7)-th column in PSF file.
	float mass;				// (8)-th column in PSF file.
	std::string segname;	// (2)-th column in PSF file.
	std::string resname;	// (4)-th column in PSF file.
	std::string atomname;	// (5)-th column in PSF file.
	std::string atomtype;	// (6)-th column in PSF file.
};

struct PSF {
	size_t n_atom = 0;
	size_t n_CAs   = 0;
	std::string psf_header;
	std::vector<size_t> index_CAs;
	std::vector<PSFatom> psf_atoms;	// Only CA atoms in this case.
	bool read_psf(std::string inp_name);
	bool write_psf(std::string out_name);
};

#endif
