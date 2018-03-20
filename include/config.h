#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include "main.h"

struct Config {
	std::string job_name;
	std::string psf_name;
	std::string dcd_name;
	integer	num_of_pc = 10;
	bool read_config(std::string inp_name);
};

std::string concat_str(std::string inp_str);


#endif
