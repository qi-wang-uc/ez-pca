#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "../include/config.h"

std::string concat_str(std::string inp_str) {
	auto npos_comment = inp_str.find_first_of('#');
	return inp_str.substr(0, npos_comment);
}

bool Config::read_config(std::string inp_name) {
        std::cout << "ReadConfig> Reading Configuration from file ["
	          << inp_name << "]" << std::endl;
        std::ifstream inp_file(inp_name);
        if(!inp_file.is_open()) {
		std::cout << "Config file not found!" << std::endl;
		return false;
	}
        std::string each_line;
	std::stringstream each_stream;
        while(std::getline(inp_file, each_line)) {
		std::string prm, arg;
                if(each_line.empty() || each_line[0]=='#') continue;
		each_stream.clear();
		each_stream.str(concat_str(each_line));
		each_stream >> prm >> arg;
                if(!prm.compare("job_name")) { this->job_name=arg;}
		if(!prm.compare("psf_name")) { this->psf_name=arg;}
		if(!prm.compare("dcd_name")) { this->dcd_name=arg;}
                if(!prm.compare("num_of_pc")){ this->num_of_pc=std::atoi(arg.c_str());}
        }
        inp_file.close();
	std::cout << "ReadConfig> Done. " << std::endl << std::endl;
	if(this->job_name.empty() || this->psf_name.empty() || this->dcd_name.empty()) {
		return false;
	}
        return true;
}