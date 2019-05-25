#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include "../include/psf.h"

// Any PSF file that can help to retrieve atom index information is ok.
bool PSF::read_psf(std::string inp_name) {
    std::cout << "ReadPSF> Reading protein structure from file ["
          << inp_name << "]" << std::endl;
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> PSF file not found." << std::endl;
        return false;
    }
    std::string each_line;
    std::stringstream each_stream;
    std::string placeholder1, placeholder2, placeholder3;
    bool is_NATOM = false;
    while(getline(inp_file, each_line)) {
        each_stream.clear();
        // skip empty lines and print comments
        if(each_line.empty()) continue;
        if(each_line[0]=='*') {
            std::cout << each_line << std::endl;
            continue;
        }
        // detect NATOM block
        each_stream.str(each_line);
        each_stream >> placeholder1 >> placeholder2;
        if(placeholder1=="PSF") {
            this->psf_header = each_line;
            continue;
        }
        if(placeholder2=="!NATOM") {
            this->n_atom=std::atoi(placeholder1.c_str());
            is_NATOM = true;
            continue;
        }
        // read NATOM block
        if(is_NATOM) {
            PSFatom psf_atom;
            each_stream.str(each_line);
            each_stream >> psf_atom.atomid  >> psf_atom.segname  >> psf_atom.resid 
                        >> psf_atom.resname >> psf_atom.atomname >> psf_atom.atomtype 
                        >> psf_atom.charge  >> psf_atom.mass     >> psf_atom.unused;
            // unsigned int index = std::atoi(psf_atom.atomid);
            if(psf_atom.atomname=="CA") {
                this->index_CAs.push_back(psf_atom.atomid);
                this->psf_atoms.push_back(psf_atom);
            } 
            if(psf_atom.atomid==this->n_atom) is_NATOM=false;
        }
    }
    this->n_CAs=this->index_CAs.size();
    inp_file.close();
    std::cout << "ReadPSF> After reading psf file, (" << this->n_atom 
              << ") atoms were recorded" << ", (" << this->n_CAs 
              << ") CAs were found." << std::endl << std::endl;
    return true;
}

bool PSF::write_psf(std::string inp_name) {
    std::string out_name = inp_name + "_CA.psf"	;
    std::cout << "WritePSF> Writing PSF of CAs to file: [" 
              << out_name << "] ..." << std::endl;
    std::ofstream out_file(out_name);
    if(!out_file.is_open()) return false;
    out_file << this->psf_header << std::endl << std::endl 
             << this->n_CAs << " !NATOM" << std::endl;
    for(size_t i=0; i<this->psf_atoms.size(); i++) {
        out_file << std::right << std::setw(8) << (i+1)
                 << " " << std::setw(4) << this->psf_atoms[i].segname 
                 << " " << std::left << std::setw(4) << this->psf_atoms[i].resid 
                 << " " << std::left << std::setw(4) << this->psf_atoms[i].resname 
                 << " " << std::left << std::setw(4) << this->psf_atoms[i].atomname 
                 << " " << std::left << this->psf_atoms[i].atomtype 
                 << "    0.00000       100.000           0   0.00000     -0.301140E-02" // quick & dirty
                 << std::endl;
    }
    out_file.close();
    std::cout << "WritePSF> Done." << std::endl << std::endl;
    return true;
}