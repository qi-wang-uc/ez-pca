#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

struct AtomsCount{
	unsigned int n_atoms;
	unsigned int n_CAs;
};

struct Params{
	string job_name;
	string psf_name;
	string dcd_name;
	unsigned int n_of_pc;
};

vector<int>   index_CAs;
vector<float> coor_buff;	// coor_buff.size() = size of each dcd frame.
vector<float> coor_CA_frame;	// coor_CA_frame.size() = 3 * N_OF_CA.
vector<float> coor_CA_avrg;	// coor_CA_avrg.size()  = 3 * N_OF_CA.
vector<vector<float> > C;	// covariance matrix, C.size() = 3*N_OF_CA x 3*N_OF_CA.
vector<vector<float> > E;	// eigenvector matrix, E.size() = 3*N_OF_CA x 3*N_OF_CA.

const float ABSERR = 1e-6;

Params read_param(string inp_name);
void ez_exit(const char* msg);
void print_banner();
AtomsCount  read_psf(const char* inp_name);
void build_corr(const char* inp_name, AtomsCount atom_num_psf);
void diag_corr();
void fit_coor();
void write_data();

int main(int argc, char* argv[]){
	print_banner();
	const char* psf_name = "step3_pbcsetup.xplor.ext.psf";
	AtomsCount atom_num_psf = read_psf(psf_name);
	const char* dcd_name = "step5_production.dcd";
	build_corr(dcd_name, atom_num_psf);
	diag_corr();
	write_data();
	return 0;
}

Params read_param(string inp_name){
        cout << "ReadParam> Reading parameters from file ";
        cout << " ["<< inp_name << "]" << endl << endl;
        std::ifstream inp_file(inp_name);
        if(!inp_file.is_open()) ez_exit("Config file not found!");
        Params result;
        string each_line;
        while(std::getline(inp_file, each_line)){
                if(!each_line.length() || each_line[0]=='#') continue;
                size_t npos = each_line.find_first_of(' ');
                string token_1 = each_line.substr(0, npos);
                npos = each_line.find_last_of(' ') + 1;
                string token_2 = each_line.substr(npos, each_line.length());
                if(token_1=="job_name"){result.job_name=token_2.c_str();}
                if(token_1=="psf_name"){result.psf_name=token_2.c_str();}
                if(token_1=="dcd_name"){result.dcd_name=token_2.c_str();}
                if(token_1=="n_of_pc"){result.n_of_pc=std::atoi(token_2.c_str());}
        }
        inp_file.close();
        return result;
}

void ez_exit(const char* msg){
        cout << "ERROR> " << msg << endl;
        exit(0);
}

void print_banner(){
	cout << "EZPCA> Principal Components Analysis of Molecular Dynamics Trajectory" << endl;
	cout << "EZPCA> version 1.0" << endl << endl;
}

AtomsCount  read_psf(const char* inp_name){
	cout << "ReadPSF> Reading protein structure from file [";
	cout << inp_name << "] ..." << endl;
	AtomsCount result;
	ifstream inp_file(inp_name);
	string each_line;
	while(getline(inp_file, each_line)){
		if(each_line.find("PSF")!=each_line.npos){
			if(each_line.find("XPLOR")==each_line.npos ||
			   each_line.find("EXT"  )==each_line.npos){
				ez_exit("Currently only EXT and XPLOR formats are supported.");
			}
		}
		if(each_line.find("NATOM")!=each_line.npos) result.n_atoms=atoi(each_line.substr(0,10).c_str());
		if(each_line[0]=='*') cout << "ReadPSF> " << each_line << endl;
		if(each_line[0]!='*' && each_line.size()==118 && each_line.substr(38,2)=="CA"){
			index_CAs.push_back(atoi(each_line.substr(0,10).c_str()));			
		}
	}
	result.n_CAs=index_CAs.size();
	inp_file.close();
	cout << "ReadPSF> After reading psf file, (" << result.n_atoms << ") atoms were recorded";
	cout << ", ("<< result.n_CAs << ") CAs were found." << endl << endl;
	return result;
}

void build_corr(const char* inp_name, AtomsCount atom_num_psf){
	cout << "BuildCorr> Building covariance matrix from trajectory ..." << endl;
	cout << "BuildCorr> Reading trajectory from file [";
	cout << inp_name << "] ..." << endl; 
	std::ifstream inp_file(inp_name, std::ios::binary);
        if(!inp_file.is_open()) ez_exit("Failed to open file.");

        char dcd_head1[100];
        char dcd_head2[16];
        char title[80];
	char test[56];

	inp_file.read(dcd_head1, 100);
        if(strncmp(&dcd_head1[4], "CORD", 4)!=0) ez_exit("Wrong dcd format.");
	cout << "BuildCorr> After reading, the following information were found:" << endl;

	int n_file = *(unsigned int*)(&dcd_head1[8]);
	int n_priv = *(unsigned int*)(&dcd_head1[12]);
	int n_savc = *(unsigned int*)(&dcd_head1[16]);
	int n_step = *(unsigned int*)(&dcd_head1[20]);
	float delta = *(float*)(&dcd_head1[44]);
	int q_cell = *(unsigned int*)(&dcd_head1[48]);
	int c24tag = *(unsigned int*)(&dcd_head1[84]);
	if(c24tag=!24) ez_exit("Currently only NAMD trajectory format is supported.");

	cout << "BuildCorr> NFILE=" << n_file << " NPRIV=" << n_priv;
	cout <<" NSAVC=" << n_savc << " NSTEP=" << n_step << endl;
	cout << "BuildCorr> DELTA="<<delta;
	cout << (q_cell==1?"  PBC cells detected.":"  PBC cells NOT detected.") << endl;
	inp_file.read(title, 80);
	cout << "BuildCorr> "<< title << endl;
	inp_file.read(title, 80);
	cout << "BuildCorr> " << title << endl;
	inp_file.read(dcd_head2, 16);
	int n_atoms = *(int*)(&dcd_head2[8]);
	cout << "BuildCorr> (" << n_atoms << ") atoms found in trajectory file." << endl;
	if(n_atoms!=atom_num_psf.n_atoms) ez_exit("Atoms numbers in PSF and DCD files do not match.");	

/* 
	Till now we have finished reading all header information of dcd file.
	Now read coordinates. If PBC cells detected, use (56) as an offset of each frame.
	Each frame is organized as follows: 
	cell_offset + coor_offset + x_coor + 2*coor_offset + y_coor + 2*coor_offset + z_coor + coor_offset.
	the corresponding size configuration is:
	(14 or 0)   +    (1)      + n_atoms+      (2)      + n_atoms+      (2)      + natom  +    (1).
	Thus, the Cartesian coordinates of the (i-1)-th atom are: (x_offset+i, y_offset+i, z_offset+i).
*/
	unsigned int frame_size;
	unsigned int coorf_size = 4*n_atoms;	//coordinates as floats.
	unsigned int dim3 = 3*atom_num_psf.n_CAs;
	unsigned int x_offset = (q_cell==1) ? 15           : 1;	
	unsigned int y_offset = (q_cell==1) ? n_atoms+17   : n_atoms + 3;
	unsigned int z_offset = (q_cell==1) ? 2*n_atoms+19 : 2*n_atoms+5;
	int coor_offset = 1;

	frame_size = (q_cell==1)?(3*(4*n_atoms+8)+56):(3*(4*n_atoms+8));
	coor_buff.resize(frame_size/4);
	coor_CA_avrg.resize(dim3);

	char frame_buffer[frame_size];
//1st read: Calculate the average coordinates.
	for(int i=0; i<n_file; i++){
		inp_file.read(reinterpret_cast<char*>(coor_buff.data()), frame_size);
		for(int j=0; j<index_CAs.size(); j++){
			coor_CA_avrg[j*3+0] += coor_buff[x_offset + index_CAs[j]];
			coor_CA_avrg[j*3+1] += coor_buff[y_offset + index_CAs[j]];
			coor_CA_avrg[j*3+2] += coor_buff[z_offset + index_CAs[j]];
		}
	}
	for(int j=0; j<dim3; j++)
		coor_CA_avrg[j] /= n_file;
	
//2nd read: Calculate the covariance matrix.
	coor_CA_frame.resize(dim3);
	C.resize(dim3);
	for(int i=0; i<dim3; i++)
		C[i].resize(dim3);

	inp_file.clear();	// seekg will clear eofbit since c++11
	inp_file.seekg(276, ios::beg);	// Rewind dcd file and skip header information.
	for(int s=0; s<n_file; s++){
		inp_file.read(reinterpret_cast<char*>(coor_buff.data()), frame_size);
//Since coor_buff contins both coordinates and garbage data, we need to put coordinates into coor_frame first.
		for(int j=0; j<index_CAs.size(); j++){
			coor_CA_frame[j*3+0] = coor_buff[x_offset + index_CAs[j]];
			coor_CA_frame[j*3+1] = coor_buff[y_offset + index_CAs[j]];
			coor_CA_frame[j*3+2] = coor_buff[z_offset + index_CAs[j]];
		}
//Then we can calculate the accumulated covariance matrix.
		for(int i=0; i<dim3; i++){
			for(int j=0; j<dim3; j++){
				C[i][j] += (coor_CA_frame[i]-coor_CA_avrg[i])*(coor_CA_frame[j]-coor_CA_avrg[j]);
			}
		}
	}
//Finally, calculate the average covariance matrix.
	for(int i=0; i<dim3; i++){
		for(int j=0; j<dim3; j++){
			C[i][j] /= n_file;
		}
	}
	inp_file.close();
}

void diag_corr(){
	cout << "DiagCorr> Diagonalizing correlation matrix ..." << endl;
        unsigned int dim3 = 3*index_CAs.size();
	unsigned int itr_counter = 0;
        float fdim = (float)(dim3*dim3);

// step0. Initialize the eigenvector array.
        E.resize(dim3);
        for(int i=0; i<dim3; i++){
                E[i].resize(dim3);
                E[i][i] = 1.0;
        }
// step1. find the sum of the square of all off-diagonal elements.
        double sum_offd = 0.0;
        for(int i=0; i<dim3; i++){
                for(int j=0; j<dim3; j++){
                        if(i!=j){
                                sum_offd+=C[i][j]*C[i][j];
                        }
                }
        }
        if(sum_offd<=ABSERR) return;
	cout << "DiagCorr> Starting with sum_offd = (" << sum_offd << ")." << endl;
// step2. We compare average for off-diagonal elements to jump over small elements.
        double avg_offd = 0.5*sum_offd/fdim;
        while(sum_offd >= ABSERR){
                itr_counter++;
                cout <<"DiagCorr> Iteration ("<< setw(8) << itr_counter <<"): ";
//                cout << "Current sum_offd :" << setw(16) << setiosflags(ios::fixed) <<setprecision(8) << sum_offd << " | ";
//                cout << "Target sum_offd :" << ABSERR << endl;
// #pragma omp parallel for
                for(int i=0; i<dim3-1; i++){
                        for(int j=i+1; j<dim3; j++){
                                if(C[j][i]*C[j][i]<=avg_offd) continue;
                                sum_offd -= 2.0*C[j][i]*C[j][i];
                                avg_offd = 0.5*sum_offd/fdim;
// step3. Calculate coefficient [c] and [s] for Givens matrix.
                                double beta = (C[j][j]-C[i][i])/(2.0*C[j][i]);
                                double coeff = 0.5*beta/sqrt(1.0+beta*beta);
                                double s = sqrt(std::max(0.5+coeff, 0.0));
                                double c = sqrt(std::max(0.5-coeff, 0.0));
// step4. Update rows [i] and [j] of Hessian matrix (Givens matrix pre-multiplies Hessian).
                        #pragma omp parallel for
                                for(int k=0; k<dim3; k++){
                                        double cs =  c*C[i][k]+s*C[j][k];
                                        double sc = -s*C[i][k]+c*C[j][k];
                                        C[i][k] = cs;
                                        C[j][k] = sc;
                                }
// step5. Givens matrix post-multiplies Hessian. Also calculate eigenvectors. 
                        #pragma omp parallel for
                                for(int k=0; k<dim3; k++){
                                        double cs =  c*C[k][i]+s*C[k][j];
                                        double sc = -s*C[k][i]+c*C[k][j];
                                        C[k][i] = cs;
                                        C[k][j] = sc;
                                        cs =  c*E[k][i]+s*E[k][j];
                                        sc = -s*E[k][i]+c*E[k][j];
                                        E[k][i] = cs;
                                        E[k][j] = sc;
                                }
                        }
                }
		cout << "Current sum_offd :" << setw(16) << setiosflags(ios::fixed) <<setprecision(8) << sum_offd << " | ";
                cout << "Target sum_offd :" << ABSERR << endl;
        }
	cout << "DiagCorr> Done." << endl << endl;
}
void write_data(){
	cout << "WriteData> Writing data ..." << endl;


	cout << "WriteData> Done." << endl;
}
