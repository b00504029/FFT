// End-bead velocity autocorrelation function w.r.t grafting points 
// Read the input only when the frames need to be calculated
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <complex>
#include <fftw3.h>
#pragma comment(lib, "libfftw3-3.lib")
using namespace std;

class ATOM
{
    public:
    // Parameters for sorting particles
	string str, name;                         				// Store the instruction, the file name
    int tot_part, tot_frame, tot_type;        				// Total # of particles, frames, and atom types
    int n_chain;                              				// # of total chains in one NOHMs particle
    int n_bead;                               				// # of beads per chain
    int n_NOHMs;                               				// # of NOHMs in the system

	// Useful vectors and variables
	vector< vector< double > > vrel;    					// Relative velocities between end beads: [frames][direction]
	
    // Parameters for VACF
    double dt;                                				// Timestep (ps) between frames
	int ti, tf;												// Initial, final frames in the trajectory
	int taui, tauf, taus;                					// Minimum, final, stride of, current span of frames
	vector< vector < double > > VACF;						// velocity autocorrelation function: [time][direction]
    void VACF_calc();										// velocity autocorrelation function calculation
	void ens_avg_calc(double []);							// Ensemble average of velocity w.r.t grafting-point

    // Initialization
    void init();

	// Input and output
	ifstream intrj;											// Load the trajectory
	ofstream outf;											// Write the VACF result

    // Skipping frames unwanted until ti
    void skip_lines(int ,ifstream* );
};

int main(int argc, char *argv[]){
    // Construct a system
    ATOM atom;

    // File name to read
	if (argc < 2) {
		//cout << "ERROR: no trajectory filename !" << endl;
		return 0;
	}
	else {
		atom.name = argv[1];
	}

    // Timestep in ps
	/*
    cout<<"Enter the timestep in ps between frames: ";
    cin>>atom.dt;
	*/
	atom.dt = 0.012;

    // Lag times t from taui to tauf by taus
	/*
    cout<<"Enter the initial lag time, the final lag time,and a stride of lag times: ";
    cin>>atom.taui>>atom.tauf>>atom.taus;
	*/
	atom.taui = 0;
	atom.tauf = 90000;
	atom.taus = 1;
	
	// For analysis
	atom.ti = 0;
	atom.tf = 100000;

    // Calculate the wall-clock time
    auto wcts = std::chrono::system_clock::now();

    // Initialization
    atom.init();

    // End-bead velocity autocorrelation function (VACF) calculation
    atom.VACF_calc();

    // End of calculation
    std::chrono::duration< double > wctduration = (std::chrono::system_clock::now() - wcts);
    std::cout<< "Finished in "<< wctduration.count() << " seconds [Wall Clock]" <<std::endl;

    return 0;
}

void ATOM::init()
{
    // Assign the initial values to all parameters
    // Atom types: 1=B 2=C 3=S
	tot_frame = tf-ti+1;
	tot_type = 3;

	// System parameters
	n_chain = 25;
	n_bead = 15;
	n_NOHMs = 256;
    tot_part = n_NOHMs*n_chain*2;

	// Allocate VACF array
	VACF.resize(tauf-taui+1);
	for (int i; i < VACF.size(); ++i) VACF[i].resize(3,0);

	// Construct the array vrel
	vrel.resize(tot_frame);
	for (int i = 0; i < vrel.size(); ++i) vrel[i].resize(3,0);
	
}

void ATOM::VACF_calc()
{
	// Useful variables for VACF calculation
	int temp_id;										// Temporary ID
	int temp_mol;										// Temporary molecule-tag
    int temp_type;										// Temporary atom type
	int ichain;											// Temporary chain index
	int ncycle = 20000;									// # of data per cycle
	double temp_x;										// Temporary coordinates
	double temp_v[2][3];								// Temporary velocities
	double v_ens[3];									// Ensemble average of velocity w.r.t grafting-point
	ofstream log[3];									// Output each end-to-end vector component
	
	vector< vector< vector< double > > > temp_vrel;		// Temporary relative velocities
	temp_vrel.resize(ncycle);
	for (int i = 0; i < temp_vrel.size(); ++i) {
		temp_vrel[i].resize(n_NOHMs*n_chain);
		for (int j = 0; j < temp_vrel[i].size(); ++j) {
			temp_vrel[i][j].resize(3,0);
		}
	}
	
	// Ensemble average calculation
	//ens_avg_calc(v_ens);
	
	// Skip the lines unwanted until the first frame
	for (int i = 0; i < ti; ++i){
		skip_lines(9+tot_part,& intrj);
	}
	
	// Data generation
	intrj.open(name,ios::in);
    for (int i = ti; i < tf; ++i){
		// Read one frame at a time
		skip_lines(9,&intrj);
		for (int j = 0; j < tot_part; ++j){
			ichain = floor(j/2);
			// One end
			if ( j%2 == 0 ) {
				intrj >> temp_id;
				intrj >> temp_mol;
				intrj >> temp_type;
				intrj >> temp_x >> temp_x >> temp_x;
				for (int k = 0; k < 3; ++k) intrj >> temp_v[j%2][k];
			}
			// The other end
			if ( j%2 == 1 ) {
				intrj >> temp_id;
				intrj >> temp_mol;
				intrj >> temp_type;
				intrj >> temp_x >> temp_x >> temp_x;
				for (int k = 0; k < 3; ++k) {
					intrj >> temp_v[j%1][k];
					temp_vrel[i%ncycle][ichain][k] = temp_v[1][k]-temp_v[0][k];
				}
			}
			getline(intrj,str);
		}
		
		if ( i%ncycle == ncycle-1 ) {
			for (int j = 0; j < n_NOHMs*n_chain; ++j) {
				log[0].open("end-to-end_vec_x_"+to_string(j+1)+".txt",ios::out|ios::app);
				log[1].open("end-to-end_vec_y_"+to_string(j+1)+".txt",ios::out|ios::app);
				log[2].open("end-to-end_vec_z_"+to_string(j+1)+".txt",ios::out|ios::app);
					
				for (int icycle = 0; icycle < ncycle; ++icycle) 
					for (int k = 0; k < 3; ++k) 
						log[k] << temp_vrel[icycle][j][k] << endl;
					
				log[0].close();
				log[1].close();
				log[2].close();
			}
		}
		else if ( i == tf-1 ) {
			for (int j = 0; j < n_NOHMs*n_chain; ++j){
				log[0].open("end-to-end_vec_x_"+to_string(j+1)+".txt",ios::out|ios::app);
				log[1].open("end-to-end_vec_y_"+to_string(j+1)+".txt",ios::out|ios::app);
				log[2].open("end-to-end_vec_z_"+to_string(j+1)+".txt",ios::out|ios::app);
				
				for (int k = 0; k < 3; ++k) log[k] << temp_vrel[0][j][k] << endl;
				
				log[0].close();
				log[1].close();
				log[2].close();
			}
		}
		
		if (i%1000 == 0) cout << i << endl;
	}	
	intrj.close();
	

/*
	// Read the required data and perform VACF calculation using FFT
	log[0].open("end-to-end_vec_x.txt",ios::in);
	log[1].open("end-to-end_vec_y.txt",ios::in);
	log[2].open("end-to-end_vec_z.txt",ios::in);
	
	for (int i = 0; i < vrel.size(); ++i) {
		cout << i << endl;
		for (int j = 0; j < vrel[i].size(); ++j) {
			log[j] >> vrel[i][j];
			getline(log[j],str);
		}
	}
			
	log[0].close();
	log[1].close();
	log[2].close();
	
	// VACF calculation using FFT
	int fft_leng = 2*tot_frame;
	fftw_complex inp[fft_leng],out[fft_leng];

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < fft_leng; ++j) {
			inp[j][0] = inp[j][1] = out[j][0] = out[j][1] = 0;
		}
		
		for (int j = 0; j < vrel.size(); ++j)
			inp[j][0] = vrel[j][i];

		// Forward FFT
		fftw_plan planf = fftw_plan_dft_1d(fft_leng, inp, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(planf);
		fftw_destroy_plan(planf);

		for (int j = 0; j < fft_leng; ++j) {
			out[j][0] = pow(out[j][0],2)+pow(out[j][1],2);
			out[j][1] = 0;
		}

		// Backward FFT
		fftw_plan planb = fftw_plan_dft_1d(fft_leng, out, inp, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(planb);
		fftw_destroy_plan(planb);

		for (int j = 0; j < VACF.size(); ++j) 
			VACF[j][i] = inp[j][0]/real(fft_leng*(tot_frame-i)*n_NOHMs*n_chain);

	}

	// Output the VACF result
    outf.open("vacf_wrt_gp.txt",ios::out);
    outf<<"t"<<" "<<"VACF_x"<<" "<<"VACF_y"<<" "<<"VACF_z"<<endl;
	for (int i = 0; i < VACF.size(); ++i) {
		outf << i*dt <<" ";
		for (int j = 0; j < VACF[i].size(); ++j) {
			outf << VACF[i][j] <<" ";
		}
		outf << endl;
	}
    outf.close();
*/
}

/*
void ATOM::ens_avg_calc(double v_ens[])
{
	// Necessary variables
	int temp_id, temp_mol, temp_type;				// Temporary attributes
	double temp_x;									// Temporary coordinates
	double vx,vy,vz;									// Relative velocity
	double v_size;									// Speed

	// Read the required data and perform VACF calculation
	intrj.open(name,ios::in);
	//outf.open("end-bead_velocity.txt",ios::out);

	for(int i = 0; i < 3; ++i) v_ens[i] = 0;

	int n = 0; // Counter
    do {
		skip_lines(9,& intrj);
		for (int i = 0; i < tot_part; ++i){
			intrj >> temp_id;
			intrj >> temp_mol;
			intrj >> temp_type;
			intrj >> temp_x >> temp_x >> temp_x;
			it_type = find(itype.begin(), itype.end(), temp_type);
			it_id = find(id[ temp_type-1 ].begin(), id[ temp_type-1 ].end(), temp_id);

			if (it_type != itype.end()) {
				for (int j = 0; j < 3; ++j) {
					// Save the coordinates of the current frame
					intrj >> v[0][ temp_type-1 ][ distance( id[ temp_type-1 ].begin(), it_id)][j];

				}
			}
			getline(intrj,str);
		}
		
		// Ensemble average calculation
		for (int i = 0; i < n_NOHMs*n_chain; ++i){
			vx = v[0][0][i][0] - v[0][2][i][0];
			vy = v[0][0][i][1] - v[0][2][i][1];
			vz = v[0][0][i][2] - v[0][2][i][2];
			v_size = sqrt(vx*vx+vy*vy+vz*vz);
			//outf<<" "<< vx <<" "<< vy <<" "<< vz <<" "<< v_size <<endl;
			
			v_ens[0] += vx;
			v_ens[1] += vy;
			v_ens[2] += vz;
		}
		++n;
	}while( !intrj.eof() );

	cout << n << endl;
	for (int j = 0; j < 3; ++j) {
		v_ens[j] /= (n*n_NOHMs*n_chain);
	}

	intrj.close();
	outf.close();
}
*/

void ATOM::skip_lines(int num_lines, ifstream* data)
{
	int i;
    for (i = 0; i < num_lines; i++){
        getline(*data,str);
    }
}
