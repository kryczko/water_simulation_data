/*The original intentions of this program is that the user enters an xyz file and lets the program know what 2 atom types
  they would like to compute frequency for those 2 atoms.*/

#include <iostream> // input-output availability
#include <fstream> // input-output file availability
#include <string> // for using strings
#include <vector> // for using vectors (more convienient than arrays)
#include <cmath> // for mathematical functions like sin, sqrt, etc.
#include <iomanip> // for printing out fixed amounts of digits
#include <fftw3.h> // For fourier transform

using namespace std;

//####################################
//function to deal with the periodic boundary conditions
int pbc_round(double input)
{
	int i  = input;

	if (abs(input - i) >= 0.5)
	{
		if (input > 0) {i += 1;}
		if (input < 0) {i -= 1;}
	}
return i;
}
//####################################

int main()
{
	//###########################################################################
	// Main menu, get info from user
	string infile, atom1, atom2;
	double xlat, ylat, zlat, timestep;
	int noa1, noa2;
	
	cout << "XYZ file\n==> ";
	cin >> infile;
	cout << "Atom names (e.g. O H)\n==> ";
	cin >> atom1 >> atom2;
	cout << "Number of each atom\n==> ";
	cin >> noa1 >> noa2;
	cout << "Lattice constants (x y z)\n==> ";
	cin >> xlat >> ylat >> zlat;
	cout << "Timestep (fs)\n==> ";
	cin >> timestep;
	cout << "\n\n########## PROGRAM RUNNING, PLEASE WAIT PATIENTLY ##########\n\n";
	//#############################################################################

	//############################################################################
	// Go through inputfile
	ifstream input;
	input.open(infile.c_str());

	vector <double> a1x, a2x, a1y, a2y, a1z, a2z;
	string dummy;
	double x, y, z;
	while (! input.eof())
	{
		input >> dummy;
		if (dummy == atom1)
		{
			input >> x >> y >> z;
			a1x.push_back(x);	
			a1y.push_back(y);
			a1z.push_back(z);
		}
		if (dummy == atom2)
                {
                        input >> x >> y >> z;
                        a2x.push_back(x);
                        a2y.push_back(y);
                        a2z.push_back(z);
                }
 	}
	cout << "########## DONE READING INPUTFILE ##########\n\n";
	//###############################################################################
	
	//###############################################################################
	// Go through the data obtained from the inputfile

	//firstly, lets create the neighbor list array; since this is general lets say there are maybe 6 neighbors at most. 
        int nframes = a1x.size() / noa1;
	
	//velocity array for every frame, atom2 and in the x,y, and z direction
	double vel[nframes - 1][noa2][3];
	ofstream output;
	output.open("freq.dat");

	for (int i = 0; i < nframes - 1; i ++)
	{
		for (int j = 0; j < noa2; j ++)
		{
			double dx = a2x[j + (i + 1)*noa2] - a2x[j + (i)*noa2];
                        double dy = a2y[j + (i + 1)*noa2] - a2y[j + (i)*noa2];
                        double dz = a2z[j + (i + 1)*noa2] - a2z[j + (i)*noa2];

			dx -= xlat*pbc_round(dx/xlat); 
			dy -= ylat*pbc_round(dy/ylat);
			dz -= zlat*pbc_round(dz/zlat);

			vel[i][j][0] = dx/timestep;
			vel[i][j][1] = dy/timestep;
			vel[i][j][2] = dz/timestep;
		}
	}

	cout << "########## COMPUTED VELOCITES FROM DATAFILE ##########\n\n";

	double *Z;
	Z = new double [nframes - 1];

	for (int m = 0; m < nframes - 1; m ++)
	{
		for (int n = 0; n < nframes - 1 - m - 1; n ++)
		{
			for (int i = 0; i < noa2; i ++)
			{
				for (int j = 0; j < 3; j ++)
				{
					Z[m] += vel[n + m][i][j]*vel[n][i][j];
				}				
			}
		}
	}
	for (int m = 0; m < nframes - 1; m ++)
	{
		Z[m] /= (nframes - m);
	}
	for (int m = 1; m < nframes - 1; m ++)
	{
		Z[m] /= Z[0];
	}
	cout << "########## COMPUTED THE VELOCITY AUTOCORRELATION FUNCTION ##########\n\n";
	
	int N = nframes - 1;
	double sigma = N / 2.5;
	
	for (int i = 0; i < N; i ++)
	{
		Z[i] *= exp(-i*i / (2*sigma*sigma)) / (sigma * sqrt(2*3.14159265359));
	}
	for (int m = 1; m < N; m ++)
	{
		Z[m] /= Z[0];
	}	
	
	//##############################################################
	// fourier transform
	fftw_complex *out;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2)+1);
	fftw_plan p;
	p = fftw_plan_dft_r2c_1d(N,Z,out,FFTW_ESTIMATE);
	fftw_execute(p);
	
	cout << "########### COMPUTED FOURIER TRANSFORMS ##########\n\n";

	for (int i = 0; i < (N/2)+1; i ++)
	{
		output << i <<"\t"<< out[i][0] << "\t" << out[i][1] << endl;	
	}
	cout << "######### OUTPUTTED DATA TO freq.dat ##########\n\n";
	
	input.close();
	output.close();
	return 0;
}
