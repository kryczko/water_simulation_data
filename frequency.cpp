/*The original intentions of this program is that the user enters an xyz file and lets the program know what 2 atom types
  they would like to compute frequency for those 2 atoms.*/

#include <iostream> // input-output availability
#include <fstream> // input-output file availability
#include <string> // for using strings
#include <vector> // for using vectors (more convienient than arrays)
#include <cmath> // for mathematical functions like sin, sqrt, etc.
#include <iomanip> // for printing out fixed amounts of digits

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
	double xlat, ylat, zlat, equildist, maxdist;
	int noa1, noa2;
	
	cout << "XYZ file\n==> ";
	cin >> infile;
	cout << "Atom names (e.g. O H)\n==> ";
	cin >> atom1 >> atom2;
	cout << "Number of each atom\n==> ";
	cin >> noa1 >> noa2;
	cout << "Lattice constants (x y z)\n==> ";
	cin >> xlat >> ylat >> zlat;
	cout << "Maximum bond length and equilibrium bondlength (in Angstrom)\n==> ";
	cin >> maxdist >> equildist;
	cout << "\n\n########## PROGRAM RUNNING, PLEASE WAIT PATIENTLY ##########\n\n"
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
        int nframes = alx.size() / noa1;

	int neighbors[noa1*nframes][6];	
	for (int i = 0; i < noa1*nframes; i ++)
	{
		for (int j = 0; j < 6; j ++)
		{
			neighbors[i][j] = -1
		}
	}	

	double ndist[noa1*nframes][6];
	for (int i = 0; i < noa1*nframes; i ++)
        {
                for (int j = 0; j < 6; j ++)
                {
                        ndist[i][j] = -1
                }
        }


	for (int i = 0; i < nframes; i ++)
	{
		for (int j = 0; j < noa1; j ++)
		{
			int ncount = 0;
			for (int k = 0; k < noa2; k ++)
			{
				double dx = a1x[j + i*noa1] - a2x[k + i*noa2];
                                double dy = a1y[j + i*noa1] - a2y[k + i*noa2];
                                double dz = a1z[j + i*noa1] - a2z[k + i*noa2];
				
				dx -= xlat*pbc_round(dx/xlat);
        	                dy -= ylat*pbc_round(dy/ylat);
	                        dz -= zlat*pbc_round(dz/zlat);
				
				double dist = sqrt ( dx*dx + dy*dy + dz*dz );
				//if dist less than the max bond length, then it must be a neighbor!
				if (dist <= maxdist)
				{
					ndist[k +i*noa1][ncount] = dist;
					neighbors[k + i*noa1][ncount] = k;
					ncount ++;
				}
			}
		}
	}	 
				
	for (int i = 0; i < nframes
