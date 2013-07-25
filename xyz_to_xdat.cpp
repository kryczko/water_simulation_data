#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

double unitcell (double coord, double lattice)
{
	if (coord > lattice)
	{
		coord = coord - lattice;
	}
return coord;
}


int main()
{
	// main menu for the program, prompts the user for required info to run the program
	string infile;	
	double xlat, ylat, zlat;
	int nframes;

	cout << "Please enter your xyz file: ";
	cin >> infile;
	cout << "Please enter the number of frames: ";
	cin >> nframes;
	cout << "Please enter the lattice constants (x y z) in Angstroms: ";
	cin >> xlat >> ylat >> zlat;
	
	// declare and open input and output files
	ifstream input1, input2;
	ofstream outputfile;

	input1.open(infile.c_str());
	input2.open(infile.c_str());

	outputfile.open("OHxdat.dat");
	
	// go through the inputted xyz file for oxygen
	vector <string> unwanted;
	vector <string> wanted;
	int count = 0;
	string dummy0;
	double dummy1, dummy2, dummy3;	
	vector <char> atomname;
	vector <double> ox, oy, oz, hx, hy, hz;
	
	while (!input1.eof())
	{
		input1 >> dummy0;
		if (dummy0 == "O")
		{
			input1 >> dummy1 >> dummy2 >> dummy3;
			dummy1 = unitcell(dummy1, xlat)/xlat;
			dummy2 = unitcell(dummy2, ylat)/ylat;
			dummy3 = unitcell(dummy3, zlat)/zlat;
			ox.push_back(dummy1);
			oy.push_back(dummy2);
			oz.push_back(dummy3);
		}
	}

	// repeat the previous while loop for hydrogen
	while (!input2.eof())
	{
                input2 >> dummy0;
                if (dummy0 == "H")
                {
                        input2 >> dummy1 >> dummy2 >> dummy3;
                        dummy1 = unitcell(dummy1, xlat)/xlat;
                        dummy2 = unitcell(dummy2, ylat)/ylat;
                        dummy3 = unitcell(dummy3, zlat)/zlat;
                	hx.push_back(dummy1);
			hy.push_back(dummy2);
			hz.push_back(dummy3);
		}
	}
	
	int nooa = ox.size()/nframes;
	int noha = hx.size()/nframes;

	for (int i = 0; i < nframes; i ++)
	{
		for (int j = 0; j < nooa; j ++)
		{
			outputfile << ox[j + i*nooa] << "  " << oy[j + i*nooa] << "  " << oz[j + i*nooa] << endl;
		}
		for (int k = 0; k < noha; k ++)
		{
			outputfile << hx[k + i*noha] << "  " << hy[k + i*noha] << "  " << hz[k + i*noha] << endl;
		}
	}
	input1.close();
	input2.close();
	outputfile.close();
	
	return 0;
}
