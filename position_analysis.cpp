#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

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


int main()
{
// file streams
ifstream inputfile;
ofstream oo_outputfile, oh_outputfile, angle_outputfile;

char oodistance, ohdistance, angle; 
string infile;
int timesteps, number_of_atoms;
double lattice;

// main menu for the program
cout << "\n\nWelcome to the bulk-water data analysis program!\n\n";
cout << "If you would like to calculate the average O-O distance between the closest adjacent atoms with respect to timesteps enter \"y\", if not enter \"n\": ";
cin >> oodistance;
cout << "If you would like to plot the average O-H distance with respect to timesteps enter \"y\", if not enter \"n\": ";
cin >> ohdistance;
cout << "If you would like to plot the average H-O-H angle with respect to timesteps enter \"y\", if not enter \"n\": ";
cin >> angle;
cout << "Please enter the filename of your file: ";
cin >> infile;
cout << "Please enter the number of atoms: ";
cin >> number_of_atoms;
cout << "Please enter the number of timesteps: ";
cin >> timesteps;
cout << "Please enter the lattice constant for your periodic cube: ";
cin >> lattice;
// end of main menu

//read the inputfile
inputfile.open(infile.c_str());

double a,b,c,x[number_of_atoms*timesteps], y[number_of_atoms*timesteps], z[number_of_atoms*timesteps];
int counter = 0;

while (!inputfile.eof())
{
	//inputfile >> x[counter] >> y[counter] >> z[counter];	
	inputfile >> a >> b >> c;
	counter ++;
}
cout << counter << endl;
// done reading inputfile

/*
//decision for oodistance
if (oodistance == 'y')
{
	oo_outputfile.open("oo_avg_distance.dat");
	double ox[number_of_atoms/3], oy[number_of_atoms/3], oz[number_of_atoms/3];
	double dx, dy, dz;
	double difference[number_of_atoms/3 - 1], last_difference[number_of_atoms/3 - 1];	

	for (int i = 0; i < timesteps; i ++)
	{
		for (int j = 0; j < number_of_atoms/3; j ++)
		{
			ox[j] = lattice*x[j + i*number_of_atoms];
			oy[j] = lattice*y[j + i*number_of_atoms];
			oz[j] = lattice*z[j + i*number_of_atoms];
		}
		
		for (int k = 0; k < number_of_atoms/3 - 1; k ++)
		{
			for (int n = 1; n < number_of_atoms/3; n ++)
			{
				dx = ox[n] - ox[k];
				dy = oy[n] - oy[k];
				dz = oz[n] - oz[k];
				
				dx -= lattice*pbc_round(dx/lattice);
				dy -= lattice*pbc_round(dy/lattice);
				dz -= lattice*pbc_round(dz/lattice);
			
				double distance = sqrt( dx*dx + dy*dy + dz*dz );
				difference[n - 1] = distance;
			}
		
			double lowest = difference[0];
			for ( int w = 1; w < number_of_atoms/3 - 1; w ++)
			{
				if ( difference[w] < lowest )
				{
					lowest = difference[w];
				}
			}
			
			last_difference[k] = lowest;
			oo_outputfile << lowest << endl;	
		}					
	}	
	cout << "\n\nYour O-O average distance data with respect to timesteps has been placed in \"oo_avg_distance.dat\" and can now be easily plotted with gnuplot.\n\n"; 
}
else
{
	cout << "\n\nYou either entered the wrong key or you do not want to plot the average O-O distance.\n\n";
}

//decision for ohdistance
if (ohdistance == 'y')
{
	oh_outputfile.open("oh_avg_distance.dat");


        cout << "\n\nYour O-H average distance data with respect to timesteps has been placed in \"oh_avg_distance.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
        cout << "\n\nYou either entered the wrong key or you do not want to plot the average O-H distance.\n\n";
}

//decision for angle
if (angle == 'y')
{
        angle_outputfile.open("avg_angle.dat");


        cout << "\n\nYour H-O-H average angle data with respect to timesteps has been placed in \"avg_angle.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
        cout << "\n\nYou either entered the wrong key or you do not want to plot the average H-O-H angle.\n\n";
}





*/
inputfile.close();



return 0;
}
