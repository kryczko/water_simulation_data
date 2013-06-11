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

double min_distance(double array[], int length)
{
double min = array[0];
for (int i = 1; i < length; i ++)
{
	if (array[i] < min && array[i] != 0)
	{
		min = array[i];
	}
}
return min;
}


int main()
{
// file streams
ifstream inputfile;
ofstream oo_outputfile, oh_outputfile, angle_outputfile, hbonds_outputfile;

char oodistance, ohdistance, angle, hbonds; 
string infile;
int timesteps, number_of_atoms;
double lattice;

// main menu for the program
cout << "\n\nWelcome to the bulk-water data analysis program!\n\n";
cout << "If you would like to create a histogram of probability with respect to neighboring O-O distance enter \"y\", if not enter \"n\": ";
cin >> oodistance;
cout << "If you would like to create a histogram of probability with respect to the O-H distance \"y\", if not enter \"n\": ";
cin >> ohdistance;
cout << "If you would like to create a histogram of probability with respect to the H-O-H angle enter \"y\", if not enter \"n\": ";
cin >> angle;
cout << "If you would like to create a histogram of the number of H-bonds with respect to the z-axis enter \"y\", if not enter \"n\": ";
cin >> hbonds;
cout << "Please enter the filename of your file: ";
cin >> infile;
cout << "Please enter the number of atoms: ";
cin >> number_of_atoms;
cout << "Please enter the number of timesteps: ";
cin >> timesteps;
cout << "Please enter the lattice constant for your periodic cube: ";
cin >> lattice;
cout << "Program running...please wait a moment.\n\n";
// end of main menu*/

//read the inputfile
inputfile.open(infile.c_str());

double x[number_of_atoms*timesteps], y[number_of_atoms*timesteps], z[number_of_atoms*timesteps];
int counter = 0;

while (!inputfile.eof())
{
	inputfile >> x[counter] >> y[counter] >> z[counter];	
	counter ++;
	
}

// done reading inputfile


//decision for oodistance
if (oodistance == 'y')
{
	oo_outputfile.open("oo_histogram.dat");
	double ox[number_of_atoms/3], oy[number_of_atoms/3], oz[number_of_atoms/3];
	double dx, dy, dz, lowest;
	int bin_number;
	double last_difference[(number_of_atoms/3-1)*timesteps], bin[100] = {};	

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
			double difference[number_of_atoms/3 - (k+1)];
			for (int n = k + 1; n < number_of_atoms/3; n ++)
			{
				dx = ox[n] - ox[k];
				dy = oy[n] - oy[k];
				dz = oz[n] - oz[k];
				
				dx -= lattice*pbc_round(dx/lattice);
				dy -= lattice*pbc_round(dy/lattice);
				dz -= lattice*pbc_round(dz/lattice);
			
				double distance = sqrt( dx*dx + dy*dy + dz*dz );
				difference[n - (k + 1)] = distance;
		
			}
		
 		last_difference[k + i*number_of_atoms/3] = min_distance(difference, number_of_atoms/3 - (k + 1));
		}

	
					
	}
	int n = 0;
	double sum(0);
	for ( int i = 0; i < (number_of_atoms/3-1)*timesteps; i ++)
	{
		if ((last_difference[i] != 0 || last_difference[i] != 0.0) && (last_difference[i] < 7.0))
		{
			bin_number = last_difference[i]*10;
			bin[bin_number] += 1;
			sum += last_difference[i];
			n ++;
		}
	}
	oo_outputfile << "# The average is " << sum/n << "\n\n";
	for(int i = 0; i < 100; i ++)
	{
		oo_outputfile << i << "\t" << bin[i]/126000 << endl;
		oo_outputfile << i + 1 << "\t" << bin[i]/126000 << endl;
	}
	cout << "\n\nYour O-O histogram data has been placed in \"oo_histogram.dat\" and can now be easily plotted with gnuplot.\n\n"; 
}
else
{
	cout << "\n\nYou either entered the wrong key or you do not want to plot the O-O distance.\n\n";
}

//decision for ohdistance
if (ohdistance == 'y')
{
	oh_outputfile.open("oh_histogram.dat");
	
	double ox[number_of_atoms/3], oy[number_of_atoms/3], oz[number_of_atoms/3], hx[2*number_of_atoms/3], hy[2*number_of_atoms/3], hz[2*number_of_atoms/3];
	double dx1, dx2, dy1, dy2, dz1, dz2, distance[2*number_of_atoms/3], final_distance[2*number_of_atoms/3*timesteps], bin[200]={}  ;	
	int bin_number;

	for (int i = 0; i < timesteps; i ++)
	{
		 for (int j = 0; j < number_of_atoms/3; j ++)
                {
                        ox[j] = lattice*x[j + i*number_of_atoms];
                        oy[j] = lattice*y[j + i*number_of_atoms];
                        oz[j] = lattice*z[j + i*number_of_atoms];
		}
		for (int k = 0; k < 2*number_of_atoms/3; k ++)
		{
			hx[k] = lattice*x[number_of_atoms/3 + k + i*number_of_atoms];
			hy[k] = lattice*y[number_of_atoms/3 + k + i*number_of_atoms];
			hz[k] = lattice*z[number_of_atoms/3 + k + i*number_of_atoms];
		}
		for (int n = 0; n < number_of_atoms/3; n ++)
		{
			dx1 = ox[n] - hx[2*n];
			dx2 = ox[n] - hx[2*n + 1];
			dy1 = oy[n] - hy[2*n];
			dy2 = oy[n] - hy[2*n + 1];
			dz1 = oz[n] - hz[2*n];
			dz2 = oz[n] - hz[2*n + 1];
			
			dx1 -= lattice*pbc_round(dx1/lattice);
                        dy1 -= lattice*pbc_round(dy1/lattice);
                        dz1 -= lattice*pbc_round(dz1/lattice);
	
			dx2 -= lattice*pbc_round(dx2/lattice);
                        dy2 -= lattice*pbc_round(dy2/lattice);
                        dz2 -= lattice*pbc_round(dz2/lattice);

			distance[2*n] = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
			distance[2*n + 1] = sqrt ( dx2*dx2 + dy2*dy2 + dz2*dz2 );
		}
		for (int g = 0; g < 2*number_of_atoms/3; g ++)
		{	
			final_distance[g + i*2*number_of_atoms/3] = distance[g];
		}
			
	}
	double sum(0);
	int n(0);	
	for (int i = 0; i < 2*number_of_atoms/3*timesteps; i ++)
	{
	
		bin_number = final_distance[i]*100;
		bin[bin_number] += 1;
		sum += final_distance[i];
		n ++;
	}	
	oh_outputfile << "# The average is " << sum/n << "\n\n";
	for (int j = 0; j < 200; j ++)
	{
		oh_outputfile << j << "\t" << bin[j]/256000. << endl;
		oh_outputfile << j + 1 << "\t" << bin[j]/256000. << endl;
	}
        cout << "\n\nYour O-H distance histogram data has been placed in \"oh_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
        cout << "\n\nYou either entered the wrong key or you do not want to plot the O-H distance.\n\n";
}

//decision for angle
if (angle == 'y')
{
        angle_outputfile.open("angle_histogram.dat");

	double ox[number_of_atoms/3], oy[number_of_atoms/3], oz[number_of_atoms/3], hx[2*number_of_atoms/3], hy[2*number_of_atoms/3], hz[2*number_of_atoms/3];
	double dx1, dx2, dy1, dy2, dz1, dz2, dot[number_of_atoms/3], distance[2*number_of_atoms/3], angle[number_of_atoms/3*timesteps], bin[20000] = {};
	int bin_number;

	for (int i = 0; i < timesteps; i ++)
	{
		 for (int j = 0; j < number_of_atoms/3; j ++)
                {
                        ox[j] = lattice*x[j + i*number_of_atoms];
                        oy[j] = lattice*y[j + i*number_of_atoms];
                        oz[j] = lattice*z[j + i*number_of_atoms];
                }
                for (int k = 0; k < 2*number_of_atoms/3; k ++)
                {
                        hx[k] = lattice*x[number_of_atoms/3 + k + i*number_of_atoms];
                        hy[k] = lattice*y[number_of_atoms/3 + k + i*number_of_atoms];
                        hz[k] = lattice*z[number_of_atoms/3 + k + i*number_of_atoms];
                }
                for (int n = 0; n < number_of_atoms/3; n ++)
                {
                        dx1 = ox[n] - hx[2*n];
                        dx2 = ox[n] - hx[2*n + 1];
                        dy1 = oy[n] - hy[2*n];
                        dy2 = oy[n] - hy[2*n + 1];
                        dz1 = oz[n] - hz[2*n];
                        dz2 = oz[n] - hz[2*n + 1];

                        dx1 -= lattice*pbc_round(dx1/lattice);
                        dy1 -= lattice*pbc_round(dy1/lattice);
                        dz1 -= lattice*pbc_round(dz1/lattice);

                        dx2 -= lattice*pbc_round(dx2/lattice);
                        dy2 -= lattice*pbc_round(dy2/lattice);
                        dz2 -= lattice*pbc_round(dz2/lattice);
			
			dot[n] = dx1*dx2 + dy1*dy2 + dz1*dz2;

                        distance[2*n] = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
                        distance[2*n + 1] = sqrt ( dx2*dx2 + dy2*dy2 + dz2*dz2 );
                }
                for (int g = 0; g < number_of_atoms/3; g ++)
                {
                        angle[g + i*number_of_atoms/3] = acos( dot[g]/(distance[2*g]*distance[2*g + 1]) ) * 57.2957795;
                }
	
	}
	double sum(0);
	int n(0);
	for (int i = 0; i < number_of_atoms/3*timesteps; i ++)
	{
		bin_number = angle[i]*100;
		bin[bin_number] += 1;
		sum += angle[i];
		n ++;
	}
	angle_outputfile << "# The average is " << sum/n << "\n\n";
	for (int j = 0; j < 20000; j ++)
	{
		angle_outputfile << j << "\t" << bin[j]/128000. << endl;
		angle_outputfile << (j + 1) << "\t" << bin[j]/128000. << endl;
	}

        cout << "\n\nYour H-O-H angle histogram data has been placed in \"angle_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
        cout << "\n\nYou either entered the wrong key or you do not want to plot the H-O-H angle.\n\n";
}

if (hbonds == 'y')
{
	hbonds_outputfile.open("hbonds_histogram.dat");
	
	double ox[number_of_atoms/3], oy[number_of_atoms/3], oz[number_of_atoms/3], hx[2*number_of_atoms/3], hy[2*number_of_atoms/3], hz[2*number_of_atoms/3];
	double ohdifference[2*number_of_atoms*number_of_atoms/9], oodifference[number_of_atoms*number_of_atoms/9], nhlist[number_of_atoms/3][5], nolist[number_of_atoms/3][5];
	double dxo, dyo, dzo, dxh1, dxh2, dyh1, dyh2, dzh1, dzh2;
		

	for (int i = 0; i < timesteps; i ++)
        {
                for (int j = 0; j < number_of_atoms/3; j ++)
                {
                        ox[j] = lattice*x[j + i*number_of_atoms];
                        oy[j] = lattice*y[j + i*number_of_atoms];
                        oz[j] = lattice*z[j + i*number_of_atoms];
                }

		 for (int k = 0; k < 2*number_of_atoms/3; k ++)
                {
                        hx[k] = lattice*x[number_of_atoms/3 + k + i*number_of_atoms];
                        hy[k] = lattice*y[number_of_atoms/3 + k + i*number_of_atoms];
                        hz[k] = lattice*z[number_of_atoms/3 + k + i*number_of_atoms];
                }

		for (int u = 0; u < number_of_atoms/3; u ++)
		{
			for (int g = 0; g < 2*number_of_atoms/3; g ++)
			{
		
				dxh1 = ox[u] - hx[2*g];
				dxh2 = ox[u] - hx[2*g + 1];
				dyh1 = oy[u] - hy[2*g];
				dyh2 = oy[u] - hy[2*g + 1];
				dzh1 = oz[u] - hz[2*g];
				dzh2 = oz[u] - hz[2*g + 1];

				dxh1 -= lattice*pbc_round(dxh1/lattice);
				dxh2 -= lattice*pbc_round(dxh2/lattice);
				dyh1 -= lattice*pbc_round(dyh1/lattice);
				dyh2 -= lattice*pbc_round(dyh2/lattice);
				dzh1 -= lattice*pbc_round(dzh1/lattice);
				dzh2 -= lattice*pbc_round(dzh2/lattice);
				
				ohdifference[2*g + u*2*number_of_atoms/3] = sqrt ( dxh1*dxh1 + dyh1*dyh1 + dzh1*dzh1 );
				ohdifference[2*g + 1 + u*2*number_of_atoms/3] = sqrt ( dxh2*dxh2 + dyh2*dyh2 + dzh2*dzh2 );
				
			}
		}
		
		for (int q = 0; q < number_of_atoms/3; q ++)
		{
			int c = 0, d = 1;
			nhlist[q][c] = q;			

			for (int w = 0; w < 2*number_of_atoms/3; w ++)
			{
				if (ohdifference[w + q*2*number_of_atoms/3] < 1.2)
				{	
					nhlist[q][d] = w + q*2*number_of_atoms/3;
					d ++;
				}
			}
		}
		
		for ( int p = 0; p < number_of_atoms/3; p ++)
		{
			for (int m = 0; m < number_of_atoms/3; m ++)
			{
				dxo = ox[p] - ox[m];
				dyo = ox[p] - ox[m];
				dzo = ox[p] - ox[m];					
				
				dxo -= lattice*pbc_round(dxo/lattice);
                                dyo -= lattice*pbc_round(dyo/lattice);
                                dzo -= lattice*pbc_round(dzo/lattice);
				
				oodifference[m + p*number_of_atoms/3] = sqrt ( dxo*dxo + dyo*dyo + dzo*dzo );
			}
		}

		for (int p = 0; p < number_of_atoms/3 ; p ++)
		{
			for (int m = 0; m < number_of_atoms/3; m ++)
			{
				if (oodifference[m + p*number_of_atoms/3] > 0.0 && oodifference[m + p*number_of_atoms/3] < 3.6)
				{
					hbondo_distance = 
				}
			}
		
		}
        }	
	
	cout << "\n\nYour H-bond histogram data has been placed in \"hbonds_histogram.dat\" and can now easily be plotted with gnuplot.\n\n";
}
else
{
	cout << "\n\nYou either entered the wrong key or you do not want to plot the H-bond histogram.\n\n";
}





inputfile.close();
oo_outputfile.close();
oh_outputfile.close();
angle_outputfile.close();
hbonds_outputfile.close();

return 0;
}
