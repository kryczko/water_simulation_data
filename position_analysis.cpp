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

/*bool ohdecision(double ox1, double oy1, double oz1, double ox2, double oy2, double oz2, double hx[], double hy[], double hz[], int h_indices[], double lattice)
{
        bool ohtruth(false), angletruth(false), finaltruth;
        for (int i = 0; i < 3; i ++)
        {
                if (h_indices[i] != 0)
                {       
                        double dxoh = ox1 - hx[h_indices[i]];
                        double dyoh = oy1 - hy[h_indices[i]];
                        double dzoh = oz1 - hz[h_indices[i]];
        
                        double dxoo = ox1 - ox2;
                        double dyoo = oy1 - oy2;
                        double dzoo = oz1 - oz2;                                
        
                        dxoh -= lattice*pbc_round(dxoh/lattice);
                        dyoh -= lattice*pbc_round(dyoh/lattice);
                        dzoh -= lattice*pbc_round(dzoh/lattice);
                        dxoo -= lattice*pbc_round(dxoo/lattice);
                        dyoo -= lattice*pbc_round(dyoo/lattice);
                        dzoo -= lattice*pbc_round(dzoo/lattice);

                        double inner_prod = dxoh*dxoo + dyoh*dyoo + dzoh*dzoo;
                        double oodistance = sqrt (dxoo*dxoo + dyoo*dyoo + dzoo*dzoo ); 
                        double ohdistance = sqrt ( dxoh*dxoh + dyoh*dyoh + dzoh*dzoh );
                        double angle = acos ( inner_prod / (oodistance*ohdistance) ) * 57.2957795;
			

                        if (ohdistance > 0. && ohdistance < 2.4 && angle < 30.0 )
                        {
                                ohtruth = true;
                        }
                       
                }
        }
        return ohtruth;
}*/



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
ofstream oo_outputfile, oh_outputfile, angle_outputfile, hbonds_outputfile, rmsd_outputfile;

char oodistance, ohdistance, angle, hbonds, rmsd; 
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
cout << "If you would like to create a plot of the mean squared displacement with respect to time enter \"y\", if not enter \"n\": ";
cin >> rmsd;
//cout << "If you would like to create a histogram of the number of H-bonds with respect to the z-axis enter \"y\", if not enter \"n\": ";
//cin >> hbonds;
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
		oo_outputfile << i/10. << "\t" << bin[i]/126000 << endl;
		oo_outputfile << (i + 1)/10. << "\t" << bin[i]/126000 << endl;
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
		oh_outputfile << j/100. << "\t" << bin[j]/256000. << endl;
		oh_outputfile << (j + 1)/100. << "\t" << bin[j]/256000. << endl;
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
		angle_outputfile << j/100. << "\t" << bin[j]/128000. << endl;
		angle_outputfile << (j + 1)/100. << "\t" << bin[j]/128000. << endl;
	}

        cout << "\n\nYour H-O-H angle histogram data has been placed in \"angle_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
        cout << "\n\nYou either entered the wrong key or you do not want to plot the H-O-H angle.\n\n";
}

//COMPUTES THE MEAN SQUARE DISTANCE OF THE WATER MOLECULES

if (rmsd == 'y')
{
	rmsd_outputfile.open("msd.dat");

	int nooa = number_of_atoms/3;
	double oxyz[nooa*timesteps][3], distance[nooa*timesteps], rmsdistance[timesteps], sum[timesteps], COM[timesteps][3];

	for (int i = 0; i < timesteps; i ++)
	{
		for (int j = 0; j < nooa; j ++)
		{
			oxyz[j + i*nooa][0] = lattice*x[j + i*number_of_atoms];
                        oxyz[j + i*nooa][1] = lattice*y[j + i*number_of_atoms];
                        oxyz[j + i*nooa][2] = lattice*z[j + i*number_of_atoms];
                }
	}
	
	for(int i = 0; i < timesteps; i ++)
	{
		double dx(0), dy(0), dz(0);
		for (int j = 0; j < nooa; j ++)
		{
			dx += oxyz[j + i*nooa][0];
			dy += oxyz[j + i*nooa][1];
			dz += oxyz[j + i*nooa][2];
		}
		
		dx /= nooa;
		dy /= nooa;
		dz /= nooa;
	
		COM[i][0] = dx;
		COM[i][1] = dy;
		COM[i][2] = dz;
	}
	
	for (int i = 0; i < timesteps; i ++)
	{	
		for (int j = 0; j < nooa; j ++)
		{
			oxyz[j + i*nooa][0] -= COM[i][0];
			oxyz[j + i*nooa][1] -= COM[i][1];
			oxyz[j + i*nooa][2] -= COM[i][2];
	
//			oxyz[j + i*nooa][0] -= lattice*pbc_round(oxyz[j + i*nooa][0]/lattice);
  //                      oxyz[j + i*nooa][1] -= lattice*pbc_round(oxyz[j + i*nooa][1]/lattice);
    //                    oxyz[j + i*nooa][2] -= lattice*pbc_round(oxyz[j + i*nooa][2]/lattice);
		}
	}
	
 	for (int i = 0; i < timesteps; i ++)
	{
		for (int j = 0; j < nooa; j ++)
		{
			double dxa = oxyz[j + i*nooa][0] - oxyz[0][0];
			double dya = oxyz[j + i*nooa][1] - oxyz[0][1];
			double dza = oxyz[j + i*nooa][2] - oxyz[0][2];
			
			distance[j + i*nooa] = dxa*dxa + dya*dya + dza*dza ;
		}
	}
		
	for (int i = 0; i < timesteps; i ++)
	{
		double add(0);	
		for (int j = 0; j < nooa; j ++)
		{
			add += distance[j + i*nooa];
			
		}
		sum[i] = add/nooa;
	}
	for (int i = 0; i < timesteps; i ++)
	{
		rmsd_outputfile <<  sum[i] << endl;
	} 

cout << "\n\nYour rms displacement data has been placed in \"msd.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
cout << "\n\nYou either entered the wrong key or you do not want to plot the rms displacement.\n\n";
}
/*if (hbonds == 'y')
{
	hbonds_outputfile.open("hbonds_histogram.dat");
	
	double oxyz[number_of_atoms/3][4], hx[2*number_of_atoms/3], hy[2*number_of_atoms/3], hz[2*number_of_atoms/3];
	double ohdifference[2*number_of_atoms*number_of_atoms/9], oodifference[number_of_atoms*number_of_atoms/9], nhlist[number_of_atoms/3][5], nolist[number_of_atoms/3][5], bin[130];
	double dxo, dyo, dzo, dxh1, dxh2, dyh1, dyh2, dzh1, dzh2, nhbs = 0, sum = 0;
	int bin_number;	

	for (int i = 0; i < timesteps; i ++)
        {
                for (int j = 0; j < number_of_atoms/3; j ++)
                {
                        oxyz[j][0] = j;
			oxyz[j][1] = lattice*x[j + i*number_of_atoms];
                        oxyz[j][2] = lattice*y[j + i*number_of_atoms];
                        oxyz[j][3] = lattice*z[j + i*number_of_atoms];
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
		
				dxh1 = oxyz[u][1] - hx[2*g];
				dxh2 = oxyz[u][1] - hx[2*g + 1];
				dyh1 = oxyz[u][2] - hy[2*g];
				dyh2 = oxyz[u][2] - hy[2*g + 1];
				dzh1 = oxyz[u][3] - hz[2*g];
				dzh2 = oxyz[u][3] - hz[2*g + 1];

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
				if (ohdifference[w + q*2*number_of_atoms/3] < 1.15)
				{	
					nhlist[q][d] = w ;
					d ++;
				}
			}
		}
		
		for ( int p = 0; p < number_of_atoms/3; p ++)
		{
			for (int m = 0; m < number_of_atoms/3; m ++)
			{
				dxo = oxyz[p][1] - oxyz[m][1];
				dyo = oxyz[p][2] - oxyz[m][2];
				dzo = oxyz[p][3] - oxyz[m][3];					
				
				dxo -= lattice*pbc_round(dxo/lattice);
                                dyo -= lattice*pbc_round(dyo/lattice);
                                dzo -= lattice*pbc_round(dzo/lattice);
				
				oodifference[m + p*number_of_atoms/3] = sqrt ( dxo*dxo + dyo*dyo + dzo*dzo );
			}
		}
		
		for (int oindex = 0; oindex < number_of_atoms/3 ; oindex ++)
		{
			nhbs = 0;
			for (int hindex = 0; hindex < number_of_atoms/3; hindex ++)
			{
				if (oodifference[hindex + oindex*number_of_atoms/3] > 0.0 && oodifference[hindex + oindex*number_of_atoms/3] < 3.6)
				{
					int h_indices[4] = { nhlist[hindex][1], nhlist[hindex][2], nhlist[hindex][3], nhlist[hindex][4] };
					
					nhbs += ohdecision(oxyz[oindex][1], oxyz[oindex][2], oxyz[oindex][3], oxyz[hindex][1], oxyz[hindex][2], oxyz[hindex][3], hx, hy, hz, h_indices, lattice);
				}
			}
			if (nhlist[oindex][3] == 0 && nhlist[oindex][4] == 0)
                                        {
                                                nhbs += 2;
                                        }
                                        if (nhlist[oindex][3] != 0 && nhlist[oindex][4] == 0)
                                        {
                                                nhbs += 3;
                                        }
                                        if (nhlist[oindex][4] != 0)
                                        {
                                                nhbs += 4;
                                        }

		sum += nhbs;
//		int bin_number = oxyz[oindex][3]*10;
//		bin[bin_number] = nhbs;
		}	
        }	
	cout << sum  << endl;
//	for (int i = 0; i < 130; i ++)
//	{
//		hbonds_outputfile << i << "\t" << bin[i] << endl;
//		hbonds_outputfile << i + 1 << "\t" << bin[i] << endl;
//	}
	
	cout << "\n\nYour H-bond histogram data has been placed in \"hbonds_histogram.dat\" and can now easily be plotted with gnuplot.\n\n";
}
else
{
	cout << "\n\nYou either entered the wrong key or you do not want to plot the H-bond histogram.\n\n";
}
*/




inputfile.close();
oo_outputfile.close();
oh_outputfile.close();
angle_outputfile.close();
hbonds_outputfile.close();
rmsd_outputfile.close();
return 0;
}
