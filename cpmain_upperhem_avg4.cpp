#include <iostream>
#include <vector>
#include <string>
#include <cmath>
//#include "getfilename_upperhem_MPImock.h"
#include "extract.h"
#include "writetotxt.h"
//#include "convertcoordinates.h"
#include "cpmain_upperhem_avg.h"	
using namespace std;

/* NOTE: this macro's directory must contain two subdirectories: '\output' (empty) and '\coordinates' (with the L, B, and E .txt files) */
/* this macro's job is to take 3 files (L, B, and E of the whole upper hemisphere) and compute Q's */

/* NOTE II: This version of cp_main uses the AVERAGED E1 and E2 positions (mathematically equivalent to what I was doing before!!) */

int main() 
{
	string inprefix = "coordinates/";
	string outprefix = "output/";

	string e_filename = inprefix + "e";
	string l_filename = inprefix + "l";
	string b_filename = inprefix + "b";
	vector<double> e = extract(e_filename);
	vector<double> l = extract(l_filename);
	vector<double> b = extract(b_filename);

	/* now create the unit vectors */ 
	vector<double> x, y, z;
	makevectors(l, b, x, y, z);

	/* now we bin the data */
	newbin bin1, bin2, bin3, bin4, bin5;
	makebins(bin1, bin2, bin3, bin4, bin5, e, x, y, z, l, b);

	/* now we can actually do the statistics */
	vector<double> q_e10_40, q_e10_30, q_e10_20, q_e20_40, q_e20_30, q_e30_40;
	vector<double> delta_e10_40, delta_e10_30, delta_e10_20, delta_e20_40, delta_e20_30, delta_e30_40;
	//vector<double> n_e10_40, n_e10_30, n_e10_20, n_e20_40, n_e20_30, n_e30_40;
	
	calculateQ(bin1, bin2, bin3, bin4, bin5, q_e10_40, q_e10_30, q_e10_20, q_e20_40, q_e20_30, q_e30_40, delta_e10_40, delta_e10_30, delta_e10_20, delta_e20_40, delta_e20_30, delta_e30_40);

	/* write out values of Q to txt files */
	writetotxt(outprefix + "q_e10_40", q_e10_40);
	writetotxt(outprefix + "q_e10_30", q_e10_30);
	writetotxt(outprefix + "q_e10_20", q_e10_20);
	writetotxt(outprefix + "q_e20_40", q_e20_40);
	writetotxt(outprefix + "q_e20_30", q_e20_30);
	writetotxt(outprefix + "q_e30_40", q_e30_40);

	writetotxt(outprefix + "delta_e10_40", delta_e10_40);
	writetotxt(outprefix + "delta_e10_30", delta_e10_30);
	writetotxt(outprefix + "delta_e10_20", delta_e10_20);
	writetotxt(outprefix + "delta_e20_40", delta_e20_40);
	writetotxt(outprefix + "delta_e20_30", delta_e20_30);
	writetotxt(outprefix + "delta_e30_40", delta_e30_40);

return 1;}

void makevectors(vector<double>& l, vector<double>& b, vector<double>& x, vector<double>& y, vector<double>& z)
{
	for(int i=0, max=l.size(); i<max; i++)
	{
		x.push_back(cos(l[i]*M_PI/180.)*cos(b[i]*M_PI/180.));
		y.push_back(sin(l[i]*M_PI/180.)*cos(b[i]*M_PI/180.));
		z.push_back(sin(b[i]*M_PI/180.));
	}
}

void makebins(newbin& bin1, newbin& bin2, newbin& bin3, newbin& bin4, newbin& bin5, vector<double>& e, vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& l, vector<double>& /*b*/)
{
	double cut1 = 10000;
	double cut2 = 20000;
	double cut3 = 30000;
	double cut4 = 40000;
	double cut5 = 50000;
	double cut6 = 60000;

	for (int i=0, max=l.size(); i<max; i++)
	{
		if ((cut1 <= e[i]) && (e[i] < cut2))
		{
			bin1.x.push_back(x[i]);
			bin1.y.push_back(y[i]);
			bin1.z.push_back(z[i]);
		}
		else if ((cut2 <= e[i]) && (e[i] < cut3))
		{
			bin2.x.push_back(x[i]);
			bin2.y.push_back(y[i]);
			bin2.z.push_back(z[i]);
		}
		else if ((cut3 <= e[i]) && (e[i] < cut4))
		{
			bin3.x.push_back(x[i]);
			bin3.y.push_back(y[i]);
			bin3.z.push_back(z[i]);
		}
		else if ((cut4 <= e[i]) && (e[i] < cut5))
		{
			bin4.x.push_back(x[i]);
			bin4.y.push_back(y[i]);
			bin4.z.push_back(z[i]);
		}
		else if ((cut5 <= e[i]) && (e[i] < cut6))
		{
			bin5.x.push_back(x[i]);
			bin5.y.push_back(y[i]);
			bin5.z.push_back(z[i]);
		}
	}

	bin1.size = bin1.x.size();
	bin2.size = bin2.x.size();
	bin3.size = bin3.x.size();
	bin4.size = bin4.x.size();
	bin5.size = bin5.x.size();
}

void calculateQ(newbin& bin1, newbin& bin2, newbin& bin3, newbin& bin4, newbin& bin5, vector<double>& q_e10_40, vector<double>& q_e10_30, vector<double>& q_e10_20, vector<double>& q_e20_40, vector<double>& q_e20_30, vector<double>& q_e30_40, vector<double>& delta_e10_40, vector<double>& delta_e10_30, vector<double>& delta_e10_20, vector<double>& delta_e20_40, vector<double>& delta_e20_30, vector<double>& delta_e30_40)
{	
	double maxgamma = 60.;
	int maxradius = 30;
	newbin* binA = NULL;
	newbin* binB = NULL;
	newbin* binC = &bin5;

	for (int r=0; r<maxradius; r++)
	{

		for (int q=0; q<6; q++)  //6 Q's per radius
		{	
			int N1 = 0;
			int N2 = 0;
			int N3 = 0;
			double sum = 0.;
			double sum2 = 0.;
			double sigma = 0.;
			double delta = 0.;

			switch (q)
			{
			case 0: 
				binA = &bin1; binB = &bin4; 
				break;
			case 1: 
				binA = &bin1; binB = &bin3;
				break;
			case 2: 
				binA = &bin1; binB = &bin2;
				break;
			case 3: 
				binA = &bin2; binB = &bin4;
				break;
			case 4: 
				binA = &bin2; binB = &bin3;
				break;
			case 5: 
				binA = &bin3; binB = &bin4;
				break;
			default: 
				break;
			}
			
			for (int k=0; k<(binC->size); k++) 		/* LOOP OVER E3 */
			{
				double cx = binC->x[k];
				double cy = binC->y[k];
				double cz = binC->z[k];

				//printf("norm(n_c) = %10.7f\n",ax*ax + ay*ay + az * az );printf("norm(n_c) = %10.7f\n",cx*cx + cy*cy + cz * cz );

				double gamma = acos( cz ) * 180./M_PI;	

				/* 'gamma' is the angle b/w the z-axis and the E3 photon */

				if ( gamma <= maxgamma )
				{
					/* now that we've chosen a patch with gamma <= R, we can start averaging E1 and E2 positions */

					double ax = 0;
					double ay = 0;
					double az = 0;
					double bx = 0;
					double by = 0;
					double bz = 0;
					int N1_k = 0;
					int N2_k = 0;

					for (int i=0; i<(binA->size); i++)	/* Get avg. position of E1 w/in range of E3 */
					{
						double ax_temp = binA->x[i];
						double ay_temp = binA->y[i];
						double az_temp = binA->z[i];

						//printf("norm(n_a) = %10.7f\n",ax_temp*ax_temp + ay_temp*ay_temp + az_temp * az_temp );

						double alpha = acos( ( ax_temp * cx ) + ( ay_temp * cy ) + ( az_temp * cz ) ) * 180./M_PI;
						//cout << alpha << endl;

						if (alpha <= r+1)
						{
							ax += ax_temp;
							ay += ay_temp;
							az += az_temp;

							//printf("norm(n_c) = %10.7f\n",ax*ax + ay*ay + az * az );
							N1_k ++;
							//cout << alpha << endl;
						}
					}
					ax /= N1_k; 
					ay /= N1_k;
					az /= N1_k;
					//printf("N1_k = %10i\t norm(n_a) = %10.7f\n", N1_k, ax*ax + ay*ay + az * az );
					if (N1_k == 0)
					{
						ax = 0;
						ay = 0;
						az = 0;
					}

					for (int j=0; j<(binB->size); j++)		/* get avg. position of E2 w/in range of E3 */
					{
						double bx_temp = binB->x[j];
						double by_temp = binB->y[j];
						double bz_temp = binB->z[j];


						//printf("norm(n_b) = %10.7f\n",bx_temp*bx_temp + by_temp*by_temp + bz_temp * bz_temp );
						double beta = acos( ( bx_temp * cx ) + ( by_temp * cy ) + ( bz_temp * cz ) ) * 180./M_PI;

						if (beta <= r+1)
						{
							bx += bx_temp;
							by += by_temp;
							bz += bz_temp;
							N2_k ++; 
						}
					}

					bx /= N2_k;
					by /= N2_k;
					bz /= N2_k;
					
					if (N2_k == 0)
					{
						bx = 0;
						by = 0;
						bz = 0;
					} 

					/* calculate Q and Q^2 */

					sum += ( ax * by * cz ) - ( ax * bz * cy ) - ( ay * bx * cz ) + ( ay * bz * cx ) + ( az * bx * cy ) - ( az * by * cx );
					sum2 += (ax * by * cz - ax * bz * cy - ay * bx * cz + ay * bz * cx + az * bx * cy - az * by * cx) * (ax * by * cz - ax * bz * cy - ay * bx * cz + ay * bz * cx + az * bx * cy - az * by * cx);

					N1 += N1_k;
					N2 += N2_k;
					N3 ++;

					//cout << ax*ax + ay*ay + az * az << endl;
	
				}	/* end IF gamma */ 
			} 		/* end E3 (k)  */
					
			sum /= N3;
			sum2 /= N3;

			/* NOTE: This sigma is the population mean!!! We'll divide by sqrt(N3) again when glueing */
			sigma = sqrt( (sum2 - ( sum * sum )) / (N3-1.) );

			/* EDIT: Here's the reduced sigma (delta) */
			delta = sigma / sqrt(N3);

			if (N3 == 0) {sum = 0; sum2 = 0; sigma = 0; delta = 0;}
			if (N3 == 1) {sigma = 0; delta = 0;}

			printf("RADIUS: %i\t #%i\t N1: %i\t N2: %i\t N3: %i\t Q: %10.7g\t σ: %10.7g\t δ: %10.7g\n", r+1, q+1, N1, N2, N3, sum, sigma, delta);

			switch (q)
			{
			case 0: 
				q_e10_40.push_back(sum);
				delta_e10_40.push_back(delta);
				break;
			case 1: 
				q_e10_30.push_back(sum);
				delta_e10_30.push_back(delta);
				break;
			case 2: 
				q_e10_20.push_back(sum);
				delta_e10_20.push_back(delta);
				break;
			case 3: 
				q_e20_40.push_back(sum);
				delta_e20_40.push_back(delta);
				break;
			case 4: 
				q_e20_30.push_back(sum);
				delta_e20_30.push_back(delta);
				break;
			case 5: 
				q_e30_40.push_back(sum);
				delta_e30_40.push_back(delta);
				break;
			default: 
				cout << endl; 
				break;
			}			
		}	/* end loop over q */
		printf("===================================================================================================\n");	
	}	
}

