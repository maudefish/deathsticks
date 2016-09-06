#ifndef RANDGAUSS_H
#define RANDGAUSS_H
#include <iostream> 
#include <cmath>
#include "randnew.h"
using namespace std;

double randgauss(double mean, double var)
{

	double fac, rsq, v1, v2;
	do {
		v1 = randnew(-1,1) ;
		v2 = randnew(-1,1) ;
		rsq = v1*v1 + v2*v2;
		//cout << rsq << endl;
		} while ((rsq >= 1.0) || (rsq == 0.0));
		fac = sqrt(-2.0*log(rsq)/rsq);
		return  v2*fac*sqrt(var) + mean;
}
#endif