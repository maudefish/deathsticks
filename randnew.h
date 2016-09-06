#ifndef RANDNEW_H
#define RANDNEW_H
#include <cstdlib>	

double randnew(double a, double b) //generates random double between a and b 
{
	return a + (rand() / (RAND_MAX + 1.0))*(b-a);
}
#endif 

// a) main must include <ctime> and <cstdlib>
// b) seed in mainwith " srand(time(NULL)*getpid()); "