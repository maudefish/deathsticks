#ifndef CPMAIN_UPPERHEM_AVG_H
#define CPMAIN_UPPERHEM_AVG_H
#include <string>
#include <vector>
using namespace std;

class newbin
{
public:
	vector<double> e, x, y, z, l, b;
	int size;
};

void makevectors(vector<double>& l, vector<double>& b, vector<double>& x, vector<double>& y, vector<double>& z); 

void makebins(newbin& bin1, newbin& bin2, newbin& bin3, newbin& bin4, newbin& bin5, vector<double>& e, vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& l, vector<double>& b); 

void calculateQ(newbin& bin1, newbin& bin2, newbin& bin3, newbin& bin4, newbin& bin5, vector<double>& q_e10_40, vector<double>& q_e10_30, vector<double>& q_e10_20, vector<double>& q_e20_40, vector<double>& q_e20_30, vector<double>& q_e30_40, vector<double>& delta_e10_40, vector<double>& delta_e10_30, vector<double>& delta_e10_20, vector<double>& delta_e20_40, vector<double>& delta_e20_30, vector<double>& delta_e30_40);


#endif
