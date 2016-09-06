#ifndef EXTRACT_H
#define EXTRACT_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
 
/* This program takes in a .txt filename and returns its contents in a vector */

vector<double> extract(string filename)
{
	filename += ".txt";
	ifstream myreader(filename);
	bool open_or_not = myreader.is_open();

	if (open_or_not == false) {printf("\n\tUNABLE TO ACCESS %s\n", filename.c_str());}

	string datum;
	vector<double> data;

	while(getline(myreader, datum))
	{
		data.push_back(stod(datum));
		//cout << datum << endl;
	}
	
	if (open_or_not == true) {printf("\n\tData from %s read in successfully.\n", filename.c_str());}

	return data;
}

#endif 
