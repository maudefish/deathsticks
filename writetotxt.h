#ifndef WRITETOTXT_H
#define WRITETOTXT_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

/* This program writes data from a vector to a .txt file */
void writetotxt(string filename, const vector<double>& data)
{
	
	filename += ".txt";
	ofstream myscribe(filename);
	for (int i=0, max=data.size(); i<max; i++)
	{
		myscribe << data[i] << "\n";
	}

	/* attempt to open the file to check if it's actually there!! */
	ifstream check(filename);		
	if (check.is_open() == false)		
	{
		printf("\n\tUnable to write to %s!!\n", filename.c_str());
	}
	else {printf("\n\tData written to %s successfully.\n", filename.c_str());}
}
 
#endif
