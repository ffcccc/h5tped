//#pragma warning( C4996)
// stl
#include <vector>
#include <string>
#include <iostream>
// options

// h5 SNP file format
#include "../hdfvar/hdfTPED.h"

#define THDF5 hdfTPED<char, float>

// main
int test()
{
	// initialization
	vector<string> alls;
	alls.clear();
	alls.push_back("AA");	// -> 0
	alls.push_back("AB");	// -> 1
	alls.push_back("BB");	// -> 2
	
	vector <THDF5::T_GType> codes;
	//char delimiter('\t');
	std::string filename("bladder-metabolici-riparo.h5");
	std::string strNA("NC");
	
	hdfTPED<char, float> *hf = new hdfTPED<char, float> (filename);
	hf->loadAllele(alls);

	unsigned int nSNPs, nSamples, dimGTypes;

	// Get the extent of the first (and only) dimension in the dataset 
	dimGTypes = hf->numGTypeCodes();
	nSNPs = hf->numSnps();
	nSamples = hf->numSamples();
	cout << "file: " << filename << " array extent: " << nSNPs << " snps and " << nSamples << " samples\n";
	//cout << hf.getSNPSummary(xyz);
	hf->getSampleSummary();
	delete hf;
	
	return 1;
}

int main(int argc, char* argv[])
{
	int res = test();
	return res;
}

