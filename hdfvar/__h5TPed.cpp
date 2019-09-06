/*                                                                                */
/*   Copyright  (C) 2005 Intel Corporation.  All rights reserved.                 */
/*                                                                                */
/*   The information and source code contained herein is the exclusive property   */
/*   of Intel Corporation and may not be disclosed, examined, or reproduced in    */
/*   whole or in part without explicit written authorization from the Company.    */
/*                                                                                */
/*                                                                                */
#ifdef WIN32

#include "time/mytime.h"
#include "h5TPed.h"
#include <iostream>
#include <fstream>
#include "AvCrt.h"


using namespace std;

//-----------------------------------------------------------------------------------------------------------
//                  h5TPED implemented with use of Intel Array Viewer library
//-----------------------------------------------------------------------------------------------------------
h5TPED::h5TPED(const std::string &szFilename) {
	m_filename = szFilename;
	BOOL bOk;
	// looking for filename.h5
	bOk = avOpen(filename().c_str(), FALSE);	
	if (!bOk) {
		string msg(avGetErrorMessage());
		cout << "\nCouldn't open file: "<< szFilename << "\nmsg: " << msg;
		return;
	}
	m_nChr = 0;
	m_nSnp = avGetExtent("/SNPInfoTable", 1);	// table 1xN
	m_nSamples = avGetExtent("/SampleTable", 1);// table 1xM
	int r = avGetExtent("/SNPDataTable", 1);	// table MxN or NxM
	int c = avGetExtent("/SNPDataTable", 2);
	m_SamplesInRows = (m_nSamples == r);
	m_NA = "NA";
	T_GType **pData0 = (T_GType **)avGetArrayPointer("/SNPDataTableInv");
	if (pData0 == NULL) {
		doInvDataTable();
	}
}

h5TPED::~h5TPED(void){
	BOOL bOk = avSave(filename().c_str());
	if (bOk) {
		printf("Data saved to: %s\n", filename().c_str());
	} else {
		printf("Error saving to %s\n", filename().c_str());
	}
	avClose();
}


int h5TPED::doInvDataTable(){
	BOOL bOk;
	int dims2[2];
	T_GType **p2DArray = 0;

	// Create an inverted 2D array
	if(m_SamplesInRows){
		dims2[1] = numSamples(); // #rows
		dims2[0] = numSnps(); // #cols
	} else {
		dims2[1] = numSnps(); // #rows
		dims2[0] = numSamples(); // #cols
	}

	p2DArray = (T_GType **)avAlloc(2, dims2, AV_CHAR_TYPE, "SNPDataTableInv");
	bOk = avSave(filename().c_str());
	
	T_GType val;
	for(int i=0;i<dims2[0];i++){
		for(int j=0;j<dims2[1];j++){
			this->getData("/SNPDataTable", j, i, val);
			this->setData("/SNPDataTableInv", i, j, val);
		}
		if((i%1000)==0) 
			cout << "#";
	}
	bOk = avSave(filename().c_str());
	return bOk;
}

int h5TPED::doDataTable(){
	BOOL bOk;
	int dims2[2];	// = {300000, 96};
	T_GType **p2DArray = 0;

	// Create a 2D array
	if(m_SamplesInRows){
		dims2[0] = numSamples(); // #rows
		dims2[1] = numSnps(); // #cols
	} else {
		dims2[0] = numSnps(); // #rows
		dims2[1] = numSamples(); // #cols
	}
	cout << sizeof(T_GType) << std::endl;

	p2DArray = (T_GType **)avAlloc(2, dims2, AV_CHAR_TYPE, "SNPDataTable");
	if (!p2DArray) {
		char msg[255];
		strcpy(msg, avGetErrorMessage());
		cout <<"couldn't create data table: " << std::endl << "msg: " << msg << std::endl;
		return false;
	}
	bOk = avUpdate(p2DArray);
	bOk = avSave(filename().c_str());
	
	/*avClose();*/
	return bOk;
}


int h5TPED::doSampleTable(){
	BOOL bOk;
	int hS1Type, hStringType12;
	T_Sample *s1;  // n-element array of type s1_t
	int n = numSamples();
	
	// Create string types
	hStringType12 = avCreateAsciiStringType(12);
	assert(hStringType12 >= 0);

	// Create a compound type for s1_t
	hS1Type = avCreateType(sizeof(T_Sample));
	assert(hS1Type >= 0);

	// Add elements to the type
	// Use the offsetof macro to get the offsets into the type
	bOk = avInsertType(hS1Type, "index", offsetof(T_Sample, index), AV_LONG_TYPE, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "name", offsetof(T_Sample, name), hStringType12, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "famIndex1", offsetof(T_Sample, famIndex1), AV_LONG_TYPE, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "famIndex2", offsetof(T_Sample, famIndex2), AV_LONG_TYPE, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "phenotype", offsetof(T_Sample, phenotype), AV_FLOAT_TYPE, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "sex", offsetof(T_Sample, sex), AV_LONG_TYPE, 1);
	assert(bOk);
	s1 = (T_Sample *)avAlloc(1, &n, hS1Type, "SampleTable");

	cout << "creating the file: " << filename() << std::endl;
	bOk = avUpdate(s1);
	bOk = avSave(filename().c_str());
	if (bOk) {
		cout <<"Data saved to: " << filename() << std::endl;
	} else {
		cout <<"Error saving to " << filename() << std::endl;
	}
	/*avClose();*/
	cout <<"writeInfoTable ended." << std::endl;
	return bOk;
}


int h5TPED::doSNPTable(){
	BOOL bOk;
	int hS1Type, hStringType12, hStringType3;
	T_SNP *s1;  // n-element array of type s1_t
	int n = numSnps();
	
	// Create string types
	hStringType12 = avCreateAsciiStringType(12);
	assert(hStringType12 >= 0);
	hStringType3 = avCreateAsciiStringType(3);
	assert(hStringType3 >= 0);
	
	// Create a compound type for s1_t
	hS1Type = avCreateType(sizeof(T_SNP));
	assert(hS1Type >= 0);

	// Add elements to the type
	// Use the offsetof macro to get the offsets into the type
	bOk = avInsertType(hS1Type, "index", offsetof(T_SNP, index), AV_LONG_TYPE, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "name", offsetof(T_SNP, name), hStringType12, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "chr", offsetof(T_SNP, chr), hStringType3, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "position", offsetof(T_SNP, position), AV_LONG_TYPE, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "score", offsetof(T_SNP, score), AV_FLOAT_TYPE, 1);
	assert(bOk);
	bOk = avInsertType(hS1Type, "allele", offsetof(T_SNP, allele), hStringType3, 1);
	assert(bOk);
		
	s1 = (T_SNP *)avAlloc(1, &n, hS1Type, "SNPInfoTable");

	cout << "creating the file: " << filename() << std::endl;
	bOk = avUpdate(s1);
	bOk = avSave(filename().c_str());
	if (bOk) {
		cout <<"Data saved to: " << filename() << std::endl;
	} else {
		cout <<"Error saving to " << filename() << std::endl;
	}
	/*avClose();*/
	cout <<"writeInfoTable ended." << std::endl;
	return bOk;
}

// create a new file and its data structure

bool h5TPED::buildStruct(){
	BOOL bOk = false;
	
	// looking for filename.h5
	bOk = avOpen(filename().c_str(), TRUE);
	if (!bOk) {
		char msg[255];
		strcpy(msg, avGetErrorMessage());
		cout <<"couldn't open file: " << filename() << std::endl << "msg: " << msg << std::endl;
		return false;
	}
	//bOk = doDataTable();
	bOk = bOk && doSNPTable();
	bOk = bOk && doSampleTable();
	bOk = bOk && avSave(filename().c_str());
	/*avClose();*/
	return (bOk > 0);
}

// public members

// setters

void h5TPED::setData(const std::string &table, const int row, const int col, const T_GType val){
	//cout << "Getting pointer to " << table << "\n";
	// if we know that dset1.1.2 is one-dimensional and of type AV_LONG, 
	// so cast the pointer returned by avGetArrayPointer to a pointer to long
	//if(table == "/SNPDataTable"){
		T_GType **pData0 = (T_GType **)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		pData0[row][col] = val;
		avUpdate(pData0);
	//}
}


void h5TPED::setData(const std::string &table, const int row, const T_Sample &val){
	//cout << "Getting pointer to " << table << "\n";
	if(table == "/SampleTable"){
		T_Sample *pData0 = (T_Sample *)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		pData0[row] = val;
		avUpdate(pData0);
	}
}


void h5TPED::setData(const std::string &table, const int row, const T_SNP &val){
	//cout << "Getting pointer to " << table << "\n";
	if(table == "/SNPInfoTable"){
		T_SNP *pData0 = (T_SNP *)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		pData0[row] = val;
		avUpdate(pData0);
	}
}

void h5TPED::setData(const std::string &table, const int row, const long &val){
	//cout << "Getting pointer to " << table << "\n";
	if(table == "/EntrezInfoTable"){
		long *pData0 = (long *)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		pData0[row] = val;
		avUpdate(pData0);
	}
}
// getters

void h5TPED::getData(const std::string &table, const int row, const int col, T_GType &val) const {
	//cout << "Getting pointer to " << table << "\n";
	// if we know that dset1.1.2 is one-dimensional and of type AV_LONG, 
	// so cast the pointer returned by avGetArrayPointer to a pointer to long
	//if(table == "/SNPDataTable"){
		T_GType **pData0 = (T_GType **)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		val = pData0[row][col];
	//}
}


void h5TPED::getSnpPtr(const int row, T_GType *&val, const std::string &table) const {
	//cout << "Getting pointer to " << table << "\n";
	// if we know that dset1.1.2 is one-dimensional and of type AV_LONG, 
	// so cast the pointer returned by avGetArrayPointer to a pointer to long
	//if(table == "/SNPDataTable"){
		//table += "/train";
		T_GType **pData0 = (T_GType **)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		val = pData0[row];
	//}
}

//void h5TPED::getSnpRow(const std::string &table, const int row, valarray<T_GType> &val) const {
//	//cout << "Getting pointer to " << table << "\n";
//	// if we know that dset1.1.2 is one-dimensional and of type AV_LONG, 
//	// so cast the pointer returned by avGetArrayPointer to a pointer to long
//	//if(table == "/SNPDataTable"){
//		T_GType **pData0 = (T_GType **)avGetArrayPointer(table.c_str());
//		if (pData0 == NULL) {
//			cout << "couldn't get pointer to array\n";
//			return;
//		}
//		size_t sz = this->numSamples();
//		if(m_SamplesInRows){
//			sz = this->numSnps();
//		}
//		val.resize(sz);
//		for(unsigned int i = 0; i<sz; i++)
//			val[i] = pData0[row][i];
//	//}
//}
//
//
//void h5TPED::getSnpCol(const std::string &table, const int col, valarray<T_GType> &val) const {
//	//cout << "Getting pointer to " << table << "\n";
//	// if we know that dset1.1.2 is one-dimensional and of type AV_LONG, 
//	// so cast the pointer returned by avGetArrayPointer to a pointer to long
//	if(table == "/SNPDataTable"){
//		// questo deve puntare alla riga 'col' della tabella inversa
//		T_GType **pData0 = (T_GType **)avGetArrayPointer(table.c_str());
//		if (pData0 == NULL) {
//			cout << "couldn't get pointer to array\n";
//			return;
//		}
//		size_t sz = this->numSamples();
//		if(m_SamplesInRows){
//			sz = this->numSnps();
//		}
//		val.resize(sz); 
//		for(unsigned int i = 0; i<sz; i++)
//			val[i] = pData0[i][col];
//	}
//}

void h5TPED::getData(const std::string &table, const int row, long &val) const {
	if(table == "/EntrezInfoTable"){
		long *pData0 = (long *)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		val = pData0[row];
	}
}

void h5TPED::getData(const std::string &table, const int row, T_SNP &val) const {
	//cout << "Getting pointer to " << table << "\n";
	if(table == "/SNPInfoTable"){
		T_SNP *pData0 = (T_SNP *)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		val = pData0[row];
	}
}


void h5TPED::getData(const std::string &table, const int row, T_Sample &val) const {
	//cout << "Getting pointer to " << table << "\n";
	if(table == "/SampleTable"){
		T_Sample *pData0 = (T_Sample *)avGetArrayPointer(table.c_str());
		if (pData0 == NULL) {
			cout << "couldn't get pointer to array\n";
			return;
		}
		val = pData0[row];
	}
}
#endif
