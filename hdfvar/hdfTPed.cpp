/*                                                                                */
/*   Copyright  (C) 2005 Intel Corporation.  All rights reserved.                 */
/*                                                                                */
/*   The information and source code contained herein is the exclusive property   */
/*   of Intel Corporation and may not be disclosed, examined, or reproduced in    */
/*   whole or in part without explicit written authorization from the Company.    */
/*                                                                                */
/*                                                                                */

//#include "time/mytime.h"
#include "hdfTPED.h"
#include <iostream>
#include <fstream>
#include <H5File.h>
using namespace std;
#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif


//#include "h5cpputil.h"	// C++ utilility header file
//-----------------------------------------------------------------------------------------------------------
//                           hdfTPED implemented with direct use of hdf library
//-----------------------------------------------------------------------------------------------------------

/*******************************************************/
/* Check existence of /MyGroup/Group_A                 */
/*******************************************************/
template<typename T_GType, typename T_ExpType>
bool hdfTPED<T_GType, T_ExpType>::objExists(const string &path) const {
	hid_t       file_id = this->getId();
	herr_t      status;

  //status = H5Eset_auto(file_id, NULL, NULL);
  //status = H5Gget_objinfo (file_id, path.c_str(), 0, NULL);
	status = H5Lexists(file_id, path.c_str(), H5P_DEFAULT);
  //printf ("/MyGroup/Group_A: ");
  if (status > 0){
		//printf ("The group exists.\n");
		return true;
	}
	//printf ("The group either does NOT exist\n or some other error occurred.\n"); 
	return false;
}

template<typename T_GType, typename T_ExpType>
hdfTPED<T_GType, T_ExpType>::hdfTPED(const std::string &szFilename)	:	h5TPED_I<T_GType, T_ExpType>(szFilename) {
	tblSNPInfo=NULL;
	tblSample=NULL;
	tblSNPData=NULL;
	tblSNPDataInv=NULL;
	tblRefSeq=NULL;
	tblExpData=NULL;
	string currentPath("");
	DataSpace ds;
	hsize_t		size1[1];
	hsize_t		size2[2];
	long r = 0;
	long c = 0;
	try{
		//m_nSnp = avGetExtent("/tblSNPInfo", 1);	// table 1xN

		currentPath = this->SNPInfoPath();
		if(this->objExists(currentPath)) {
			tblSNPInfo = new DataSet(this->openDataSet(currentPath));
			ds = tblSNPInfo->getSpace();
			ds.getSimpleExtentDims(size1);
			this->setNumSnps(*size1);
			hsize_t	size = this->numSnps();
			spcSNPInfo = new DataSpace(1, &size, NULL);
		}
		else {
			cout << "warning: no dataset with SNPInfo data." << endl;
		}

		//m_nSamples = avGetExtent("/SampleTable", 1);// table 1xM
		currentPath = this->SNPInfoPath();
		if(this->objExists(currentPath)) {
			tblSample = new DataSet(this->openDataSet(currentPath));
			ds = tblSample->getSpace();
			ds.getSimpleExtentDims(size1);
			this->setNumSamples(*size1);
			hsize_t size = this->numSamples();
			spcSample = new DataSpace(1, &size, NULL);
		}
		else {
			cout << "warning: no dataset with Sample data." << endl;
		}

		currentPath = this->SNPDataPath();
		if(this->objExists(currentPath)) {tblSNPData = new DataSet(this->openDataSet(currentPath));
			ds = tblSNPData->getSpace();
			ds.getSimpleExtentDims(size2);
			r = size2[0];
			c = size2[1];
			//hsize_t tmp = size2[0];
			spcSNPData = new DataSpace(2, size2, NULL);
		}
		else {
			cout << "warning: no dataset with SNP data." << endl;
		}
	
		// linear memory dataspaces
		hsize_t	ms_size[2];
		// 1 x M
		ms_size[0] = this->numSamples();
		spaceMemNSamp = new DataSpace(1, ms_size, NULL);
		// 1 x N
		ms_size[0] = this->numSnps();
		spaceMemNSnp = new DataSpace(1, ms_size, NULL);
		// 1 x 1
		ms_size[0] = 1;
		spaceMem1 = new DataSpace(1, ms_size, NULL);

		currentPath = this->ExpDataPath();
		if(this->objExists(currentPath)) {
			tblExpData = new DataSet(this->openDataSet(currentPath));
			ds = tblExpData->getSpace();
			ds.getSimpleExtentDims(size2);
			this->setNumGenes(size2[0]);
			assert(this->numSamples() == size2[1]);
			spcExpData = new DataSpace(2, size2, NULL);
			// usa solo ind[0]
			tblRefSeq = new DataSet(this->openDataSet(this->RefSeqInfoPath()));
			spcRefSeq = new DataSpace(1, size2, NULL);
		}
		else {
			cout << "warning: no dataset with SNP Inv data." << endl;
		}
		
		currentPath = this->SNPDataInvPath();
		if(this->objExists(currentPath)) {
			tblSNPDataInv = new DataSet(this->openDataSet(currentPath));
			//tblSNPDataInv = new DataSet (createDataSet("/SNPDataTableInv", NATIVE_GTYPE, *spcSNPDataInv));
			// Add a comment to the dataset
			//setComment ("SNPDataTableInv", "This is the SNP dataset ordered by SNP-rows");
			//doInvDataTable();
		}
		else {
			cout << "warning: no dataset with SNP Inv data." << endl;
		}
	}
	catch(Exception E){
		tblSNPInfo = NULL;
		cout << "error: opening dataset:" << currentPath << endl;
	}

	size2[0] = this->numSnps();
	size2[1] = this->numSamples();
	spcSNPDataInv = new DataSpace(2, size2, NULL);

	#define SZ_HBUFF 1000
	tconv_buf = new char [SZ_HBUFF];
	xfer.setBuffer (SZ_HBUFF, tconv_buf, NULL);

	this->setSamplesByRow(this->numSamples() == r);
	this->setNA("NA");

	buildHdfTypes();
}

template<typename T_GType, typename T_ExpType>
hdfTPED<T_GType, T_ExpType>::hdfTPED(const string &szFilename, const int nsnp, const int nsamp, const int nexpr, const bool SamplesInRows)
	:	h5TPED_I<T_GType, T_ExpType>(szFilename, nsnp, nsamp, SamplesInRows) {
	
	// set sizes
	this->setNumSamples(nsamp);
	this->setNumSnps(nsnp);
	this->setNumGenes(nexpr);
	this->setSamplesByRow(SamplesInRows);
	this->setNA("NA");
	
	// init tables
	tblSNPInfo=NULL;
	tblSample=NULL;
	tblSNPData=NULL;
	tblSNPDataInv=NULL;
	tblRefSeq=NULL;
	tblExpData=NULL;
	
	// buffers & complex types
	buildHdfTypes();
	
	//
	hsize_t		size1[1];
	hsize_t		size2[2];
	hsize_t		size2Inv[2];
	// genotype by col
	if(SamplesInRows){
		size2[0] = nsamp;
		size2[1] = nsnp;
	} else {
		size2[1] = nsamp;
		size2[0] = nsnp;
	}
	spcSNPData = new DataSpace(2, size2, NULL);

	// genotype by row
	size2Inv[1] = size2[0];
	size2Inv[0] = size2[1];
	spcSNPDataInv = new DataSpace(2, size2Inv, NULL);

	// snp info
	size1[0] = nsnp;
	spcSNPInfo = new DataSpace(1, size1, NULL);

	// samples
	size1[0] = nsamp;
	spcSample = new DataSpace(1, size1, NULL);

	if(this->numGenes() > 0){
		size2[0] = this->numGenes();
		size2[1] = this->numSamples();
		spcExpData = new DataSpace(2, size2, NULL);
		spcRefSeq = new DataSpace(1, size2, NULL);
	}
	// create tables
	buildStruct();

	// gene expression
	DataSpace ds;
	if(this->numGenes() > 0){
		tblExpData = new DataSet(this->openDataSet(this->ExpDataPath()));
		tblRefSeq = new DataSet(this->openDataSet(this->RefSeqInfoPath()));
	} else {
		tblExpData = NULL;
		tblRefSeq = NULL;
		cout << "warning: no dataset with expression data." << endl;
	}

	// set tables
	tblSNPInfo = new DataSet(this->openDataSet(this->SNPInfoPath()));
	tblSample = new DataSet(this->openDataSet(this->SamplePath()));
	tblSNPData = new DataSet(this->openDataSet(this->SNPDataPath()));
	tblSNPDataInv = new DataSet(this->openDataSet(this->SNPDataInvPath()));
	
	// set linear memory dataspaces & buffers
	hsize_t	ms_size[2];
	// 1 x M
	ms_size[0] = this->numSamples();
	spaceMemNSamp = new DataSpace(1, ms_size, NULL);
	// 1 x N
	ms_size[0] = this->numSnps();
	spaceMemNSnp = new DataSpace(1, ms_size, NULL);
	// 1 x 1
	ms_size[0] = 1;
	spaceMem1 = new DataSpace(1, ms_size, NULL);
	

	// xfer
	#define SZ_HBUFF 1000
	tconv_buf = new char [SZ_HBUFF];
	xfer.setBuffer (SZ_HBUFF, tconv_buf, NULL);
}

template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::buildHdfTypes(){

	// Allocate local buffers
	buf_snpRow = new T_GType[this->numSamples()];
	buf_sampleRow = new T_GType[this->numSnps()];
	buf_expRow = new T_ExpType[this->numSamples()];

	// Create string types
	StrType hStringType24(PredType::C_S1, 24);
	StrType hStringType20(PredType::C_S1, 20);
	StrType hStringType12(PredType::C_S1, 12);
	StrType hStringType3(PredType::C_S1, 3);

	// Create TSample as a complex type
	TSampleType = new CompType(sizeof(T_Sample));

	// Add elements to the type
	// Use the offsetof macro to get the offsets into the type
	TSampleType->insertMember("index", HOFFSET(T_Sample, index), PredType::NATIVE_LONG);
	TSampleType->insertMember("name", HOFFSET(T_Sample, name), hStringType24);
	TSampleType->insertMember("famIndex1", HOFFSET(T_Sample, famIndex1), PredType::NATIVE_LONG);
	TSampleType->insertMember("famIndex2", HOFFSET(T_Sample, famIndex2), PredType::NATIVE_LONG);
	TSampleType->insertMember("score", HOFFSET(T_Sample, score), PredType::NATIVE_FLOAT);
	TSampleType->insertMember("phenotype", HOFFSET(T_Sample, phenotype), PredType::NATIVE_FLOAT);
	TSampleType->insertMember("sex", HOFFSET(T_Sample, sex), PredType::NATIVE_INT);

	TSNPType = new CompType(sizeof(T_SNP));
	TSNPType->insertMember("index", HOFFSET(T_SNP, index), PredType::NATIVE_LONG);
	TSNPType->insertMember("name", HOFFSET(T_SNP, name), hStringType24);
	TSNPType->insertMember("chr", HOFFSET(T_SNP, chr), hStringType3);
	TSNPType->insertMember("position", HOFFSET(T_SNP, position), PredType::NATIVE_LONG);
	TSNPType->insertMember("score", HOFFSET(T_SNP, score), PredType::NATIVE_FLOAT);
	TSNPType->insertMember("allele", HOFFSET(T_SNP, allele), hStringType3);

	TRefSeqType = new CompType(sizeof(T_RefSeq));
	TRefSeqType->insertMember("target", HOFFSET(T_RefSeq, target), hStringType20);
	TRefSeqType->insertMember("symbol", HOFFSET(T_RefSeq, symbol), hStringType20);
	TRefSeqType->insertMember("entrezID", HOFFSET(T_RefSeq, entrezID), PredType::NATIVE_LONG);
}

template<typename T_GType, typename T_ExpType>
hdfTPED<T_GType, T_ExpType>::~hdfTPED(void) {
	//spcExpData->close();
	//delete spcExpData;

	if(spcSample){
		spcSample->close();
		delete spcSample;
	}
	if(spcSNPData){
			spcSNPData->close();
			delete spcSNPData;
	}
	if(spcSNPDataInv){
		spcSNPDataInv->close();
		delete spcSNPDataInv;
	}
	if(spcSNPInfo){
		spcSNPInfo->close();
		delete spcSNPInfo;
	}
	if(tblSNPInfo){
		tblSNPInfo->close();
		delete tblSNPInfo;
	}
	if(tblSample){
		tblSample->close();
		delete tblSample;
	}
	if(tblSNPData){
		tblSNPData->close();
		delete tblSNPData;
	}
	if(tblSNPDataInv){
		tblSNPDataInv->close();
		delete tblSNPDataInv;
	}

	//tblExpData->close();
	//delete tblExpData;
	
	this->close();

	delete [] tconv_buf;
	H5garbage_collect();
}

template<typename T_GType, typename T_ExpType>
int hdfTPED<T_GType, T_ExpType>::doInvDataTable(){
	int bOk=true;
	int dims2[2];
	//T_GType **p2DArray = 0;

	// Create an inverted 2D array
	if(this->SamplesByRow()){
		dims2[1] = this->numSamples(); // #rows
		dims2[0] = this->numSnps(); // #cols
	} else {
		dims2[1] = this->numSnps(); // #rows
		dims2[0] = this->numSamples(); // #cols
	}

	//p2DArray = (T_GType **)avAlloc(2, dims2, AV_CHAR_TYPE, "SNPDataTableInv");
	//bOk = avSave(filename().c_str());
	
	T_GType val;
	for(int i=0;i<dims2[0];i++){
		for(int j=0;j<dims2[1];j++){
			this->getData(this->SNPDataPath(), j, i, val);
			this->setData(this->SNPDataInvPath(), i, j, val);
		}
		if((i%1000)==0) 
			cout << "#";
	}
	//bOk = avSave(filename().c_str());
	return bOk;
}

template<typename T_GType, typename T_ExpType>
bool hdfTPED<T_GType, T_ExpType>::doDataTable(){
    try {
		// Create a dataset using the default dataset creation properties.
		// We're not sure what they are, so we won't check.
		tblSNPData		= new DataSet (this->createDataSet(this->SNPDataPath(), NATIVE_GTYPE, *spcSNPData));
		tblSNPDataInv	= new DataSet (this->createDataSet(this->SNPDataInvPath(), NATIVE_GTYPE, *spcSNPDataInv));
		// Add a comment to the dataset
		this->setComment ("SNPDataTable", "This is the SNP dataset ordered by SNP-columns");
		this->setComment ("SNPDataTableInv", "This is the SNP dataset ordered by SNP-rows");
	}
	catch (Exception E) {
		// clean up and return with failure
		cout << "error: can't builb SNP data table." << endl;
		return false;
    }
	return true;
}

template<typename T_GType, typename T_ExpType>
bool hdfTPED<T_GType, T_ExpType>::doSampleTable(){
	try {
		// Create a dataset using the default dataset creation properties.
		// We're not sure what they are, so we won't check.
		tblSample = new DataSet(this->createDataSet(this->SamplePath(), *TSampleType, *spcSample));
		// Add a comment to the dataset
		this->setComment ("SampleTable", "This is the Sample dataset");
	}
	catch (Exception E) {
		// clean up and return with failure
		cout << "error: can't builb sample data table." << endl;
		return false;
    }
	return true;
}

template<typename T_GType, typename T_ExpType>
bool hdfTPED<T_GType, T_ExpType>::doSNPTable(){
	try {
		// Create a dataset using the default dataset creation properties.
		// We're not sure what they are, so we won't check.
		tblSNPInfo = new DataSet(this->createDataSet(this->SNPInfoPath(), *TSNPType, *spcSNPInfo));
		// Add a comment to the dataset
		this->setComment ("SNPInfoTable", "This is the SNP information dataset");
	}
	catch (Exception E) {
		// clean up and return with failure
		cout << "error: can't builb SNP info table." << endl;
		return false;
    }
	return true;
}

template<typename T_GType, typename T_ExpType>
bool hdfTPED<T_GType, T_ExpType>::doExpressionTable()
{
	try {
		// Create a dataset using the default dataset creation properties.
		// We're not sure what they are, so we won't check.
		tblExpData = new DataSet (this->createDataSet(this->ExpDataPath(), PredType::NATIVE_FLOAT, *spcExpData));
		tblExpData = new DataSet (this->createDataSet(this->RefSeqInfoPath(), *TRefSeqType, *spcRefSeq));
		// Add a comment to the dataset
		this->setComment ("ExpDataTable", "This is the gene expression dataset");
	}
	catch (Exception E) {
		// clean up and return with failure
		cout << "error: can't builb expression data table." << endl;
		return false;
    }
	return true;
}

// create a new file and its data structure
template<typename T_GType, typename T_ExpType>
bool hdfTPED<T_GType, T_ExpType>::buildStruct(){
	bool bOk = false;
	
	bOk = doDataTable();
	bOk = bOk && doSNPTable();
	bOk = bOk && doSampleTable();
	if(this->numGenes() > 0)
		bOk = bOk && doExpressionTable();
	return (bOk);
}

// public members

// setters
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::setData(const std::string &table, const unsigned row, const unsigned col, const T_ExpType val){
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];
	if(table == this->ExpDataPath()){
		ds_offset[0] = row;
		ds_offset[1] = col;
		ds_size[0] = 1;
		ds_size[1] = 1;
		spcExpData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblExpData->write(&val, NATIVE_ExpTYPE, *spaceMem1, *spcExpData, xfer);
	} else {
		throw 10;
	}
	return;
}

template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::setData(const std::string &table, const unsigned row, const unsigned col, const T_GType val){
	DataSpace	*dummySpc;
	DataSet		*dummySet;
	if(table == this->SNPDataPath()){
		dummySet=tblSNPData;
		dummySpc=spcSNPData;
	}
	else if(table == this->SNPDataInvPath()){
		dummySet=tblSNPDataInv;
		dummySpc=spcSNPDataInv;
	}
	else {
		throw 10;
	}
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];

	ds_offset[0] = row;
	ds_offset[1] = col;
	ds_size[0] = 1;
	ds_size[1] = 1;
	dummySpc->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
	dummySet->write(&val, NATIVE_GTYPE, *spaceMem1, *dummySpc, xfer);
}

template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::setData(const std::string &table, const unsigned row, const T_Sample &val){
	hsize_t	ds_size[1];
	hsize_t	ds_offset[1];
	ds_offset[0] = row;
	ds_size[0] = 1;
	spcSample->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
	tblSample->write(&val, *TSampleType, *spaceMem1, *spcSample, xfer);
}

template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::setData(const std::string &table, const unsigned row, const T_RefSeq &val){
	hsize_t	ds_size[1];
	hsize_t	ds_offset[1];

	ds_offset[0] = row;
	ds_size[0] = 1;
	spcRefSeq->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );

	tblRefSeq->write(&val, *TRefSeqType, *spaceMem1, *spcRefSeq, xfer);
}

template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::setData(const std::string &table, const unsigned row, const T_SNP &val){
	hsize_t	ds_size[1];
	hsize_t	ds_offset[1];
	ds_offset[0] = row;
	ds_size[0] = 1;
	spcSNPInfo->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
	tblSNPInfo->write(&val, *TSNPType, *spaceMem1, *spcSNPInfo, xfer);
}

//
// getters
//
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getData(const std::string &table, const unsigned row, const unsigned col, T_GType &val) const {
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];

	ds_offset[0] = row;
	ds_offset[1] = col;
	ds_size[0] = 1;
	ds_size[1] = 1;
	spcSNPData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
	tblSNPData->read(&val, NATIVE_GTYPE, *spaceMem1/*NSamp*/, *spcSNPData, xfer);
}

// carica tutti i genotipi dello SNP snpInd nel buffer buf_snpRow e fa puntare val a questo buffer
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getSnpPtr(const unsigned snpInd, T_GType *&val, const std::string &table) const {
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];

	if(table == this->SNPDataPath()){
		ds_offset[0] = 0;
		ds_offset[1] = snpInd;
		ds_size[0] = this->numSamples();
		ds_size[1] = 1;
		spcSNPData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblSNPData->read(buf_snpRow, NATIVE_GTYPE, *spaceMemNSamp, *spcSNPData, xfer);
	} else if(table == this->SNPDataInvPath()){
		ds_offset[1] = 0;
		ds_offset[0] = snpInd;
		ds_size[1] = this->numSamples();
		ds_size[0] = 1;
		spcSNPDataInv->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblSNPDataInv->read(buf_snpRow, NATIVE_GTYPE, *spaceMemNSamp, *spcSNPDataInv, xfer);
	} else {
		val = NULL;
		return;
	}
	val = buf_snpRow;
	return;
}

// carica tutti i genotipi dello SNP snpInd nel buffer buf_snpRow e fa puntare val a questo buffer
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getSnpSubsetMem(const unsigned snpInd, T_GType *val, const size_t mask_sz, const hsize_t *mask, const std::string &table) const {
	hsize_t	ds_size[2];
	hsize_t	ds1_size[1];
	hsize_t	ds_offset[2];

	if(table == this->SNPDataPath()){
		ds_offset[0] = 0;
		ds_offset[1] = snpInd;
		ds_size[0] = this->numSamples();
		ds_size[1] = 1;
		spcSNPData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblSNPData->read(val, NATIVE_GTYPE, *spaceMemNSamp, *spcSNPData, xfer);
	} else if(table == this->SNPDataInvPath()){
		ds_offset[1] = 0;
		ds_offset[0] = snpInd;
		ds_size[1] = this->numSamples();
		ds_size[0] = 1;
		spcSNPDataInv->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		ds1_size[0] = mask_sz;//this->numSamples();
		H5::DataSpace memspace(1, ds1_size, NULL);
		memspace.selectElements(H5S_SELECT_SET, mask_sz, mask);
		//tblSNPDataInv->read(val, NATIVE_GTYPE, *spaceMemNSamp, *spcSNPDataInv, xfer);
		tblSNPDataInv->read(val, NATIVE_GTYPE, memspace, *spcSNPDataInv, xfer);
	} else {
		val = NULL;
	}
	return;
}

// carica tutti i genotipi dello SNP snpInd nel buffer buf_snpRow e fa puntare val a questo buffer
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getSnpMem(const unsigned snpInd, T_GType *val, const std::string &table) const {
	assert(val != NULL);
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];

	if(table == this->SNPDataPath()){
		ds_offset[0] = 0;
		ds_offset[1] = snpInd;
		ds_size[0] = this->numSamples();
		ds_size[1] = 1;
		spcSNPData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblSNPData->read(val, NATIVE_GTYPE, *spaceMemNSamp, *spcSNPData, xfer);
	} else if(table == this->SNPDataInvPath()){
		ds_offset[1] = 0;
		ds_offset[0] = snpInd;
		ds_size[1] = this->numSamples();
		ds_size[0] = 1;
		spcSNPDataInv->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblSNPDataInv->read(val, NATIVE_GTYPE, *spaceMemNSamp, *spcSNPDataInv, xfer);
	} else {
		val = NULL;
	}
	return;	
}

template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getGxpPtr(const unsigned geneInd, T_ExpType *&val, const std::string &table) const{
assert(val != NULL);
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];

	if(table == this->ExpDataPath()){
		ds_offset[1] = 0;
		ds_offset[0] = geneInd;
		ds_size[1] = this->numSamples();
		ds_size[0] = 1;
		spcExpData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
	} else {
		val = NULL;
		return;
	}
	tblExpData->read(buf_expRow, NATIVE_ExpTYPE, *spaceMemNSamp, *spcExpData, xfer);
	val = buf_expRow;
}

// carica tutti i genotipi dell'individuo sampInd
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getSamplePtr(const unsigned sampInd, T_GType *&val, const std::string &table) const {
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];

	if(table == this->SNPDataPath()){
		ds_offset[0] = sampInd;
		ds_offset[1] = 0;
		ds_size[0] = 1;
		ds_size[1] = this->numSnps();
		spcSNPData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblSNPData->read(buf_sampleRow, NATIVE_GTYPE, *spaceMemNSnp, *spcSNPData, xfer);
		val = buf_sampleRow;
	} else {
		val = NULL;
	}
	return;
}

// carica tutti i genotipi dell'individuo sampInd
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getSampleMem(const unsigned sampInd, T_GType *val, const string &table) const {
	assert(val != NULL);
	hsize_t	ds_size[2];
	hsize_t	ds_offset[2];

	if(table == this->SNPDataPath()){
		ds_offset[0] = sampInd;
		ds_offset[1] = 0;
		ds_size[0] = 1;
		ds_size[1] = this->numSnps();
		spcSNPData->selectHyperslab( H5S_SELECT_SET, ds_size, ds_offset );
		tblSNPData->read(val, NATIVE_GTYPE, *spaceMemNSnp, *spcSNPData, xfer);
	}
}

// legge uno snp-info alla volta
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getData(const std::string &table, const unsigned row, T_SNP &val) const {
	assert(table == this->SNPInfoPath());
	const hsize_t	hs_size[1] = {1};
	hsize_t	hs_offset[1]; 
	hs_offset[0] = row;

	spcSNPInfo->selectHyperslab( H5S_SELECT_SET, hs_size, hs_offset );
	tblSNPInfo->read (&val, *TSNPType, *spaceMem1/*memSnp*/, *spcSNPInfo, xfer);
}

// legge un sample alla volta
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getData(const std::string &table, const unsigned row, T_Sample &val) const {
	assert(table == this->SamplePath());
	const hsize_t	hs_size[1] = {1};
	hsize_t			hs_offset[1];
	hs_offset[0] = row;

	spcSample->selectHyperslab( H5S_SELECT_SET, hs_size, hs_offset );
	tblSample->read(&val, *TSampleType, *spaceMem1/*memSample*/, *spcSample, xfer);
}

// legge un refseq alla volta
template<typename T_GType, typename T_ExpType>
void hdfTPED<T_GType, T_ExpType>::getData(const std::string &table, const unsigned row, T_RefSeq &val) const {
	assert(table == this->RefSeqInfoPath());
	const hsize_t	hs_size[1] = {1};
	hsize_t			hs_offset[1];
	hs_offset[0] = row;

	spcRefSeq->selectHyperslab( H5S_SELECT_SET, hs_size, hs_offset );
	tblRefSeq->read(&val, *TRefSeqType, *spaceMem1, *spcRefSeq, xfer);
}

template class hdfTPED<char, float>;