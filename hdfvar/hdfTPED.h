#ifndef HDFTPED_CLASS_H
#define HDFTPED_CLASS_H

#pragma once

#include "hdfvar_I.h"

template<typename _GType, typename _ExpType>
class hdfTPED : public h5TPED_I<_GType, _ExpType> {

public:
	typedef typename h5TPED_I<_GType, _ExpType>::T_ExpType T_ExpType;
	typedef typename h5TPED_I<_GType, _ExpType>::T_GType T_GType;
	typedef typename h5TPED_I<_GType, _ExpType>::T_Sample T_Sample;
	typedef typename h5TPED_I<_GType, _ExpType>::T_SNP T_SNP;
	typedef typename h5TPED_I<_GType, _ExpType>::T_RefSeq T_RefSeq;

private:
		// datasets path
		inline string SNPInfoPath()	const { return("/SNPInfoTable"); };
		inline string SamplePath()	const { return("/SampleTable"); };
		inline string SNPDataPath()	const { return("/SNPDataTable"); };
		inline string ExpDataPath()	const { return("/ExpDataTable"); };
		inline string RefSeqInfoPath()	const { return("/RefSeqInfoTable"); };
		inline string SNPDataInvPath()	const { return("/SNPDataTableInv"); };
		// datasets
		H5::DataSet *tblSNPInfo;
		H5::DataSet *tblSample;
		H5::DataSet *tblSNPData;
		H5::DataSet *tblSNPDataInv;
		H5::DataSet *tblRefSeq;
		H5::DataSet *tblExpData;
		// dataspaces
		H5::DataSpace *spcSample;//(2, size, NULL);
		H5::DataSpace *spcSNPInfo;
		H5::DataSpace *spcSNPData;
		H5::DataSpace *spcSNPDataInv;
		H5::DataSpace *spcExpData;
		H5::DataSpace *spcRefSeq;
		// dataypest
		H5::CompType *TSampleType;
		H5::CompType *TSNPType;
		H5::CompType *TRefSeqType;
		// memory-spaces
		H5::DataSpace *spaceMemNSamp;
		H5::DataSpace *spaceMemNSnp;
		H5::DataSpace *spaceMem1;
		// transfer buffer
		H5::DSetMemXferPropList xfer;
		char* tconv_buf;// = new char [1000];
		// global buffers
		T_GType	*buf_snpRow;
		T_GType	*buf_sampleRow;
		T_ExpType *buf_expRow;


	bool _SnpByRow;

	virtual bool doDataTable();
	virtual bool doSampleTable();
	virtual bool doSNPTable();
	virtual bool doExpressionTable();
	// create a new file and its data structure
	virtual bool buildStruct();
	// creates hdf complex types
	void buildHdfTypes();
	// setters
	virtual void setData(const std::string &table, const unsigned row, const unsigned col, const T_ExpType val);
	virtual void setData(const std::string &table, const unsigned row, const unsigned col, const T_GType val);
	virtual void setData(const std::string &table, const unsigned row, const T_Sample &val);
	virtual void setData(const std::string &table, const unsigned row, const T_SNP &val);
	virtual void setData(const std::string &table, const unsigned row, const T_RefSeq &val);
	//virtual void setData(const std::string &table, const unsigned row, const long &val);
	// getters
	virtual void getData(const std::string &table, const unsigned row, const unsigned col, T_ExpType &val)const{};
	virtual void getData(const std::string &table, const unsigned row, const unsigned col, T_GType &val)	const;
	virtual void getData(const std::string &table, const unsigned row, T_Sample &val)const ;
	virtual void getData(const std::string &table, const unsigned row, T_SNP &val)	const ;
	virtual void getData(const std::string &table, const unsigned row, T_RefSeq &val)const ;
	//virtual void getData(const std::string &table, const unsigned row, long &val) const ;

public:

	// Empty constructor
	hdfTPED() : h5TPED_I<T_GType, T_ExpType>() {
		buildHdfTypes();
	};
	
	// Constructor from existing file
	hdfTPED(const std::string &szFilename);
	
	// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables
	hdfTPED(const std::string &szFilename, const int nsnp, const int nsamp, const int nexpr = 0, const bool SamplesInRows = true);
	
	// Destructor
	~hdfTPED(void);
	
	void flush() {
		tblSNPInfo->flush(H5F_SCOPE_LOCAL);
		tblSample->flush(H5F_SCOPE_LOCAL);
		tblSNPData->flush(H5F_SCOPE_LOCAL);
		if(tblSNPDataInv != NULL)	
			tblSNPDataInv->flush(H5F_SCOPE_LOCAL);
		if(tblRefSeq != NULL)	
			tblRefSeq->flush(H5F_SCOPE_LOCAL);
		if(tblExpData != NULL)	
			tblExpData->flush(H5F_SCOPE_LOCAL);
	};
	bool objExists(const string &path) const;
	//
	virtual void getGxpPtr(const unsigned row, T_ExpType *&val, const std::string &table = "/ExpDataTable") const;
	//
	virtual void getSnpPtr(const unsigned row, T_GType *&val, const std::string &table = "/SNPDataTableInv") const;
	virtual void getSnpSubsetMem(const unsigned row, T_GType *val, const size_t mask_sz, const hsize_t *mask, const std::string &table = "/SNPDataTableInv") const;
	//
	virtual void getSnpMem(const unsigned row, T_GType *val, const std::string &table = "/SNPDataTableInv") const;
	//
	virtual void getSamplePtr(const unsigned sampInd, T_GType *&val, const std::string &table = "/SNPDataTable") const;
	virtual void getSampleMem(const unsigned sampInd, T_GType *val, const std::string &table = "/SNPDataTable") const;
				 
	//
	virtual int doInvDataTable();
};

#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
