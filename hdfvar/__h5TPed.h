#ifndef H5TPED_CLASS_H
#define H5TPED_CLASS_H

#pragma once
#include "h5TPed_I.h"

class h5TPED : public h5TPED_I {
protected:
	virtual int doDataTable();
	virtual int doSampleTable();
	virtual int doSNPTable();
	// create a new file and its data structure
	virtual bool buildStruct();
	// setters
	virtual void setData(const std::string &table, const int row, const int col, const T_GType val);
	virtual void setData(const std::string &table, const int row, const T_Sample &val);
	virtual void setData(const std::string &table, const int row, const T_SNP &val);
	virtual void setData(const std::string &table, const int row, const long &val);
	// getters
	virtual void getData(const std::string &table, const int row, const int col, T_GType &val)	const ;
	virtual void getData(const std::string &table, const int row, T_Sample &val)	const ;
	virtual void getData(const std::string &table, const int row, T_SNP &val)	const ;
	virtual void getData(const std::string &table, const int row, long &val) const ;

public:

	// Empty constructor
	h5TPED() : h5TPED_I() {};
	// Constructor from existing file
	h5TPED(const std::string &szFilename);
	// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables
	h5TPED(const std::string &szFilename, const int nsnp, const int nsamp, const bool SamplesInRows = true) : h5TPED_I(szFilename, nsnp, nsamp, SamplesInRows) {
		buildStruct();
	};
	~h5TPED(void);

	virtual void getSnpRowP(const std::string &table, const int row, T_GType *&val) const;
	//void getDataRow(const std::string &table, const int row, std::valarray<T_GType> &val) const;
	virtual void getSnpRow(const std::string &table,const int row, std::valarray<T_GType> &valp) const;
	//void getDataCol(const std::string &table, const int col, std::valarray<T_GType> &val) const;
	virtual void getSnpCol(const std::string &table, const int col, std::valarray<T_GType> &valp) const;
	//
	virtual int doInvDataTable();
	//
	//virtual int doTestTrainTable();
};

#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
