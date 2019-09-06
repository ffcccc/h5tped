#ifndef HDFVAR_I_CLASS_H
#define HDFVAR_I_CLASS_H

#pragma once
#include <string.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

#define NOT_IMPLEMENTED 0

#define T_ID 0
#define T_NAME 1
#define T_PHENO 2
#define S_NAME 0
#define S_CHR 1
#define S_POS 2
#define S_ALLELE 2
#define S_SCORE 4

#ifndef H5_NO_NAMESPACE
	#ifndef H5_NO_STD
		using std::cerr;
		using std::endl;
	#endif
#endif

//#include "testhdf5.h"	// C test header file
//int test();
/*
 * Include required headers.  This file tests internal library functions,
 * so we include the private headers here.
 */
#include "hdf5.h"
#include "H5public.h"
//#include "H5private.h"
//#include "H5Eprivate.h"
#include "H5Cpp.h"	// C++ API header file
#include "H5PredType.h"
#define NATIVE_ExpTYPE PredType::NATIVE_FLOAT 
#define NATIVE_GTYPE PredType::NATIVE_CHAR 


using namespace std;

//---------------------------------
//class h5TPED_I;
//---------------------------------

//typedef void(h5TPED_I::*pSnpPicker)(const unsigned row, T_GType *val, const std::string &table) const;

template<typename _GType, typename _ExpType>
class h5TPED_I : public H5::H5File {
	
public:
	typedef _ExpType T_ExpType;
	//typedef short int T_GType;
	//typedef char T_GType;
	typedef _GType T_GType;

		typedef enum {ALL_USERDEF, ALL_AB, ALL_A_B, ALL_ACGT} TEnumAlleleCode;
		static double sex_threshold_male()	{return 0.8; };
		static double sex_threshold_female(){return 0.2; };
		static const int Male 	= 1;	// plink = 1
		static const int Female	= 0;	// plink = 2
		static const int Undef	= -1;	// plink = any other than {1,2}
	
private:
	#ifdef _DEBUG
		std::vector<unsigned> testVec;
		std::vector<unsigned> trainVec;
	#endif


	typedef std::map <std::string, T_GType> Tmap_GType;
	typedef std::map <T_GType, std::string> Tmap_Alls;
	typedef typename Tmap_Alls::const_iterator Tmap_ITER;

	typedef struct {
		char	name[8];
		int		nProtCodingGenes;
		int		nBases;
		int		haploid;
		std::vector<int>	snpData;
	} T_Chr;

	typedef std::map <std::string, T_Chr > TChrMap;
	typedef std::map<std::string, int>	M_StrIndex;
	typedef M_StrIndex::iterator		IT_M_StrIndex;
	typedef M_StrIndex::const_iterator	ITc_M_StrIndex;

	typedef void(h5TPED_I::*pSnpPicker)(const unsigned row, T_GType *val, const std::string &table) const;

		// transposed flag 
		// switch action:	true -> TPed-like	false -> Ped-like )
		bool m_SamplesInRows;
		// h5 file name
		std::string m_filename;
		//
		std::string m_NA;
		// number of chrs in the data
		TChrMap m_Chr;
		void loadChrInfo();
		// number of chrs present in dataset
		long m_nChr;
		// number of variables/snps
		long m_nSnp;
		// number of samples
		long m_nSamples;
		// number of gene expression measurements
		long m_nGenes;
		// list of phenotype/covariates available
		std::vector <std::string> m_pheno;
		
		// codifica degli alleli
		TEnumAlleleCode m_codifica;
		Tmap_GType m_alleles;
		Tmap_Alls m_gtypes;
		//std::map <T_GType, std::string> m_gtypes;

	protected:

		typedef struct {
			long	index;
			char	name[24];
			char	chr[3];
			long	position;
			float	score;
			char	allele[3];
		} T_SNP;

		typedef struct {
			long	index;
			char	name[24];
			long	famIndex1;
			long	famIndex2;
			float	score;
			float	phenotype;
			int		sex;
		} T_Sample;

		typedef struct {
			long	entrezID;
			char	symbol[16];
		} T_Gene;

		typedef struct {
			char	target[20];
			char	symbol[20];
			long	entrezID;
		} T_RefSeq;
		// file metadata



		// create a new file and its data structure
		virtual bool buildStruct()=0;
		// depndent build methods
		virtual bool doDataTable()=0;
		virtual bool doSampleTable()=0;
		virtual bool doSNPTable()=0;
		virtual bool doExpressionTable()=0;

		// setters
		virtual void setData(const std::string &table, const unsigned row, const unsigned col, const T_ExpType val)=0;
		virtual void setData(const std::string &table, const unsigned row, const unsigned col, const T_GType val)=0;
		virtual void setData(const std::string &table, const unsigned row, const T_Sample &val)	=0;
		virtual void setData(const std::string &table, const unsigned row, const T_SNP &val)		=0;
		virtual void setData(const std::string &table, const unsigned row, const T_RefSeq &val)	=0;
		//virtual void setData(const std::string &table, const unsigned row, const long &val)=0;
		// getters
		virtual void getData(const std::string &table, const unsigned row, const unsigned col, T_ExpType &val)const =0;
		virtual void getData(const std::string &table, const unsigned row, const unsigned col, T_GType &val)	const =0;
		virtual void getData(const std::string &table, const unsigned row, T_Sample &val)const =0;
		virtual void getData(const std::string &table, const unsigned row, T_SNP &val)	const =0;
		virtual void getData(const std::string &table, const unsigned row, T_RefSeq &val)const =0;
		//virtual void getData(const std::string &table, const unsigned row, long &val) const =0;

		// memory indexes
		M_StrIndex _sampleIndex;
		M_StrIndex _snpIndex;
		M_StrIndex _geneIndex;
		
		// function to build indexes
		virtual bool buildIndex(){
			bool result = true;
			if(_snpIndex.empty()){
				std::string dummyname("");
				for(unsigned i=0; i<numSnps(); i++){
					getSNPName(i, dummyname);
					// check duplicates here!
					_snpIndex[dummyname] = i;
				}
			}
			result = result && _snpIndex.size() == numSnps();
			if(_sampleIndex.empty()){
				std::string dummyname("");
				for(unsigned i=0; i<numSamples(); i++){
					getSampleName(i, dummyname);
					// check duplicates here!
					_sampleIndex[dummyname] = i;
				}
			}
			result = result && _sampleIndex.size() == numSamples();
			return result;
		};

	public:
		// Empty constructor
		//h5TPED_I();
		// Constructor from existing file
		//h5TPED_I(const std::string &szFilename);
		// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables
		//h5TPED_I(const std::string &szFilename, const int nsnp, const int nsamp, const bool SamplesInRows = true);
		virtual ~h5TPED_I(void) {
			m_filename = "";
		};
		//
		virtual void flush()=0;
		// val points to memory buffer in which SNP is loaded
		virtual void getSnpPtr(const unsigned row, T_GType *&val, const std::string &table = "/SNPDataTableInv") const = 0;
		// val is a pre-allocated buffer in which the SNP is loaded
		virtual void getSnpMem(const unsigned row, T_GType *val, const std::string &table = "/SNPDataTableInv") const = 0;
		//						
		virtual void getSamplePtr(const unsigned sampInd, T_GType *&val, const std::string &table = "/SNPDataTable") const = 0;
		//
		virtual void getSampleMem(const unsigned sampInd, T_GType *val, const std::string &table = "/SNPDataTable") const = 0;
		//
		virtual void getGxpPtr(const unsigned row, T_ExpType *&val, const std::string &table = "/ExpDataTable") const =0;
		//
		virtual int doInvDataTable()=0;
		//
		void getSampleSummary() const {
			cout << this->filename() << " contains " << numSamples() << " samples\n";
			//int *vcc = new int[numSamples()];
			//getSamplePheno(vcc);
			//valarray<int> phe(vcc, numSamples());
			//cout << phe.sum() << " cases, " << numSamples() - phe.sum()	<< " controls phenotypes\n";
		}
		//
		void getSNPSummary(const string &rs) {
			//cout << "Name\tChr\tNumSamples\tMAF\tA_Freq\tB_Freq\tAA_Freq\tAB_Freq\tBB_Freq\n";
			int snpIndex(-1);
			int nA(0), nB(0), nAA(0), nAB(0), nBB(0), nNA(0);
			double maf = 0.0;
			float N(1.0);
			snpIndex = lookupSNPRow(rs);
			if (snpIndex != -1) {
				getSNPMaf(snpIndex, maf, nA, nB, nAA, nAB, nBB, nNA);
				//getSNPFreq(snpIndex,nA,nB,nAA,nAB,nBB,nNA);
				N = (float) nAA + nAB + nBB;
			}
			cout << rs << "\t" << "--" << "\t" << N + nNA << "(" << nNA << ")\t" << maf
					<< "\t" << nA / N << "\t" << nB / N << "\t" << nAA / N << "\t"
					<< nAB / N << "\t" << nBB / N << "\n";
		}

		// General Info ------------------------------------------------------------------------------------------------------------------------
		std::string filename() const { return m_filename; };
		// #samples
		inline unsigned numSamples() const { return m_nSamples; };
		inline void setNumSamples(unsigned nsamps) {m_nSamples = nsamps; };
		inline bool SamplesByRow() const { return m_SamplesInRows; };
		inline void setSamplesByRow(bool byrow) { m_SamplesInRows = byrow; };
		// # snps
		inline unsigned numSnps() const { return m_nSnp; };
		inline void setNumSnps(unsigned nsnps) {m_nSnp = nsnps; };
		// chrs
		inline unsigned numChrs() const { return m_nChr; };
		inline bool chr_sex(int snpInd) const { std::string chr; getSNPChr(snpInd, chr); return chr=="X";};
		inline bool chr_haploid(int snpInd) const { std::string chr; getSNPChr(snpInd, chr); return chr=="Y" || chr=="MT";};
		// genes
		inline unsigned numGenes()const { return m_nGenes; };
		inline void setNumGenes(unsigned ngenes) { m_nGenes = ngenes; };
		
		// default value for NA data
		inline T_GType NA() const { return -1; }
		inline void setNA(const string &NA) { m_NA = NA; }
		
		int getFirstNotNullSNP() const { 
			std::string rs("");
			for(unsigned i=0; i<numSnps(); i++){
				this->getSNPName(i,rs);
				if(!rs.empty())
					return i;
			}
		};
		// Sample Info -------------------------------------------------------------------------------------------------------------------------
		// Sample row index
		int lookupSampleRow(const std::string & name, bool useIndex=true) /*const*/ {
			int result(-1);	
			bool indexOK = false;
			if(useIndex){
				if(_sampleIndex.empty())
					indexOK = buildIndex();
				else indexOK = true;

				if(indexOK){
					ITc_M_StrIndex _idxIter = _sampleIndex.find(name);
					if(_idxIter != _sampleIndex.end())
						result = _idxIter->second;
				}
			}

			if((!useIndex )|| (!indexOK)){
				T_Sample sam;
				for(unsigned i=0; i<numSamples(); i++){
					getData("/SampleTable", i, sam);
					if(name == sam.name){
						result = i;
						break;
					}
				}
			}
			return result;
		};
		// sampling population with a specific comb. of allele
		void getSamples(const unsigned snpIndex, const T_GType geno, std::vector<unsigned> &samps){
			samps.clear();
			T_GType *buff;
			getSnpPtr(snpIndex, buff);
			for(unsigned row=0;row<numSamples();row++){
				if(buff[row]==geno){
					samps.push_back(row);
				}
			}
		}

		//
		void getSampleId(const unsigned rowInd, int &id) const { T_Sample sam; getData("/SampleTable", rowInd, sam); id = sam.index;}
		void setSampleId(const unsigned rowInd, const int id) { T_Sample sam; getData("/SampleTable", rowInd, sam); sam.index = id; setData("/SampleTable", rowInd, sam);}
		// Sample Name/code
		void setSampleName(const unsigned rowInd, const std::string & name) {
			T_Sample sam;
			getData("/SampleTable", rowInd, sam);
			strcpy(sam.name, name.c_str());
			setData("/SampleTable", rowInd, sam);
		}
		void getSampleName(const unsigned rowInd, std::string &name) const { T_Sample sam; getData("/SampleTable", rowInd, sam); name = sam.name;}
		// Family id
		void setSampleFamily(const unsigned rowInd, const long name) {
			T_Sample sam;
			getData("/SampleTable", rowInd, sam);
			sam.famIndex1 = name;
			setData("/SampleTable", rowInd, sam);
		}
		void getSampleFamily(const unsigned rowInd, long &name) const { T_Sample sam; getData("/SampleTable", rowInd, sam); name = sam.famIndex1;}
		// Sex
		void setSampleSex(const unsigned rowInd, const int val) {
			T_Sample sam;
			getData("/SampleTable", rowInd, sam);
			sam.sex = val;
			setData("/SampleTable", rowInd, sam);
		}
		void getSampleSex(const unsigned rowInd, int &val) const { T_Sample sam; getData("/SampleTable", rowInd, sam); val = sam.sex;}
		// Freq.
		//void getSampleFreq(unsigned index, int &nAA, int &nAB, int &nBB, int &nNA) const;
		// Inbreeding
		//void getSampleInbreeding(int index, double &het, bool check_sex=false);
		// Phenotype
		void setSamplePheno(const unsigned rowInd, int phe) { T_Sample sam; getData("/SampleTable", rowInd, sam); sam.phenotype = (float)phe; setData("/SampleTable", rowInd, sam);}
		void getSamplePheno(const unsigned rowInd, int &phe) const { 
			T_Sample sam; getData("/SampleTable", rowInd, sam); 
			phe = (int)sam.phenotype;
		}
		//
		void getSamplePheno(int *phe) const { 
			for(unsigned i=0;i<numSamples();i++){ 
				getSamplePheno(i, phe[i]); 
			} 
		};
		//
		void setSamplePhenoQT(const unsigned rowInd, float phe) { T_Sample sam; getData("/SampleTable", rowInd, sam); sam.phenotype = phe; setData("/SampleTable", rowInd, sam);}
		void getSamplePhenoQT(const unsigned rowInd, float &phe) const { T_Sample sam; getData("/SampleTable", rowInd, sam); phe = sam.phenotype;}
		bool isCase(const unsigned rowInd) const { T_Sample sam; getData("/SampleTable", rowInd, sam); return (sam.phenotype == 1.0);}
		
		// Genotypes Info ------------------------------------------------------------------------------------------------------------------------------
		void setGType(const unsigned rowInd, int colInd, T_GType gt)
		{
			//setData("/SNPDataTable", rowInd, colInd, gt);
			setData("/SNPDataTableInv", colInd, rowInd, gt);
		}
		void getGType(const unsigned rowInd, int colInd, T_GType &gt) const { getData("/SNPDataTable", rowInd, colInd, gt);}

		unsigned int numGTypeCodes() const { return m_gtypes.size(); }
		const std::string & decodeGType(T_GType iAll) const {
			Tmap_ITER it = m_gtypes.find(iAll);
			if(it == m_gtypes.end()){
				return m_NA;
			}
			return it->second;
		}
		// Acceptedd values
		const bool validAllele(char all0) const {
			return	(all0 == 'A')||
					(all0 == 'B')||
					(all0 == 'C')||
					(all0 == 'G')||
					(all0 == 'T')||
					(all0 == '1')||
					(all0 == '2')||
					(all0 == '3')||
					(all0 == '4');
		}
		
		const T_GType decodeAllele(std::string strAll) const {
		typename Tmap_GType::const_iterator it = m_alleles.find(strAll);
		if(it == m_alleles.end()){
			return NA();
		}
		return it->second;
	}
	
		const T_GType decodeAllele(int index, char all0, char all1) {
			bool valid0			= validAllele(all0);
			bool valid1			= validAllele(all1);
			bool valid			= valid0 && valid1;
			bool omozygote		= valid && (all0 == all1);
			bool etherozygote	= valid && (all0 != all1);
			char currAll0, currAll1;
			T_GType result(-1);

			// update current coding
			getSNPAllele(index, currAll0, currAll1);
			bool h5valid0	= validAllele(currAll0);
			bool h5valid1	= validAllele(currAll1);
			bool h5valid	= h5valid0 && h5valid1;
							
			if(!valid){
				result = T_GType(-1);
				if(!h5valid){
					// -C A-
					// -- A-
					if(!h5valid0 && valid0){
						if(all0 != currAll0)
							setSNPAllele(index, 0, all0);
					}
					// -C -A
					// -- -A
					else if(!h5valid0 && valid1){
						if(all1 != currAll0)
							setSNPAllele(index, 0, all1);
					}
					// A- C-
					// -- C-
					else if(!h5valid1 && valid0){
						if(all0 != currAll0)
							setSNPAllele(index, 1, all0);
					}
					// A- -C
					// -- -C
					else if(!h5valid1 && valid1){
						if(all1 != currAll0)
							setSNPAllele(index, 1, all1);
					}
				}
			}
			else if(valid){
				if(etherozygote){
					// ** AC
					result = T_GType(1);
					// -- AC
					if(!h5valid0 && !h5valid1){
						setSNPAllele(index, 0, all0);
						setSNPAllele(index, 1, all1);
					}
					// -C AC
					else if(!h5valid0 && (all0 != currAll1))
						setSNPAllele(index, 0, all0);
					// -C CA
					else if(!h5valid0 && (all1 != currAll1))
						setSNPAllele(index, 0, all1);
					// A- CA
					else if(!h5valid1 && (all0 != currAll0))
						setSNPAllele(index, 1, all0);
					// A- AC
					else if(!h5valid1 && (all1 != currAll0))
						setSNPAllele(index, 1, all1);
				}
				else if(omozygote && h5valid){	//omozyg H5 OK
					// AC AA
					if(all0 == currAll0)
						result = T_GType(0);
					// AC CC
					else if(all0 == currAll1)
						result = T_GType(2);
				}
				else if(omozygote && !h5valid){	//omozyg H5 TO UPDATE
					// -- AA
					if(!h5valid0 && !h5valid1){
						result = T_GType(0);
						setSNPAllele(index, 0, all0);
					}
					// A- AA
					else if(h5valid0 && (all0 == currAll0)){
						result = T_GType(0);
						//setSNPAllele(index, 1, all0);
					}
					// A- CC
					else if(h5valid0 && (all0 != currAll0)){
						result = T_GType(2);
						setSNPAllele(index, 1, all0);
					}
					// -C CC
					else if(h5valid1 && (all0 == currAll1)){
						result = T_GType(2);
						//setSNPAllele(index, 0, all0);
					}
					// -C AA
					else if(h5valid1 && (all0 != currAll1)){
						result = T_GType(0);
						setSNPAllele(index, 0, all0);
					}
				}
			}
			return result;
		}

		TEnumAlleleCode codifica() const { return m_codifica;};
		void loadAllele(const TEnumAlleleCode val);
		void loadAllele(const std::vector <std::string> &allList, const TEnumAlleleCode codif=ALL_USERDEF){
			m_codifica = codif;
			for(unsigned int i=0;i<allList.size();i++){
				m_gtypes.insert(typename Tmap_Alls::value_type(i, allList[i]));
				m_alleles.insert(typename Tmap_GType::value_type(allList[i], i));
			}
		}
		void loadAllele(const std::vector <std::string> &allList, const std::vector <T_GType> &codeList, const TEnumAlleleCode codif=ALL_USERDEF){
			m_codifica = codif;
			for(unsigned int i=0;i<allList.size();i++){
				m_gtypes.insert(typename Tmap_Alls::value_type(codeList[i], allList[i]));
				m_alleles.insert(typename Tmap_GType::value_type(allList[i], codeList[i]));
			}
		}

		// Gene & RefSeq ------------------------------------------------------------------------------------------------------------------------
		// return the name of the refseq (~gene) corresponding to position rowInd in table RefSeqInfo. Beware index starting at 0 !!
		void getGxpSymbol(const unsigned rowInd, std::string & name) const { 
			T_RefSeq rfs;
			getData("/RefSeqInfoTable", rowInd, rfs);
			name = rfs.symbol;
		}
		void setGxpSymbol(const unsigned rowInd, const std::string & name) { T_RefSeq s; getData("/RefSeqInfoTable", rowInd, s); strcpy(s.symbol, name.c_str()); setData("/RefSeqInfoTable", rowInd, s);}
		//
		void getGxpId(const unsigned rowInd, std::string & name) const { 
			T_RefSeq rfs;
			getData("/RefSeqInfoTable", rowInd, rfs);
			name = rfs.target;
		}
		void setGxpId(const unsigned rowInd, const std::string & name) { T_RefSeq s; getData("/RefSeqInfoTable", rowInd, s); strcpy(s.target, name.c_str()); setData("/RefSeqInfoTable", rowInd, s);}
		//
		void setGExp(const unsigned rowInd, int colInd, T_ExpType et) { setData("/ExpDataTable", rowInd, colInd, et);}
		void getGExp(const unsigned rowInd, int colInd, T_ExpType &et) const { getData("/ExpDataTable", rowInd, colInd, et);}


		// SNP Info ------------------------------------------------------------------------------------------------------------------------------
		//
		void setSNPName(const unsigned rowInd, const std::string & name) { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); strcpy(snp.name, name.c_str()); setData("/SNPInfoTable", rowInd, snp);}
		// Chromosome
		void getSNPChr(const unsigned rowInd, std::string & name) { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); name = snp.chr;}
		void setSNPChr(const unsigned rowInd, const std::string & name) { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); strcpy(snp.chr, name.c_str()); setData("/SNPInfoTable", rowInd, snp);}
		// SNP Base Position
		void setSNPPos(const unsigned rowInd, int pos) { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); snp.position = pos; setData("/SNPInfoTable", rowInd, snp);}
		//void setSNPGene(const unsigned rowInd, long ind) { setData("/EntrezInfoTable", rowInd, ind);}
		// todo verificare la size del campo allele
		int getSNPNumAllele(const unsigned rowInd) const { T_SNP snp; int result(0); getData("/SNPInfoTable", rowInd, snp);
			if(snp.allele[0] != '-') result++;
			if(snp.allele[1] != '-') result++;
			return result;
		};
		void setSNPAllele(const unsigned rowInd, const std::string & name) { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); strcpy(snp.allele, name.c_str()); setData("/SNPInfoTable", rowInd, snp);}
		void setSNPAllele(const unsigned rowInd, int allInd, char all)		{ T_SNP snp; getData("/SNPInfoTable", rowInd, snp); snp.allele[allInd] = all; setData("/SNPInfoTable", rowInd, snp);}
		//
		void setSNPScore(const unsigned rowInd, float score) { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); snp.score = score; setData("/SNPInfoTable", rowInd, snp);}
		// return the name of the SNP (eg "rs1234") corresponding to position rowInd. Beware index starting at 0 !!
		void getSNPName(const unsigned rowInd, std::string & name) const { 
			T_SNP snp; 
			getData("/SNPInfoTable", rowInd, snp); 
			name = snp.name;
		}
		//
		virtual void getSNPChr(const unsigned rowInd, std::string & name) const { 
			T_SNP snp; 
			getData("/SNPInfoTable", rowInd, snp); 
			name = snp.chr;}
		// set & get an index 'associated' with the rowInd-th SNP (e.g. ncbi gene-id...)
		void setSNPIndex(const unsigned rowInd, int ent_ind) { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); snp.index = ent_ind; setData("/SNPInfoTable", rowInd, snp);}
		void getSNPIndex(const unsigned rowInd, int &ent_ind) const { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); ent_ind = snp.index;};
		int lookupSNPRow(const std::string &rs, bool useIndex=true) {
			int result(-1);	
			bool indexOK = false;
			if(useIndex){
				if(_snpIndex.empty())
					indexOK = buildIndex();
				else indexOK = true;
				
				if(indexOK){
					ITc_M_StrIndex _idxIter = _snpIndex.find(rs);
					if(_idxIter != _snpIndex.end()){
						result = _idxIter->second;
					}
				}
			}
			if((!useIndex )|| (!indexOK)){
				T_SNP snp;
				for(unsigned i=0; i<numSnps(); i++){
					getData("/SNPInfoTable", i, snp);
					if(rs == snp.name){
						result = i;
						break;
					}
				}
			}
			return result;
		}
		//
		//void getSNPSummary(std::string rs);
		//
		void getSNPPos(const unsigned rowInd, long &pos) const { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); pos = snp.position;}
		//
		void getSNPAllele(const unsigned rowInd, std::string & name) const { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); name = snp.allele; }
		void getSNPAllele(const unsigned rowInd, char &all0, char &all1) const { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); all0 = snp.allele[0]; all1 = snp.allele[1];}
		//
		void getSNPScore(const unsigned rowInd, float &score) const { T_SNP snp; getData("/SNPInfoTable", rowInd, snp); score = snp.score;}
		//
		//void getSNPFreq(unsigned rsindex, int &nA, int &nB, int &nAA, int &nAB, int &nBB, int &nNA) const;
		//
		void getSNPMaf(unsigned rsindex, double &maf) const {
			int nA,nB,nAA,nAB,nBB,nNA;
			getSNPMaf(rsindex, maf, nA,nB,nAA,nAB,nBB,nNA);
		}
		void getSNPMaf(unsigned rsindex, double &maf, int &nA, int &nB, int &nAA, int &nAB, int &nBB, int &nNA) const {
			getSNPFreq(rsindex,nA,nB,nAA,nAB,nBB,nNA);
			if(nA<nB) maf = (double)nA / (double)(nA+nB);
			else	  maf = (double)nB / (double)(nA+nB);
		};
		void getSNPHWE(double &p_hwe, int obs_hets, int obs_hom1, int obs_hom2);



	//-----------------------------------------------------------------------------------------------------------
	//                                   TPED Interface class h5TPED_I
	//-----------------------------------------------------------------------------------------------------------

	// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables
	h5TPED_I();
// {
//		m_filename = "";
//		m_nSnp = 0; // table 1xN
//		m_Chr.clear();
//		m_nSamples = 0;// table 1xM
//		m_SamplesInRows = true;
//		m_NA = "NA";
//	}

	//
	h5TPED_I(const std::string &szFilename);
// : H5::H5File(szFilename, H5F_ACC_RDWR) {
//		m_filename = szFilename;
//		m_nSnp = 0; // table 1xN
//		m_Chr.clear();
//		m_nSamples = 0;// table 1xM
//		m_SamplesInRows = true;
//		m_NA = "NA";
//	}

	// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables
	h5TPED_I(const std::string &szFilename, const int nsnp,
			const int nsamp, const bool szSampInRows);
//  : H5::H5File(szFilename, H5F_ACC_RDWR) {
//
//		m_filename = szFilename;
//		m_nSnp = nsnp;
//		m_Chr.clear();
//		m_nSamples = nsamp;
//		m_SamplesInRows = szSampInRows;
//		m_NA = "NA";
//		//buildStruct();
//	}

	// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables

	void getSampleFreq(unsigned index, int &nAA, int &nAB, int &nBB, int &nNA) const {
		T_GType *gVal;
		getSamplePtr(index, gVal);
		nAA = 0;
		nAB = 0;
		nBB = 0;
		nNA = 0;
		for (unsigned i = 0; i < numSnps(); i++) {
			switch (gVal[i]) {
			case -1:
				nNA++;
			case 0:
				nAA++;
				break;
			case 1:
				nAB++;
				break;
			case 2:
				nBB++;
				break;
			}
		}
		//free(gVal);
	}

	void getSNPFreq(unsigned rsindex, int &nA, int &nB, int &nAA, int &nAB, int &nBB, int &nNA) const {
		T_GType *gVal;
		getSnpPtr(rsindex, gVal);
		
		bool haplo(false), isX(false), haploF(false), inc1(false);
		if(chr_haploid(rsindex)){
			inc1 = true;
		} else if(chr_sex(rsindex)){
			isX = true;
		};
		int sx;
		nA = 0;
		nB = 0;
		nAA = 0;
		nAB = 0;
		nBB = 0;
		nNA = 0;
		for (unsigned i = 0; i < numSamples(); i++) {
			haploF = false;
			if(haplo){
				getSampleSex(i, sx);
				// set if Y Female
				haploF= (sx == h5TPED_I::Female);		
			}
			if(isX){
				getSampleSex(i, sx);
				// set if X Male
				inc1 = (sx == h5TPED_I::Male);
			}

			switch (gVal[i]) {
			case -1:
				if(haploF)
					// do not count NA in Y-Female samples
					break;
				nNA++;
				break;
			case 0:
				if(haplo)
					nA += 1;
				else
					nA += 2;
				nAA++;
				break;
			case 1:
				if(haplo)
					cout << "log: haploid heterozigote at SNP:" << rsindex << " and Sample:" << i << endl;
				nA++;
				nB++;
				nAB++;
				break;
			case 2:
				if(haplo)
					nA += 1;
				else
					nB += 2;
				nBB++;
				break;
			}
		}
	}

};
#endif
