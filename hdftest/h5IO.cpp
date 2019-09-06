#ifdef _MSC_VER
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#endif

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <set>
#include "hugeDataFile.h"
#include "../hdfvar/hdfTPED.h"
//#include <mysql++.h>
//#include "bio/biodata.h"

using namespace std;

// read from illumina fulldata table by column... very inefficient with IRFile object
int fullData2h5_bycol(IRFile &df, h5TPED_I &outFile)
{
	int nsamp = 0;
	int nsnp = 0;
	
	std::string strVal("");
	int strPhe;
	int dVal(100);
	for (unsigned j=0; j<outFile.numSamples(); j++) {
		std::string strVal("");
		outFile.setSampleId(j,j);

		df.GetVariableName(5+j*2, strVal);
		outFile.setSampleName(j,strVal);

		// pheno
		strPhe = 0; // case
		if(	(strVal.compare(0,5,"MES M")==0)  || (strVal.compare(0,6,"MES CA")==0) || (strVal.compare(0,5,"MES 1")==0)){
		//if(	(strVal.compare(0,5,"MES C")==0) ||	(strVal.compare(0,3,"CTR")==0))
		//if(	(strVal.compare(0,2,"CO")==0) ){
			strPhe = 1; // control
		}
		outFile.setSamplePheno(j,strPhe);
		
		
		// geno data
		for(unsigned i=0; i<	outFile.numSnps(); i++){
			df.GetData(5+j*2, i, strVal);
			outFile.setGType(i, j, outFile.decodeAllele(strVal));
			
			// map data
			if(j==0){
				df.GetData(1, i, strVal);
				outFile.setSNPName(i, strVal);
				df.GetData(2, i, strVal);
				outFile.setSNPChr(i, strVal);
				df.GetData(3, i, strVal);
				outFile.setSNPPos(i, atoi(strVal.c_str()));
				df.GetData(5+j*2+1, i, strVal);
				outFile.setSNPAllele(i, strVal);
				df.GetData(4, i, strVal);
				outFile.setSNPScore(i, atof(strVal.c_str()));
			}
		}
	}
	return 0;
}

int simplegraph2dot(CDataFile &df, bool recodeName) {
	int limit = df.GetNumberOfSamples(0);
	int i;
	h5TPED_I *h5=NULL;

	if(recodeName){
		h5 = new hdfTPED("tbqc2.h5");
	}
	
	for(i=0;i<limit;i++){
		string val0("");
		df.GetData(0, i, val0);
		if(recodeName){
			h5->getSNPName(atoi(val0.c_str())-1, val0);
		}
		val0 += " --";
		df.SetData(0, i, val0.c_str());

		df.GetData(1, i, val0);
		if(recodeName){
			h5->getSNPName(atoi(val0.c_str())-1, val0);
			df.SetData(1, i, val0.c_str());
		}

		string valw("");
		df.GetData(2, i, valw);
		valw = "[label=\""+valw+"\"]";
		df.SetData(2, i, valw.c_str());
	}
	string newname = df.m_szFilename + ".dot";
	if(recodeName){
		cout << df.m_szFilename << "Recoded SNP index with rs labels";
	}
	cout << "Converted " << df.m_szFilename << " to .dot format: " << newname;
	df.WriteFile(newname.c_str(), ' ');
	return i;
}

// reads from Illumina fulldata table by rows with much less disk I/O.
typedef struct {
	// fixed fields
	int fixedCols;
	int SNPName;
	int SNPChr;
	int SNPPos;
	int SNPScore;
	// repeated fields
	int repeatCols;
	int SNPAllele;
	int GType;
	int Score;
	int LogR;
	std::set<string> excluded1;
} FullDataPos;
	// impostazioni x meso500 full (> 4GB)
	//FullDataPos tOffset;
	//tOffset.fixedCols=6;
	//tOffset.GType=0;
	//tOffset.LogR=3;
	//tOffset.repeatCols=4;
	//tOffset.Score=1;
	//tOffset.SNPAllele=2;
	//tOffset.SNPChr=3;
	//tOffset.SNPName=1;
	//tOffset.SNPPos=4;
	//tOffset.SNPScore=5;

int fullData2h5(IRFile &df, h5TPED_I &outFile, unsigned totVars)
{
	int nsamp = 0;
	int nsnp = 0;
	// impostazioni x meso500 solo genotipi
	FullDataPos tOffset;
	tOffset.fixedCols=5;
	tOffset.SNPAllele=-1;
	tOffset.SNPChr=2;
	tOffset.SNPName=1;
	tOffset.SNPPos=3;
	tOffset.SNPScore=4;

	tOffset.repeatCols=1;
	tOffset.GType=0;
	tOffset.LogR=-1;
	tOffset.Score=-1;
	tOffset.excluded1.insert("MES 138.GType");
	tOffset.excluded1.insert("CTR 208.GType");
	tOffset.excluded1.insert("MES C 169.GType");
	tOffset.excluded1.insert("CTR 202.GType");
	tOffset.excluded1.insert("CTR 201.GType");
	tOffset.excluded1.insert("MES M 42.GType");
	tOffset.excluded1.insert("MES M 36.GType");
	tOffset.excluded1.insert("MES M 62.GType");
	tOffset.excluded1.insert("MES M 61.GType");
	tOffset.excluded1.insert("MES M 24.GType");
	tOffset.excluded1.insert("MES M 123.GType");
	tOffset.excluded1.insert("CTR 247.GType");	
	std::string strVal("");
	int strPhe,cont=0;
	int dVal(100);

	// geno data
	set<unsigned> excludedCols;
	for(unsigned i=0; i<	outFile.numSnps(); i++){
		
		if(i==0){
			// File Header -> Sample data j,0
			unsigned jInsert = 0;
			for (unsigned j=0; j<totVars/*outFile.numSamples()*/; j++) {
				std::string strVal("");
				
				df.GetVariableName(tOffset.fixedCols + j * tOffset.repeatCols, strVal);
				if(tOffset.excluded1.find(strVal)==tOffset.excluded1.end()){
					// ok...
					outFile.setSampleId(jInsert,jInsert);				
					outFile.setSampleName(jInsert,strVal);

					// pheno
					strPhe = 0; // case
					if(	(strVal.compare(0,5,"MES M")==0)  || (strVal.compare(0,6,"MES CA")==0) ||
						(strVal.compare(0,5,"MES 1")==0)  || (strVal.compare(0,6,"Mes TO")==0)){
						strPhe = 1; // control
					}
					outFile.setSamplePheno(jInsert,strPhe);
					jInsert++;
				} else {
					// individuo j da escludere
					tOffset.excluded1.erase(tOffset.excluded1.find(strVal));
					excludedCols.insert(j);
				}
			}
		}
		
		//Genotype i,j
		for (unsigned j=0; j<outFile.numSamples(); j++) {
			unsigned jInsert = 0;
			if(excludedCols.find(j)==excludedCols.end()){
				unsigned jInsert = 0;
				df.GetData(tOffset.fixedCols + j * tOffset.repeatCols + tOffset.GType, i, strVal);
				outFile.setGType(jInsert, i, outFile.decodeAllele(strVal));
			
				// map data
				if(jInsert==0){
					if(tOffset.SNPName >= 0){
						df.GetData(tOffset.SNPName, i, strVal);
						outFile.setSNPName(i, strVal);
					}
					if(tOffset.SNPChr >= 0){
						df.GetData(tOffset.SNPChr, i, strVal);
						outFile.setSNPChr(i, strVal);
					}
					if(tOffset.SNPPos >= 0){
						df.GetData(tOffset.SNPPos, i, strVal);
						outFile.setSNPPos(i, atoi(strVal.c_str()));
					}
					if(tOffset.SNPAllele >= 0){
						df.GetData(tOffset.fixedCols + j * tOffset.repeatCols, i, strVal);
						outFile.setSNPAllele(i, strVal);
					}
					if(tOffset.SNPScore >= 0){
						df.GetData(tOffset.SNPScore, i, strVal);
						outFile.setSNPScore(i, atof(strVal.c_str()));
					}
				}
			}
		}
		if(cont==10000){
			cont=0;
			cout << "#";
		} else cont++;
	}
	return 1;
}

// reads from Illumina fulldata table by rows with much less disk I/O.
int csv2h5(IRFile &df, h5TPED_I &outFile)
{
	int nsamp = 0;
	int nsnp = 0;
	
	std::string strVal(""), strChr(""), strGene("");
	int phe,bp,cont=0;
	int dVal(100);
	char buff[255];

	// Sample data j,0
	for (unsigned j=0; j<outFile.numSamples(); j++) {
		std::string strVal("item");
		outFile.setSampleId(j,j);
		outFile.setSampleName(j,strVal+_itoa(j,buff,10));
	}

	// geno data
	for (unsigned j=0; j<outFile.numSamples(); j++) {

		for(unsigned i=0; i<	outFile.numSnps(); i++){
			// map data
			if(j==0){
				df.GetVariableName(i, strVal);
				outFile.setSNPName(i, strVal);
				strChr=""; strGene=""; bp=0;
				//lookupSNP(strVal, strChr, bp);
				outFile.setSNPChr(i, strChr);
				outFile.setSNPPos(i, bp);
			/*	df.GetData(5+j*2+1, i, strVal);
				outFile.setSNPAllele(i, strVal);
				df.GetData(4, i, strVal);
				outFile.setSNPScore(i, atof(strVal.c_str()));*/
			}
			df.GetData(i, j, strVal);
			outFile.setGType(j, i, outFile.decodeAllele(strVal));
		}

		// sample phenotype
		df.GetData(outFile.numSnps(), j,strVal);
		phe = atoi(strVal.c_str());
		outFile.setSamplePheno(j, phe);

		if(cont==10000){
			cont=0;
			cout << "#";
		} else cont++;
	}
	return 1;
}

extern bool testInversion(h5TPED_I &h5);

// reads from infinium 15K
int infinium2h5(IRFile &df, unsigned nf1, IRFile &df2, unsigned nf2, h5TPED_I &outFile)
{
	std::string strVal(""), strChr(""), strGene("");
	int bp;
	vector<string> alls;
	alls.push_back("AA");	// -> 0
	alls.push_back("AC");	// -> 1
	alls.push_back("AG");	// -> 2
	alls.push_back("AT");	// -> 3
	alls.push_back("CA");	// -> 1
	alls.push_back("CC");	// -> 5
	alls.push_back("CG");	// -> 6
	alls.push_back("CT");	// -> 7
	alls.push_back("GA");	// -> 2
	alls.push_back("GC");	// -> 6
	alls.push_back("GG");	// -> 10
	alls.push_back("GT");	// -> 11
	alls.push_back("TA");	// -> 3
	alls.push_back("TC");	// -> 7
	alls.push_back("TG");	// -> 11
	alls.push_back("TT");	// -> 15
	vector <T_GType> codes;
	codes.push_back(0);	// -> 0
	codes.push_back(1);	// -> 1
	codes.push_back(1);	// -> 2
	codes.push_back(1);	// -> 3
	codes.push_back(1);	// -> 4
	codes.push_back(5);	// -> 5
	codes.push_back(1);	// -> 6
	codes.push_back(1);	// -> 7
	codes.push_back(1);	// -> 8
	codes.push_back(1);	// -> 9
	codes.push_back(10);	// -> 10
	codes.push_back(1);	// -> 11
	codes.push_back(1);	// -> 12
	codes.push_back(1);	// -> 13
	codes.push_back(1);	// -> 14
	codes.push_back(15);	// -> 15
	outFile.loadAllele(alls, codes);
	bool dbok(false);
	const char* db = "gene_info", *server = "localhost", *user = "root", *pass = "";
	// Connect to the sample database.
	mysqlpp::Connection conn(false);
	if (conn.connect(db, server, user, pass)) {
		dbok = true;
	}
	else {
		cerr << "DB connection failed: " << conn.error() << endl;
	}
	set<string> excluded1;
	set<string> excludedSNP;
	IRFile dsnp("Infinium_BC_excluded_snp.txt", DF::RF_READ_AS_MIXED, ' ', "NA", 2859, 0);
	for (int j=0; j<dsnp.GetNumberOfSamples(0); j++) {
		dsnp.GetData(0,j,strVal);
		excludedSNP.insert(strVal);
	}
	// MS
	// annotated
	//Samples with low sequenom calls:
	//NONE
	//Samples with call rates < 0.93:
	//WTCCC92740   0.8676
	//WTCCC92826   0.8905  
	//WTCCC93289   0.9059 
	//WTCCC93131   0.9162 
	//WTCCC93336   0.9236
	//WTCCC93130   0.9283
	//Suspected related samples:
	//NONE
	// #Exclude putatively related individuals
	excluded1.insert("WTCCC93211");
	excluded1.insert("WTCCC92862");
	excluded1.insert("WTCCC93384");
	excluded1.insert("WTCCC93283");
	// Exclude individuals with questionable ancestry
	excluded1.insert("WTCCC92958");
	excluded1.insert("WTCCC92905");
	excluded1.insert("WTCCC92997");
	excluded1.insert("WTCCC93217");
	// Exclude individuals who have >10% missing genotypes
	excluded1.insert("WTCCC92740");
	excluded1.insert("WTCCC92826");
	excluded1.insert("WTCCC93023");
	excluded1.insert("WTCCC93336");
	excluded1.insert("WTCCC92603");
	excluded1.insert("WTCCC92995");
	excluded1.insert("WTCCC93349");
	excluded1.insert("WTCCC92836");
	excluded1.insert("WTCCC93289");
	excluded1.insert("WTCCC93130");
	excluded1.insert("WTCCC93131");
	// 58C annotated
	// call rates < 0.93:
	//excluded1.insert('WTCCC66386');
	// low sequenom calls:
	//excluded1.insert('WTCCC66283');
	// BC
	// #Exclude putatively related individuals
	excluded1.insert("WTCCC93744");
	excluded1.insert("WTCCC93793 ");
	excluded1.insert("WTCCC93709 ");
	excluded1.insert("WTCCC93714 ");
	excluded1.insert("WTCCC93889 ");
	excluded1.insert("WTCCC94065 ");
	excluded1.insert("WTCCC93954 ");
	excluded1.insert("WTCCC94149 ");
	excluded1.insert("WTCCC94134 ");
	//("Exclude individuals with questionable ancestry");
	excluded1.insert("WTCCC94336");
	excluded1.insert("WTCCC93927");
	excluded1.insert("WTCCC94729");
	excluded1.insert("WTCCC94442");
	excluded1.insert("WTCCC93690");
	excluded1.insert("WTCCC94283");
	excluded1.insert("WTCCC94262");
	excluded1.insert("WTCCC94424");
	excluded1.insert("WTCCC94632");
	excluded1.insert("WTCCC94301");
	//("Exclude individuals who have >10% missing genotypes");
	excluded1.insert("WTCCC93997");
	excluded1.insert("WTCCC94060");
	excluded1.insert("WTCCC94556");
	excluded1.insert("WTCCC94645");
	//("Remove BRCA2 +ves");
	excluded1.insert("WTCCC94365");
	excluded1.insert("WTCCC94045");
	excluded1.insert("WTCCC93991");
	excluded1.insert("WTCCC93684");
	excluded1.insert("WTCCC93841");
	excluded1.insert("WTCCC93981");
	excluded1.insert("WTCCC93968");
	excluded1.insert("WTCCC94359");
	excluded1.insert("WTCCC94271");
	excluded1.insert("WTCCC94257");
	excluded1.insert("WTCCC94272");
	excluded1.insert("WTCCC94277");
	excluded1.insert("WTCCC94406");
	excluded1.insert("WTCCC93824");
	excluded1.insert("WTCCC93980");
	excluded1.insert("WTCCC93845");
	excluded1.insert("WTCCC94227");
	excluded1.insert("WTCCC94500");
	// AS
	excluded1.insert("WTCCC94962");
	excluded1.insert("WTCCC95203");
	excluded1.insert("WTCCC95088");
	excluded1.insert("WTCCC95105");
	excluded1.insert("WTCCC95433");
	excluded1.insert("WTCCC95926");
	excluded1.insert("WTCCC95751");
	excluded1.insert("WTCCC95229");
	excluded1.insert("WTCCC95736");
	excluded1.insert("WTCCC95356");
	excluded1.insert("WTCCC95977");
	excluded1.insert("WTCCC95754");
	excluded1.insert("WTCCC95726");
	excluded1.insert("WTCCC95771");
	// ("Exclude individuals with questionable ancestry");
	excluded1.insert("WTCCC95088");
	excluded1.insert("WTCCC95427");
	excluded1.insert("WTCCC96008");
	excluded1.insert("WTCCC96028");
	excluded1.insert("WTCCC95059");
	excluded1.insert("WTCCC95068");
	excluded1.insert("WTCCC95095");
	excluded1.insert("WTCCC95817");
	excluded1.insert("WTCCC95624");
	excluded1.insert("WTCCC95640");
	excluded1.insert("WTCCC95662");
	excluded1.insert("WTCCC95677");
	excluded1.insert("WTCCC96085");
	excluded1.insert("WTCCC96108");
	excluded1.insert("WTCCC95739");
	excluded1.insert("WTCCC95758");
	excluded1.insert("WTCCC95629");
	//("Exclude individuals who have >10% missing genotypes");
	excluded1.insert("WTCCC95625");
	excluded1.insert("WTCCC95671");
	excluded1.insert("WTCCC95749");
	excluded1.insert("WTCCC94937");
	excluded1.insert("WTCCC95126");
	// ATD
	// ...

	// geno data
	int rsIndex = -1;
	int sampIndex = -1;
	int nSkipped(0);
	long entrez;
	map<string,int> rspos;
	// add data from case dataset
	bool sampleOk = false;
	bool snpOk = false;
	unsigned max = nf1*outFile.numSnps();
	for (unsigned j=0; j<max; j++) {
		// Sample data j,0
		rsIndex = (j%outFile.numSnps());
		df.GetData(1,j,strVal);	// leggo sempre per mantenere continuita su datafile
		if(rsIndex == 0){
			//df.GetData(1,j,strVal);
			if(excluded1.find(strVal)==excluded1.end()){
				sampIndex++;
				outFile.setSampleId(sampIndex, sampIndex);
				outFile.setSampleName(sampIndex,strVal);
				outFile.setSamplePheno(sampIndex, 1);
				cout << "Cases: added sample: " << strVal << " #"<< j/outFile.numSnps() << " at line "<< j << endl;
				sampleOk = true;
			} else {
				nSkipped++;
				cout << "Cases: skipped sample: " << strVal << " #"<< j/outFile.numSnps() << " at line "<< j << endl;
				sampleOk = false;
			}
		}

		// snp name [0 - numSnp]
		if(j < outFile.numSnps()){
			df.GetData(0, rsIndex, strVal);
			if(excludedSNP.find(strVal)==excludedSNP.end()){
				rspos[strVal]=rsIndex;
				outFile.setSNPName(rsIndex, strVal);
				strChr=""; strGene=""; bp=0;
				if(dbok)
					lookupGeneBySNP(strVal, strChr, strGene, bp, entrez, conn);
				outFile.setSNPChr(rsIndex, strChr);
				outFile.setSNPPos(rsIndex, bp);
				outFile.setSNPGene(rsIndex, entrez);
				//df.GetData(5+j*2+1, i, strVal);
				//outFile.setSNPAllele(i, strVal);
				//df.GetData(4, i, strVal);
				//outFile.setSNPScore(i, atof(strVal.c_str()));
			}
		}
		//df.GetData(3, j, strVal);
		//if(atof(strVal.c_str()) > 0.0){
		if(sampleOk && snpOk){
			df.GetData(2, j, strVal);
			short gt = outFile.decodeAllele(strVal);
			outFile.setGType(sampIndex, rsIndex, gt);
		}
	}
	// add data from control dataset
	max = nf2*outFile.numSnps();
	for (unsigned j=0; j<max; j++) {
		// Sample data j,0
		df2.GetData(1,j,strVal);
		if((j%outFile.numSnps()) == 0){
			//df2.GetData(1,j,strVal);
			if(excluded1.find(strVal)==excluded1.end()){
				sampIndex++;
				outFile.setSampleId(sampIndex,sampIndex);
				outFile.setSampleName(sampIndex,strVal);
				outFile.setSamplePheno(sampIndex, 0);
				cout << "Controls: added sample: " << strVal << " #"<< j/outFile.numSnps() << " at line "<< j << endl;
				sampleOk = true;
			} else {
				nSkipped++;
				cout << "Controls: skipped sample: " << strVal << " #"<< j/outFile.numSnps() << endl;
				sampleOk = false;
			}
		}
		//df.GetData(3, j, strVal);
		//if(atof(strVal.c_str()) > 0.0){
		if(sampleOk){
			string rsVal;
			df2.GetData(0, j, rsVal);
			df2.GetData(2, j, strVal);
			outFile.setGType(sampIndex, rspos[rsVal], outFile.decodeAllele(strVal));
		}
	}
	cout << "skipped samples: " << nSkipped << endl;
	cout << "added samples: " << sampIndex+1 << endl;
	// ricodifica 0151010 -> 012 e inversione tabella snp
	//encode data 012
	T_GType gt;

	map<T_GType,T_GType> code;
	for(unsigned j=0; j<outFile.numSnps(); j++){
		//outFile.getSnpRowP("/SNPDataTableInv", j, v3);
		int done = 0;
		size_t i = 0;
		string gtStr("");
		code.clear();
		int indgt=2;
		code[1] = 1;
		do {
			outFile.getGType(i++,j,gt); //gt = v3[i++];
			if((code.find(gt)==code.end())&& (gt!=1)){	// etero sono gia  ok
				if(code.size()==1)
					code[gt]=0;
				else
					code[gt]=indgt++;
			}
		} while(i<outFile.numSamples());
		string str;
		outFile.setSNPIndex(j, code.size());
		outFile.getSNPName(j, str);
		if((code.size()>3) || (code.size()<3)){
			//cout << "snp: " << str << " has " << code.size() << " genotypes" << endl;
		} else
			cout << "snp: " << str << " has " << code.size() << " genotypes" << endl;
		// aggiorno
		for(unsigned k=0; k<outFile.numSamples(); k++){
			outFile.getGType(k,j,gt);
			//if(gt!=1)
				//if(code[gt]==0)
			outFile.setGType(k,j,code[gt]);
				//else
				//	outFile.setGType(k,j,code[gt]+1);	// 1-->2 x non confondere con etero gia marcati
			/*if(gt == code[0]) outFile.setGType(k,j,0);
			else if(gt == code[1]) outFile.setGType(k,j,1);
			else if(gt == code[2]) outFile.setGType(k,j,2);*/
		}
	}
	// inversione
	outFile.doInvDataTable();
	return 1;
}

// reads from ARFF dataset (WEKA et al.) by rows.
// and import into h5file
int arff2h5(ifstream &df, h5TPED_I &outFile)
{
	int nsamp = 0;
	int nsnp = 0;
	
	std::string strVal("");
	int cont=0;
	char buffer[255];
	bool found = false;
	while(!found){
		df.getline(buffer, 255);
		found = !strncmp(buffer,"@relation",9);
		if(found){
			df.getline(buffer, 255);
		}
	}
	// info
	for (unsigned j=0; j<outFile.numSnps(); j++) {
		df.getline(buffer, 255, ' ');
		assert(!strcmp(buffer, "@attribute"));
		df.getline(buffer, 255, ' ');
		outFile.setSNPName(j, buffer);
		df.getline(buffer, 255);
		// unknown data
		//outFile.setSNPChr(i, strVal);
		//outFile.setSNPPos(i, atoi(strVal.c_str()));
		//outFile.setSNPAllele(i, strVal);
		//outFile.setSNPScore(i, atof(strVal.c_str()));
	}

	strVal = "samp_";
	for (unsigned j=0; j<outFile.numSamples(); j++) {
		outFile.setSampleId(j,j);
		outFile.setSampleName(j,strVal+_itoa(j,buffer,10));
	}
	found = false;
	while(!found){
		df.getline(buffer, 255);
		found = !strncmp(buffer,"@data",5);
	}

	// geno data
	for(unsigned i=0; i<	outFile.numSamples(); i++){

		//Genotype i,j
		for (unsigned j=0; j<outFile.numSnps(); j++) {
			df.getline(buffer, 255, ',');
			if(!strcmp(buffer,"?")){
				outFile.setGType(i, j, -1);
			} else {
				outFile.setGType(i, j, atoi(buffer));
			}
		}
		// pheno
		df.getline(buffer, 255);
		outFile.setSamplePheno(i,atoi(buffer));

		if(cont==10000){
			cont=0;
			cout << "#";
		} else cont++;
	}
	return 1;
}




// function specific for the file coming from the merge of different turin bladder genotypes
// ISI, Decode, Xi Feng done with the help of the relative .mdb file
// remember to clean the report file from spourious dots !
// Resulting format is the one adopted by Brinza with no separation char
int encodeTurinBladder(CDataFile &df)
{
	std::string strVal("");
	std::string strVal2("");
	int cont=df.GetNumberOfSamples(0);

	typedef struct {
		char snp1;
		char snp2;
	} locus;

#define TAQMAN 1	// A/C
#define XIFENG 2	// Undetermined/Both/FAM....
#define DECODE 3	// AC

	// alloc and reset
	vector <int> coding(df.GetNumberOfVariables(), 0);
	for (int j=0; j<df.GetNumberOfVariables()-1; j++) {
		int k=0;
		do{
			df.GetData(j, k++, strVal);
		} while((strVal == df.NA()) && (k < cont));

		if((strVal.find("A/A")!= string::npos)||(strVal.find("A/C")!= string::npos)||(strVal.find("A/G")!= string::npos)||(strVal.find("A/T")!= string::npos)
			||(strVal.find("C/C")!= string::npos)||(strVal.find("C/G")!= string::npos)||(strVal.find("C/T")!= string::npos)
			||(strVal.find("G/G")!= string::npos)||(strVal.find("G/T")!= string::npos)||(strVal.find("T/T")!= string::npos)
			||(strVal.find("D/D")!= string::npos)){
			coding.at(j) = TAQMAN;
		} else if((strVal.find("Undetermined")!= string::npos)||(strVal.find("Probe")!= string::npos)||(strVal.find("Both")!= string::npos)||(strVal.find("both")!= string::npos)||(strVal.find("Ch")!= string::npos)){
			coding.at(j) = XIFENG;
		} else if((strVal.find("AA")!= string::npos)||(strVal.find("AC")!= string::npos)||(strVal.find("AG")!= string::npos)||(strVal.find("AT")!= string::npos)
			||(strVal.find("CC")!= string::npos)||(strVal.find("CG")!= string::npos)||(strVal.find("CT")!= string::npos)
			||(strVal.find("GG")!= string::npos)||(strVal.find("GT")!= string::npos)||(strVal.find("TT")!= string::npos)){
			coding.at(j) = DECODE;
		}
	}

	vector <char> locus1(139);
	vector <char> locus2(139);
	// Sample data from csv file
	for (int j=0; j<df.GetNumberOfSamples(0); j++) {
		for (int i=0; i<df.GetNumberOfVariables()-1; i++) {	
			if(coding.at(i) == TAQMAN){
				strVal = "";
				bool found = false;
				int k=0;
				while((!found) && (0==j)){
					df.GetData(i, k++, strVal);
					if(!strVal.empty()){
						if(strVal.at(0) != strVal.at(2)){
							locus1.at(i) = strVal.at(0);
							locus2.at(i) = strVal.at(2);
							found=true;
						}
					}
					if((!found)&&(k==cont-1)){
						locus1.at(i) = 'A';
						locus2.at(i) = 'A';
						found = true;
					}
				}
				strVal2="3";
				df.GetData(i, j, strVal);
				if(!strVal.empty()){
					if((strVal.at(0) == strVal.at(2)) &&(strVal.at(0)==locus1.at(i))){
						strVal2 = "0";
					} else if((strVal.at(0) == strVal.at(2)) &&(strVal.at(0)==locus2.at(i))){
						strVal2 = "2";
					} else if(strVal.at(0) != strVal.at(2)){
						strVal2 = "1";
					}
				}
				df.SetData(i, j, strVal2.c_str());
			} 
			else if(coding.at(i) == DECODE) { 
				std::string strVal("");
				bool found = false;
				int k=0;
				while((!found) && (0==j)){
					df.GetData(i, k++, strVal);
					if(!strVal.empty()){
						if(strVal.at(0) != strVal.at(1)){
							locus1.at(i) = strVal.at(0);
							locus2.at(i) = strVal.at(1);
							found=true;
						}
					}
					if((!found)&&(k==cont-1)){
						locus1.at(i) = 'A';
						locus2.at(i) = 'A';
						found = true;
					}
				}
				strVal2="3";
				df.GetData(i, j, strVal);
				if(!strVal.empty()){
					if((strVal.at(0) == strVal.at(1)) &&(strVal.at(0)==locus1.at(i))){
						strVal2 = "0";
					} else if((strVal.at(0) == strVal.at(1)) &&(strVal.at(0)==locus2.at(i))){
						strVal2 = "2";
					} else if(strVal.at(0) != strVal.at(1)){
						strVal2 = "1";
					}
				}
				df.SetData(i, j, strVal2.c_str());
			}
			else if(coding.at(i) == XIFENG){
				std::string strVal("");
				if(0==j){
					locus1.at(i) = 'A';
					locus2.at(i) = 'B';
				}
				strVal2="3";
				df.GetData(i, j, strVal);
				if((!strVal.empty()) && (strVal != "Undetermined")){
					if(strVal.find("- FAM")!= string::npos ){
						strVal2 = "0";
					} else if(strVal.find("- VIC")!= string::npos ){
						strVal2 = "2";
					} else if(strVal.find("Both")!= string::npos ){
						strVal2 = "1";
					}
				}
				df.SetData(i, j, strVal2.c_str());
			}
		}
	}
	std::string fname=df.m_szFilename+".norm.csv";
	df.WriteFile(fname.c_str());
	return 1;
}

