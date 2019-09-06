/*                                                                                */
/*   Copyright  (C) 2005 Intel Corporation.  All rights reserved.                 */
/*                                                                                */
/*   The information and source code contained herein is the exclusive property   */
/*   of Intel Corporation and may not be disclosed, examined, or reproduced in    */
/*   whole or in part without explicit written authorization from the Company.    */
/*                                                                                */
/*                                                                                */

//#include "time/mytime.h"
#include "hdfvar_I.h"


//-----------------------------------------------------------------------------------------------------------
//                           hdfvar implemented with direct use of hdf5 library
//-----------------------------------------------------------------------------------------------------------
// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables
template<typename _GType, typename _ExpType>
h5TPED_I<_GType, _ExpType>::h5TPED_I() {
		m_filename = "";
		m_nSnp = 0; // table 1xN
		m_Chr.clear();
		m_nSamples = 0;// table 1xM
		m_SamplesInRows = true;
		m_NA = "NA";
	}

//
template<typename _GType, typename _ExpType>
h5TPED_I<_GType, _ExpType>::h5TPED_I(const std::string &szFilename) : H5::H5File(szFilename, H5F_ACC_RDWR) {
		m_filename = szFilename;
		m_nSnp = 0; // table 1xN
		m_Chr.clear();
		m_nSamples = 0;// table 1xM
		m_SamplesInRows = true;
		m_NA = "NA";
	}

// Misc. constructor.  Instantiates an instance of h5Gtype file and build folder and tables
template<typename _GType, typename _ExpType>
h5TPED_I<_GType, _ExpType>::h5TPED_I(const std::string &szFilename, const int nsnp,
			const int nsamp, const bool szSampInRows)  : H5::H5File(szFilename, H5F_ACC_RDWR) {

		m_filename = szFilename;
		m_nSnp = nsnp;
		m_Chr.clear();
		m_nSamples = nsamp;
		m_SamplesInRows = szSampInRows;
		m_NA = "NA";
		//buildStruct();
	}
	
template class h5TPED_I<char, float>;