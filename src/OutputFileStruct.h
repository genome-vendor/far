/*
 * OutputFileStruct.h
 *
 *  Created on: Jun 7, 2011
 *      Author: mat
 */

#ifndef OUTPUTFILESTRUCT_H_
#define OUTPUTFILESTRUCT_H_

#include "SequenceOutputFilter.h"
template <typename TString, typename TIDString>
/*
 * Structure to store statistics for each generated FASTQ file (how many reads were discarded due to beeing to short in advance, etc.)
 */
class OutputFileStruct {
public:

	typedef SequenceOutputFilter<TString,TIDString> TOutputFilter;
	unsigned int m_cnt_short_1;
	unsigned int m_cnt_short_2;
	TOutputFilter *f1;
	TOutputFilter *f2;
	TOutputFilter *single1;
	TOutputFilter *single2;

	OutputFileStruct(){
		m_cnt_short_1 = 0;
		m_cnt_short_2 = 0;
	};
	virtual ~OutputFileStruct(){
	};
};

#endif /* OUTPUTFILESTRUCT_H_ */
