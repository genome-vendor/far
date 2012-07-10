/*
 * AlignmentFilter.h
 *
 *  Created on: Jun 29, 2010
 *      Author: mat
 */

#ifndef AdapterLoader_H
#define AdapterLoader_H

#include <sstream>
#include "SequencingRead.h"
#include <seqan/graph_align.h>
#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>
#include "Enums.h"
#include "SequenceConverter.h"

#include <seqan/basic.h>

template <class TString, class TIDString>
class AdapterLoader : public tbb::filter{
private:
	typedef ::std::pair<SequencingRead<seqan::CharString,seqan::CharString>*,unsigned int> TAdapter;
	tbb::concurrent_vector<TAdapter> adapters;


public:
	AdapterLoader() : tbb::filter(serial) {

	};


	virtual ~AdapterLoader(){

	};

	void* operator()( void* item ){
		std::stringstream ss;
		SequencingRead<TString,TIDString> *myRead = static_cast< SequencingRead<TString,TIDString>* >(item);

		TAdapter adap;
		adap.first = myRead;
		adap.second = 0;
		adapters.push_back(adap);

		return NULL;
	};


	tbb::concurrent_vector<TAdapter> getAdapters(){
		return adapters;
	}



};

#endif /* ALIGNMENTFILTER_H_ */
