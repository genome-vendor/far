/*
 * MultiplexedAlignmentFilter.h
 *
 *  Created on: Jun 1, 2011
 *      Author: mat
 */

#ifndef MULTIPLEXEDALIGNMENTFILTER_H_
#define MULTIPLEXEDALIGNMENTFILTER_H_

#include <sstream>
#include "MultiplexedRead.h"
#include "NeedlemanWunschAlignmentAlgorithm.h"
#include "ExactStringMatchingAlgorithm.h"
#include "AlignmentFilter.h"
#include <seqan/graph_align.h>
#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>
#include "Enums.h"
#include "SequenceConverter.h"


/**This class processes a MultiplexedRead. It will assign a barcode-ID to the read (if it is processing a barcoded run)
and remove adapter sequences. */
template <class TString, class TIDString, class TAlgorithm >
class MultiplexedAlignmentFilter : public tbb::filter{
private:
	tbb::atomic<unsigned long> m_sumLength;
	int m_minOverlap;
	tbb::atomic<FlexAR::TrimEnd> m_trimEnd;
	FlexAR::LogLevel m_verbose;
	tbb::atomic<int> m_minOverlapLength;
	tbb::atomic<int> m_maxOverlapLength;
	tbb::atomic<int> m_modified;
	int m_overlaps[::FlexAR::MAX_READLENGTH];
	TString m_tag;
	int m_cutOff;
	tbb::concurrent_vector<int> overlaps;
	int m_scoreMatch;
	int m_scoreMismatch;
	int m_scoreGapInsert;
	int m_scoreGapExtend;
	int m_bThresh;
	FlexAR::TrimEnd m_bEnd;
	tbb::atomic<bool> m_writeTag;
	bool m_isAdaptiveOverlap;
	bool m_demultiplexOnly;

	typedef ::std::pair<SequencingRead<seqan::CharString,seqan::CharString>*,unsigned int> TAdapter;
	tbb::concurrent_vector<TAdapter> *m_adapters;
	tbb::concurrent_vector<TAdapter> *m_barcodes;
	AlignmentFilter<TString,TIDString,NeedlemanWunschAlignmentAlgorithm<TString> > *m_filter;
	AlignmentFilter<TString,TIDString,NeedlemanWunschAlignmentAlgorithm<TString> > *m_bfilter;
public:
	MultiplexedAlignmentFilter(tbb::concurrent_vector<TAdapter> *adapters, tbb::concurrent_vector<TAdapter> *barcodes,int bThres,int barcode_min_overlap, FlexAR::TrimEnd bend, int scoreMatch, int scoreMismatch, int scoreGapInsert, int scoreGapExtend, FlexAR::TrimEnd end,FlexAR::LogLevel log = FlexAR::NONE) : tbb::filter(parallel) {
		m_verbose = log;
		m_bEnd = bend;
		m_trimEnd = end;
		m_adapters = adapters;
		m_barcodes = barcodes;
		m_minOverlap = 10;
		m_cutOff = 1;
		m_modified = 0;
		m_sumLength = 0;
		m_maxOverlapLength=0;
		m_minOverlapLength=1000;
		m_scoreMatch = scoreMatch;
		m_scoreMismatch = scoreMismatch;
		m_scoreGapInsert = scoreGapInsert;
		m_scoreGapExtend = scoreGapExtend;
		m_tag = "";
		m_minOverlap = 10;
		m_writeTag = false;
		m_bThresh = bThres;
		m_isAdaptiveOverlap = false;
		m_filter = new AlignmentFilter<TString,TIDString,NeedlemanWunschAlignmentAlgorithm<TString> >(m_adapters,m_scoreMatch,m_scoreMismatch,m_scoreGapInsert,m_scoreGapExtend, m_trimEnd, m_verbose);
		m_filter->setCutOff(m_cutOff);
		m_filter->setMinOverlap(m_minOverlap,m_isAdaptiveOverlap);
		m_bfilter = new AlignmentFilter<TString,TIDString, NeedlemanWunschAlignmentAlgorithm<TString> >(m_barcodes,3,-3,-5,-1, m_bEnd, m_verbose);
		m_bfilter->setMinOverlap(barcode_min_overlap, false);
		m_bEnd = bend;
		m_demultiplexOnly = false;
	}
	virtual ~MultiplexedAlignmentFilter(){};

	void setCutOff(float cutoff){
		m_cutOff = static_cast<int>(cutoff);
		m_filter->setCutOff(m_cutOff);
	}

	void setDemultiplexOnly(){
		m_demultiplexOnly = true;
	}

	void setMinOverlap(int min_overlap, bool isAdaptiveOverlap){
		m_minOverlap = min_overlap;
		m_isAdaptiveOverlap = isAdaptiveOverlap;
		m_filter->setMinOverlap(m_minOverlap,m_isAdaptiveOverlap);
	}

	void* operator()( void* item ){
		if(item!=NULL){
			MultiplexedRead<TString,TIDString> *myRead = static_cast< MultiplexedRead<TString,TIDString>* >(item);

			if(m_barcodes->size()!=0){
				//do barcode classification
				if(m_bEnd==FlexAR::OFF)myRead->m_barcode_id = m_bfilter->align(myRead->m_b,false);
				else myRead->m_barcode_id = m_bfilter->align(myRead->m_r1,true);
				if(m_verbose==FlexAR::ALL)std::cout << std::endl << "Barcode classification - ID: " << myRead->m_barcode_id << std::endl;
			}

			if(!m_demultiplexOnly){
				if(m_verbose==FlexAR::ALL)std::cout << std::endl << "Alignment of adapters: " << std::endl;
				m_filter->align(myRead->m_r1,true);
				//if it's a paired end run...
				if(myRead->m_r2!=NULL)m_filter->align(myRead->m_r2,true);
			}
			return myRead;
		}
		else return NULL;
	}

	int getNrModifiedReads()
	{
		return m_filter->getNrModifiedReads();
	}

	void printAlignmentSummary(){
		if(m_filter->getNrModifiedReads()>0)std::cout << "Min-/Max-/Mean-/Median-overlap length: "<< m_filter->getMinOverlapLength()  << " / "<< m_filter->getMaxOverlapLength() << " / " << m_filter->getMeanOverlapLength() << " / " << m_filter->getMedianOverlapLength() << std::endl;
	}
};

#endif /* MULTIPLEXEDALIGNMENTFILTER_H_ */
