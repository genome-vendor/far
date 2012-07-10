/*
 * AlignmentFilter.h
 *
 *  Created on: Jun 29, 2010
 *      Author: mat
 */

#ifndef ALIGNMENTFILTER_H_
#define ALIGNMENTFILTER_H_

#include <sstream>
#include "SequencingRead.h"
#include <seqan/graph_align.h>
#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>
#include "Enums.h"

#include "SequenceConverter.h"

/**This class does the actual alignment via the passed Algorithm. It has an internal vector of adapter sequences and will align each adapter to the read. The one with the highest score will be used for the final alignment. */
template <class TString, class TIDString, class TAlgorithm >
class AlignmentFilter {
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
	//TString m_adapter;
	int m_cutOff;
	tbb::concurrent_vector<int> overlaps;
	int m_scoreMatch;
	int m_scoreMismatch;
	int m_scoreGapInsert;
	int m_scoreGapExtend;
	tbb::atomic<bool> m_writeTag;
	//TAlgorithm *m_alg;

	typedef ::std::pair<SequencingRead<seqan::CharString,seqan::CharString>*,unsigned int> TAdapter;
	tbb::concurrent_vector<TAdapter> *m_adapters;



public:
	AlignmentFilter(tbb::concurrent_vector<TAdapter> *adapters, int scoreMatch, int scoreMismatch, int scoreGapInsert, int scoreGapExtend, FlexAR::TrimEnd end,FlexAR::LogLevel log = FlexAR::NONE) {//: tbb::filter(parallel) {
		m_verbose = log;
		m_trimEnd = end;
		m_adapters = adapters;
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
		m_writeTag = false;
		//if(m_verbose==FlexAR::TAB)std::cout << "read-ID\tadapter-ID\tadapter-start\tadapter-end\toverlap-length\tmismatches\tindels\tallowed-errors" << std::endl;
	};

	void setModifiedTag(TString tag){
		m_tag = tag;
		m_writeTag = true;
	}

	void setCutOff(float cutoff){
		m_cutOff = static_cast<int>(cutoff);
	}

	float getCutOff(){
		return static_cast<float>(m_cutOff);
	}

	void setMinOverlap(int min_overlap, bool isAdaptiveOverlap){
		//std::cout << "Overlap-length vector: " << std::endl;
		m_minOverlap = min_overlap;
		if((isAdaptiveOverlap)&&((m_trimEnd==FlexAR::RIGHT_TAIL)||(m_trimEnd==FlexAR::RIGHT))){
			for(unsigned int n=0;n<FlexAR::MAX_READLENGTH;++n){
				if(n < ::FlexAR::MAX_READLENGTH-min_overlap+1){
					m_overlaps[n] = min_overlap;
				}
				else {
					--min_overlap;
					m_overlaps[n] = min_overlap;
			//		std::cout << n << ":" << m_overlaps[n] << ",";
				}

			}
			//std::cout << std::endl;
		}

		if((isAdaptiveOverlap)&&((m_trimEnd==FlexAR::LEFT_TAIL)||(m_trimEnd==FlexAR::LEFT))){
			for(int n=FlexAR::MAX_READLENGTH - 1;n>0;--n){
				if(n < min_overlap+1){
					m_overlaps[n] = min_overlap;
					--min_overlap;
				//	std::cout << n << ":" << m_overlaps[n] << ",";
				}
				else {

					m_overlaps[n] = min_overlap;
				}

			}
			//std::cout << std::endl;
		}

	if(!isAdaptiveOverlap){
			for(unsigned int n=0;n<FlexAR::MAX_READLENGTH;++n){
				m_overlaps[n] = min_overlap;
				//std::cout << m_overlaps[n] << ",";
			}
			//std::cout << std::endl;
		}
	}

	int getMinOverlap(){
		return m_minOverlap;
	}

	int getMinOverlapLength(){
		return m_minOverlapLength;
	}

	int getMaxOverlapLength(){
		return m_maxOverlapLength;
	}

	int getMedianOverlapLength(){
		if(overlaps.size()>0){
			std::sort(overlaps.begin(), overlaps.end());

			return overlaps.at(overlaps.size()/2);

		}
		else return 0;
	}

	int getMeanOverlapLength(){
		if(getNrModifiedReads() > 0)return m_sumLength / getNrModifiedReads();
		else return -1;

	}

	int getNrModifiedReads()
	{
		return m_modified;
	}

	virtual ~AlignmentFilter(){

	};

	//void* operator()( void* item ){
	/*! This function aligns all adapter sequences passed to the class via constructor to the read that is passed as an argument.
	The read will be cut if wanted*/
	int align( void* item , bool postCutRead = false){
		std::stringstream ss;
		std::string finalAlString;
		SequencingRead<TString,TIDString> &myRead = *static_cast< SequencingRead<TString,TIDString>* >(item);

		using namespace seqan;
		using namespace std;

		int fmismatches=0,readLength=0,freadLength=0,fstartPos=0,fstartPosA=0,fstartPosS=0,fendPos=0,fendPosS=0,fendPosA=0,foverlapLength=0,adapterLength=0,overlapLength=0,fgaps=0,fIndex=-1;
		float fallowedErrors = 0,scoreMax=-10000,scoreMin=10000;
		TString read="",bpRead="",quality="";

		///Preprocessing of the read: only align suffixes if colorspace format...
		switch(myRead.getFormat()){
			case FlexAR::CSFASTQ:
				read = suffix<TString>(myRead.getSequence(),2);
				quality = suffix<TString>(myRead.getQuality(),1);
				break;
			case FlexAR::FASTQ:
				read = myRead.getSequence();
				quality = myRead.getQuality();
				break;
			case FlexAR::CSFASTA:
				read = suffix<TString>(myRead.getSequence(),2);
				quality = "";
				break;
			case FlexAR::FASTA:
				read = myRead.getSequence();
				quality = "";
				break;
		}

		//common calculations
		readLength = length(read);

		TString sequence="",squality = "";
		TString current_adapter;
		TAlgorithm alg(m_scoreMatch,m_scoreMismatch,m_scoreGapInsert,m_scoreGapExtend);

		///Iterate over all passed adapter sequences, remember the best one
		for(unsigned int i=0;i<m_adapters->size();++i)
		{
			current_adapter = m_adapters->at(i).first->getSequence();
			adapterLength = length(current_adapter);
			//In tail mode trim read before alignment
			if((m_trimEnd==FlexAR::LEFT_TAIL)||(m_trimEnd==FlexAR::RIGHT_TAIL))
			{
				//only search adapter length
				if(adapterLength<readLength){
					//only search the last part of sequence (m_adaptor length) - speedup
					if(m_trimEnd == FlexAR::LEFT_TAIL){
						sequence = seqan::prefix<TString>(read,adapterLength);
						if(myRead.getFormat()==FlexAR::FASTQ)squality = seqan::prefix<TString>(quality,adapterLength);
						if(myRead.getFormat()==FlexAR::CSFASTQ)squality = seqan::prefix<TString>(quality,adapterLength-1);
					}
					else {

						sequence = seqan::suffix<TString>(read,readLength - adapterLength);
						if(myRead.getFormat()==FlexAR::FASTQ)squality = seqan::suffix<TString>(quality,readLength - adapterLength);
						if(myRead.getFormat()==FlexAR::CSFASTQ)squality = seqan::suffix<TString>(quality,readLength - adapterLength-1);
					}
				} else {
					// if m_adaptor is longer than sequence use prefix/suffix of m_adaptor

					current_adapter = seqan::prefix<TString>(m_adapters->at(i).first->getSequence(),readLength);

					sequence = read;
					if(myRead.getFormat()==FlexAR::FASTQ)squality = seqan::prefix<TString>(quality,adapterLength);
					if(myRead.getFormat()==FlexAR::CSFASTQ)squality = seqan::prefix<TString>(quality,adapterLength-1);

					if((m_verbose==FlexAR::ALL)||(m_verbose==FlexAR::CHANGED)){
						ss << "Adapter was trimmed to readlength:" << current_adapter << std::endl;
						ss << endl << endl << "Sequence ID:" << myRead.getSequenceTag() << std::endl;
					}
				}
			} else {
				///If trim mode was NOT set to --trim-end tail_left or tail_right:
				sequence = read;
				squality = quality;
			}

			int startPos=0,endPos=0,startPosA=0,endPosA=0,startPosS=0,endPosS=0,alignmentScore=0,mismatches=0,gaps=0;
			std::stringstream alString;

			///Align the current_adapter via the passed algorithm
			alg.align(current_adapter,sequence, squality, gaps, mismatches, startPos, endPos, startPosA, endPosA, startPosS, endPosS, alignmentScore, alString);

			//remember min/max alignment score
			//if(static_cast<float>(alignmentScore)>scoreMax)scoreMax = static_cast<float>(alignmentScore);
			if(static_cast<float>(alignmentScore)<scoreMin)scoreMin = static_cast<float>(alignmentScore);

			//calculate overlap
			overlapLength = endPos - startPos;

			//calculate allowed errors in current sequence based on cutoff (errors per 10 bases)
			float allowedErrors = m_cutOff * static_cast<float>(overlapLength/10.0f);

			//check if allowed errors and overlapLength are below thresh
			if((static_cast<float>(mismatches+gaps) <= allowedErrors)&&(this->m_overlaps[FlexAR::MAX_READLENGTH - readLength + startPos] <= overlapLength)&&(alignmentScore > scoreMax))
			{
				//ss << "Current overlap threshold: " << this->m_overlaps[FlexAR::MAX_READLENGTH - readLength + startPos] << endl;
				//Remeber the highest alignmentScore
				scoreMax = static_cast<float>(alignmentScore);
				fstartPos = startPos;
				fstartPosA = startPosA;
				fstartPosS = startPosS;
				fendPos = endPos;
				fendPosA = endPosA;
				fendPosS = endPosS;
				freadLength = readLength;
				fallowedErrors = allowedErrors;
				fgaps = gaps;
				fmismatches = mismatches;
				finalAlString = alString.str();
				fIndex = i;
			}
		}

		++fIndex;

		///////////////////////////////////////////////////////
		///	Cutting the read								///
		///////////////////////////////////////////////////////

		bool mod = false;

		if((postCutRead)&&(fIndex>0)){

			//Now cut the read depending on the best adapter match
			FlexAR::TrimEnd trimEnd = m_trimEnd;

			if(trimEnd==FlexAR::ANY){
				if(fstartPos > (freadLength - fendPos)){

					trimEnd = FlexAR::RIGHT;
					ss << "Trimming from right..." << endl;
				}
				else {

					trimEnd = FlexAR::LEFT;
					ss << "Trimming from left..." << endl;
				}
			}

			switch(trimEnd)
			{
					//break;
				case FlexAR::LEFT_TAIL:
					sequence = read;
					squality = quality;
					fendPos = overlapLength;
					//Now do the same as trimming from LEFT
				case FlexAR::LEFT:
					if(fendPosA<=fendPosS){
						//translate coordinates in alignment structure
						if(fstartPosS > 0)fendPos -= fstartPosS;
						mod = true;
						//Assume adapter covers read in the interval adapterpos < endRead
						if((myRead.getFormat()==FlexAR::FASTA)||(myRead.getFormat()==FlexAR::FASTQ))
						{
							if(fendPos > freadLength)fendPos = freadLength;
							erase(sequence,0,fendPos);
							myRead.setSequence(sequence);

							if(myRead.getFormat()==FlexAR::FASTQ){
								erase(squality,0,fendPos);
								myRead.setQuality(squality);
							}
						}
						else
						{
							//colorspace
							//First 2 Bases are ignored in colorspace. if trimmed from left  they have to be

							TString result;
							//Get full read as basepair sequence
							bpRead = SequenceConverter<TString>::getInstance()->colorSpaceToBasepairSpace(prefix<TString>(myRead.getSequence(),fendPos + 3));
							//delete one more base in colorspace
							if(fendPos + 1 < freadLength)
							{
								erase(sequence,0,fendPos+1);
								//Build prefix of read (TX)
								result = SequenceConverter<TString>::getInstance()->getColorcodeFromT(bpRead[fendPos + 2]);
								TString colRead = "T";
								append(colRead,result);
								result = "";
								append(result,colRead);

								append(result, sequence);
								myRead.setSequence(result);
								// Trim quality string
								if(myRead.getFormat()==FlexAR::CSFASTQ){
									erase(squality,0,fendPos);
									myRead.setQuality(squality);

								}
							}
							else {
								myRead.setQuality("");
								myRead.setSequence("");
							}

						}
					}
					break;
				case FlexAR::RIGHT_TAIL:
					sequence = read;
					squality = quality;
					freadLength = length(read);
					fstartPos = freadLength - overlapLength;

				case FlexAR::RIGHT:
					if((fstartPosA>=0)&&(fstartPosS==0)){
						mod = true;
						if((myRead.getFormat()==FlexAR::FASTA)||(myRead.getFormat()==FlexAR::FASTQ))
						{
							erase(sequence,fstartPos,freadLength);
							myRead.setSequence(sequence);
							if(myRead.getFormat()==FlexAR::FASTQ){
								erase(squality,fstartPos,freadLength);
								myRead.setQuality(squality);
							}
						}
						else
						{
							if(fstartPos-1 >= 0){
								//colorspace
								erase(sequence,fstartPos-1,freadLength);
								//Append TX
								TString result;
								result = prefix(myRead.getSequence(),2);
								append(result, sequence);
								myRead.setSequence(result);
								if(myRead.getFormat()==FlexAR::CSFASTQ){
									erase(squality,fstartPos,freadLength);
									myRead.setQuality(squality);
								}
							} else {
								myRead.setQuality("");
								myRead.setSequence("");
							}
						}
					}
					break;
				case FlexAR::ANY:
					break;
			}

			if(mod){
				if(m_writeTag){
					TString newTag = myRead.getSequenceTag();
					append(newTag,":");
					append(newTag,m_tag);
					myRead.setSequenceTag(newTag);
				}
				++m_modified;
				//count for each adapter how often it was removed... d
				m_adapters->at(fIndex-1).second += 1;

				//output alignment
				if((m_verbose==FlexAR::CHANGED)||(m_verbose==FlexAR::ALL)){
					ss << endl;
					ss << "read-tag     : " << myRead.getSequenceTag() << endl;
					ss << "adapter-tag	: " << m_adapters->at(fIndex-1).first->getSequenceTag()<< endl;
					ss << "read        : " << read << endl;
					ss << "quality     : " << quality << endl;

					ss << "beginR       : " << fstartPosS << endl;
					ss << "endR         : " << fendPosS << endl;

					ss << "beginA       : " << fstartPosA << endl;
					ss << "endA         : " << fendPosA << endl;

					ss << finalAlString << endl;
					ss << "overlap(bp) : " << overlapLength << endl;
					ss << "score       : " << scoreMax << endl;

					ss << "indels        : " << fgaps << endl;
					ss << "mismatches    : " << fmismatches << endl;
					ss << "allowed errors: " << fallowedErrors << endl;
					ss << "cutOff        : " << m_cutOff << endl;

					ss << endl;
					ss << "Adapter removed !" << endl << endl;
					//cout << ss.str() << std::endl;
				}

				if(m_verbose==FlexAR::TAB)ss << myRead.getSequenceTag() << "\t" << m_adapters->at(fIndex-1).first->getSequenceTag() << "\t" << fstartPosA << "\t"<< fendPosA << "\t" << overlapLength << "\t" << fmismatches << "\t" << fgaps << "\t" << fallowedErrors << std::endl;

				//calculate summed overlap length
				m_sumLength += overlapLength;

				//store for median calculation
				overlaps.push_back(overlapLength);

				//calculate max/min overlap length
				if(overlapLength>m_maxOverlapLength)m_maxOverlapLength=overlapLength;
				if(overlapLength<m_minOverlapLength)m_minOverlapLength=overlapLength;

				if((m_verbose==FlexAR::ALL)||(m_verbose==FlexAR::CHANGED))
				{
					ss << "Trimmed sequence : " << myRead.getSequence() << std::endl;
					if((myRead.getFormat()==FlexAR::FASTQ)||(myRead.getFormat()==FlexAR::CSFASTQ))
					{
						ss << "Trimmed quality  : " << myRead.getQuality() << std::endl<< std::endl;
					}
				}

				//break;//for loop
			}//if(mod) end

		}//if(postCutRead...)
		
		//the read was not modified - output the alignment		
		if(!mod&&(m_verbose==FlexAR::ALL)){
			ss << endl;
			ss << "read-tag     : " << myRead.getSequenceTag() << endl;
			ss << "read        : " << read << endl;
			ss << "quality     : " << quality << endl;
			ss << endl;

			ss << "No adapter aligned with the given thresholds for this read. " << std::endl;
		}

		//Bundeled output for multithreading
		if(m_verbose!=FlexAR::NONE){
			if(fIndex!=-1)cout << ss.str();
		}

		return fIndex;
		//return &myRead;

	};



};

#endif /* ALIGNMENTFILTER_H_ */
