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
template <class TString, class TIDString, class TAlgorithm >
class AlignmentFilter : public tbb::filter{
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
	TAlgorithm *m_alg;

	typedef ::std::pair<SequencingRead<seqan::CharString,seqan::CharString>*,unsigned int> TAdapter;
	tbb::concurrent_vector<TAdapter> *m_adapters;



public:
	AlignmentFilter(tbb::concurrent_vector<TAdapter> *adapters, int scoreMatch, int scoreMismatch, int scoreGapInsert, int scoreGapExtend, FlexAR::TrimEnd end,FlexAR::LogLevel log = FlexAR::NONE) : tbb::filter(parallel) {
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
		if(m_verbose==FlexAR::TAB)std::cout << "read-ID\tadapter-start\tadapter-end\toverlap-length\tmismatches\tindels\tallowed-errors" << std::endl;
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
		std::cout << "Overlap-length vector: " << std::endl;
		m_minOverlap = min_overlap;
		if((isAdaptiveOverlap)&&((m_trimEnd==FlexAR::RIGHT_TAIL)||(m_trimEnd==FlexAR::RIGHT))){
			for(unsigned int n=0;n<FlexAR::MAX_READLENGTH;++n){
				if(n < ::FlexAR::MAX_READLENGTH-min_overlap+1){
					m_overlaps[n] = min_overlap;
				}
				else {
					--min_overlap;
					m_overlaps[n] = min_overlap;
					std::cout << n << ":" << m_overlaps[n] << ",";
				}

			}
			std::cout << std::endl;
		}

		if((isAdaptiveOverlap)&&((m_trimEnd==FlexAR::LEFT_TAIL)||(m_trimEnd==FlexAR::LEFT))){
			for(unsigned int n=FlexAR::MAX_READLENGTH - 1;n>0;--n){
				if(n < min_overlap+1){
					m_overlaps[n] = min_overlap;
					--min_overlap;
					std::cout << n << ":" << m_overlaps[n] << ",";
				}
				else {

					m_overlaps[n] = min_overlap;
				}

			}
			std::cout << std::endl;
		}

	if(!isAdaptiveOverlap){
			for(unsigned int n=0;n<FlexAR::MAX_READLENGTH;++n){
				m_overlaps[n] = min_overlap;
				std::cout << m_overlaps[n] << ",";
			}
			std::cout << std::endl;
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

	void* operator()( void* item ){
		std::stringstream ss;
		SequencingRead<TString,TIDString> &myRead = *static_cast< SequencingRead<TString,TIDString>* >(item);

		using namespace seqan;
		using namespace std;

		int readLength=0,adapterLength=0,overlapLength=0;

		TString read="",bpRead="",quality="";

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
				sequence = read;
				squality = quality;
			}

			int startPos=0,endPos=0,startPosA=0,endPosA=0,startPosS=0,endPosS=0,alignmentScore=0,mismatches=0,gaps=0;
			float scoreMax=-10000,scoreMin=10000;

			std::stringstream alString;


			alg.align(current_adapter,sequence, squality, gaps, mismatches, startPos, endPos, startPosA, endPosA, startPosS, endPosS, alignmentScore, alString);

			//remember min/max alignment score
			if(static_cast<float>(alignmentScore)>scoreMax)scoreMax = static_cast<float>(alignmentScore);
			if(static_cast<float>(alignmentScore)<scoreMin)scoreMin = static_cast<float>(alignmentScore);

			//calculate overlap
			overlapLength = endPos - startPos;

			//calculate allowed errors in current sequence based on cutoff (errors per 10 bases)
			float allowedErrors = m_cutOff * static_cast<float>(overlapLength/10.0f);

			ss << "Current overlap threshold: " << this->m_overlaps[FlexAR::MAX_READLENGTH - readLength + startPos] << endl;
			//check if allowed errors and overlapLength are below thresh
			if((static_cast<float>(mismatches+gaps) <= allowedErrors)&&(this->m_overlaps[FlexAR::MAX_READLENGTH - readLength + startPos] <= overlapLength))
			{

				bool mod = false;

				FlexAR::TrimEnd trimEnd = m_trimEnd;

				if(trimEnd==FlexAR::ANY){
					if(startPos > (readLength - endPos)){

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
						endPos = overlapLength;
						//Now do the same as trimming from LEFT
					case FlexAR::LEFT:
						if(endPosA<=endPosS){
							//translate coordinates in alignment structure
							if(startPosS > 0)endPos -= startPosS;
							mod = true;
							//Assume adapter covers read in the interval adapterpos < endRead
							if((myRead.getFormat()==FlexAR::FASTA)||(myRead.getFormat()==FlexAR::FASTQ))
							{
								if(endPos > readLength)endPos = readLength;
								erase(sequence,0,endPos);
								myRead.setSequence(sequence);

								if(myRead.getFormat()==FlexAR::FASTQ){
									erase(squality,0,endPos);
									myRead.setQuality(squality);
								}
							}
							else
							{
								//colorspace
								//First 2 Bases are ignored in colorspace. if trimmed from left  they have to be

								TString result;
								//Get full read as basepair sequence
								bpRead = SequenceConverter<TString>::getInstance()->colorSpaceToBasepairSpace(prefix<TString>(myRead.getSequence(),endPos + 3));
								//delete one more base in colorspace
								if(endPos + 1 < readLength)
								{
									erase(sequence,0,endPos+1);
									//Build prefix of read (TX)
									result = SequenceConverter<TString>::getInstance()->getColorcodeFromT(bpRead[endPos + 2]);
									TString colRead = "T";
									append(colRead,result);
									result = "";
									append(result,colRead);

									append(result, sequence);
									myRead.setSequence(result);
									// Trim quality string
									if(myRead.getFormat()==FlexAR::CSFASTQ){
										erase(squality,0,endPos);
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
						readLength = length(read);
						startPos = readLength - overlapLength;

					case FlexAR::RIGHT:
						if((startPosA>=0)&&(startPosS==0)){
							mod = true;
							if((myRead.getFormat()==FlexAR::FASTA)||(myRead.getFormat()==FlexAR::FASTQ))
							{
								erase(sequence,startPos,readLength);
								myRead.setSequence(sequence);
								if(myRead.getFormat()==FlexAR::FASTQ){
									erase(squality,startPos,readLength);
									myRead.setQuality(squality);
								}
							}
							else
							{
								if(startPos-1 >= 0){
									//colorspace
									erase(sequence,startPos-1,readLength);
									//Append TX
									TString result;
									result = prefix(myRead.getSequence(),2);
									append(result, sequence);
									myRead.setSequence(result);
									if(myRead.getFormat()==FlexAR::CSFASTQ){
										erase(squality,startPos,readLength);
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
					m_adapters->at(i).second += 1;

					//output alignment
					if((m_verbose==FlexAR::CHANGED)||(m_verbose==FlexAR::ALL)){
						ss << endl;
						ss << "read-tag     : " << myRead.getSequenceTag() << endl;
						ss << "adapter-tag	: " << m_adapters->at(i).first->getSequenceTag()<< endl;
						ss << "read        : " << read << endl;
						ss << "quality     : " << quality << endl;

						ss << "beginR       : " << startPosS << endl;
						ss << "endR         : " << endPosS << endl;

						ss << "beginA       : " << startPosA << endl;
						ss << "endA         : " << endPosA << endl;

						ss << alString.str() << endl;
						ss << "overlap(bp) : " << overlapLength << endl;
						ss << "score       : " << alignmentScore << endl;

						ss << "indels        : " << gaps << endl;
						ss << "mismatches    : " << mismatches << endl;
						ss << "allowed errors: " << allowedErrors << endl;
						ss << "cutOff        : " << m_cutOff << endl;

						ss << endl;
						ss << "Adapter removed !" << endl << endl;

						//cout << ss.str() << std::endl;
					}

					if(m_verbose==FlexAR::TAB)ss << myRead.getSequenceTag() << "\t" << startPosA << "\t"<< endPosA << "\t" << overlapLength << "\t" << mismatches << "\t" << gaps << "\t" << allowedErrors << std::endl;

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
						if((myRead.getFormat()==FlexAR::FASTQ)||(myRead.getFormat()==FlexAR::CSFASTQ)){
							ss << "Trimmed quality  : " << myRead.getQuality() << std::endl<< std::endl;
						}
					}

					break;//for loop
				}
				else
				{
					if(m_verbose==FlexAR::ALL){
						ss << endl;
						ss << "read-tag     : " << myRead.getSequenceTag() << endl;
						ss << "adapter-tag	: " << m_adapters->at(i).first->getSequenceTag()<< endl;
						ss << "read        : " << read << endl;
						ss << "quality     : " << quality << endl;

						ss << "beginR       : " << startPosS << endl;
						ss << "endR         : " << endPosS << endl;

						ss << "beginA       : " << startPosA << endl;
						ss << "endA         : " << endPosA << endl;

						ss << alString.str() << endl;
						ss << "overlap(bp) : " << overlapLength << endl;
						ss << "score       : " << alignmentScore << endl;

						ss << "indels        : " << gaps << endl;
						ss << "mismatches    : " << mismatches << endl;
						ss << "allowed errors: " << allowedErrors << endl;
						ss << "cutOff        : " << m_cutOff << endl;
						ss << endl;
						ss << "adapter not aligned to sequence, trim-end criteria missed." << endl;
					}
					//myRead.setDiscard(true);
				}


			}
			else
			{
				if(m_verbose==FlexAR::ALL){
					ss << endl;
					ss << "read-tag     : " << myRead.getSequenceTag() << endl;
					ss << "adapter-tag	: " << m_adapters->at(i).first->getSequenceTag()<< endl;
					ss << "read        : " << read << endl;
					ss << "quality     : " << quality << endl;

					ss << "beginR       : " << startPosS << endl;
					ss << "endR         : " << endPosS << endl;

					ss << "beginA       : " << startPosA << endl;
					ss << "endA         : " << endPosA << endl;

					ss << alString.str() << endl;
					ss << "overlap(bp) : " << overlapLength << endl;
					ss << "score       : " << alignmentScore << endl;

					ss << "indels        : " << gaps << endl;
					ss << "mismatches    : " << mismatches << endl;
					ss << "allowed errors: " << allowedErrors << endl;
					ss << "cutOff        : " << m_cutOff << endl;
					ss << endl;
					ss << "adapter not aligned to sequence, too many mismatches/gaps or overlap-length to low." << endl;
				}

			}
		}

		//Bundeled output for multithreading
		if(m_verbose!=FlexAR::NONE)cout << ss.str();

		return &myRead;

	};



};

#endif /* ALIGNMENTFILTER_H_ */
