/*
 * SequenceInputFilter.h
 *
 *  Created on: Apr 12, 2010
 *      Author: mat
 */

#ifndef SEQUENCEINPUTFILTER_H_
#define SEQUENCEINPUTFILTER_H_
#include <string>
#include <fstream>
#include "SequencingRead.h"
#include <tbb/pipeline.h>
#include "Enums.h"
#include <stdexcept>

template <class TString, class TIDString>
class SequenceInputFilter : public tbb::filter {

private:
    std::fstream fstrm;
	std::string m_filePath;
	FlexAR::FileFormat m_format;
	bool m_discard;
	int m_allowed_uncalled;
	unsigned int m_fixedPreTrim;
	int m_fixedPrePhredTrim;
	unsigned long m_count_uncalled;
	unsigned long m_count_lowphred;
	unsigned long m_nrReads;
	unsigned int m_minLength;
	TIDString m_nextTag;

	bool isUncalledSequence(TString sequence){
		int n = 0;

		typename seqan::Iterator<TString >::Type it = seqan::begin(sequence);
		typename seqan::Iterator<TString >::Type itEnd = seqan::end(sequence);
		while (it != itEnd) {
			 if((*it=='.')||(*it=='N'))++n;
			 ++it;
		}

		if(n <= m_allowed_uncalled)return false;
		return true;
	}

public:
	void setNrAllowedUncalledBases(int nr){
		m_allowed_uncalled = nr;
	}

	void setFixedPreTrim(int nr){
		if(nr > 0)this->m_fixedPreTrim = nr;
	}

	void setPrePhredTrim(int nr,::FlexAR::QualityType qtype){
		if(nr > 0)this->m_fixedPrePhredTrim = nr;
		else {
			std::cout << "Using no phred-quality trimming." << std::endl;
			return;
		}


		switch(qtype){
			case FlexAR::SANGER:{
				this->m_fixedPrePhredTrim += 33;
				std::cout << "Trimming reads from 3' end to phred quality " << nr << " ( value in SANGER: " << m_fixedPrePhredTrim << " )" << std::endl;
				break;
			}
			case FlexAR::SOLEXA:{
				this->m_fixedPrePhredTrim += 59;
				std::cout << "Trimming reads from 3' end to phred quality " << nr << " ( value in SOLEXA: " << m_fixedPrePhredTrim << " )" << std::endl;
				break;
			}
			case FlexAR::ILLUMINA13:{
				this->m_fixedPrePhredTrim += 64;
				std::cout << "Trimming reads from 3' end to phred quality " << nr << " ( value in illumina13: " << m_fixedPrePhredTrim << " )" << std::endl;
				break;
			}
			case FlexAR::ILLUMINA15:{
				this->m_fixedPrePhredTrim += 66;
				std::cout << "Trimming reads from 3' end to phred quality " << nr << " ( value in illumina15: " << m_fixedPrePhredTrim << " )" << std::endl;
				break;
			}
		}




	}

	void SetMinReadLength(unsigned int readMinLength){
		m_minLength = readMinLength;
	}

	unsigned long getNrUncalledDiscarded(){
		return m_count_uncalled;
	}

	unsigned long getNrLowPhredDiscarded(){
		return m_count_lowphred;
	}
    //static const size_t n_buffer = 8;
	SequenceInputFilter(std::string filePath, FlexAR::FileFormat format ) : filter(serial_in_order)
	{
		m_allowed_uncalled = -1;
		m_discard = true;
		m_filePath = filePath;
		m_format = format;
		m_nextTag = "";
		m_nrReads = 0;
		m_count_uncalled = 0;
		m_fixedPreTrim = 0;
		m_fixedPrePhredTrim = -1;

		fstrm.open(m_filePath.c_str(), std::ios_base::in);
		if(!fstrm.is_open())
		{
			std::cout << "Error opening File: " << m_filePath << std::endl;
			exit(0);
		}

	};

	void setDiscradUncalledReads(bool discard){
		m_discard = discard;
	}

	unsigned long getNrProcessedReads(){
		return m_nrReads;
	}

	TString readLine(){
		char line[4096];
		TString text = "";

		if(fstrm.good()){
			fstrm.getline(line,4096);

			//ignore comment lines
			//while(fstrm.good() && fstrm.gcount() > 0 && line[0] == ';'){
			//	fstrm.getline(line,4096);
			//}

			text = line;
		}
		return text;
	}

	//override
	void* operator()(void*) {
		SequencingRead<TString, TIDString> *myRead;

		TString source="",quality="",dummy="";
		TIDString sequence_tag="";

		if(!fstrm.eof()){
			try{
				//if(found)exit(0);
				if((m_format == FlexAR::FASTA)||(m_format == FlexAR::CSFASTA)){
					bool search = true;
					while(search){
						//FastA parsing
						//tag line will be read in previous iteration
						if(m_nextTag=="")sequence_tag = readLine();
						else sequence_tag = m_nextTag;

						if(seqan::length(sequence_tag) > 0){
							if(seqan::isNotEqual(getValue(sequence_tag,0), '>')){
								std::stringstream error;
								error << "Incorrect FASTA entry, missing new > line. Input: " << sequence_tag << std::endl;
								throw std::runtime_error(error.str());

							}
							else sequence_tag = seqan::suffix(sequence_tag,1);

							if(seqan::length(sequence_tag)==0){
								std::stringstream error;
								error << "Incorrect FASTA entry, missing readname after >. Input: " << sequence_tag << std::endl;
								throw std::runtime_error(error.str());
							}


						}
						else return NULL;

						source = readLine();

						if(seqan::length(source) < 1){
							std::stringstream error;
							error << "Warning, found sequence tag without read! Tag: " << sequence_tag << std::endl;
							throw std::runtime_error(error.str());
						}

						m_nextTag = readLine();

						//For fasta files wich have sequences splitted over several lines
						while(fstrm.good() && seqan::isNotEqual(getValue(m_nextTag,0), '>')){
							append(source,m_nextTag);
							m_nextTag = readLine();
						}

						if(m_allowed_uncalled!=-1){
							if(!isUncalledSequence(source))search = false;
							else ++m_count_uncalled;
						}
						else search=false;
					}

					if((this->m_fixedPreTrim > 1)&&(seqan::length(source) > m_fixedPreTrim)){
						source = seqan::prefix(source,m_fixedPreTrim);
					}

					myRead = new SequencingRead<TString,TIDString>(source,sequence_tag,m_format);
					++m_nrReads;

				}
				else{
					bool search = true;
					while(search){
						//FastQ parsing
						source = readLine();

						if(seqan::length(source) > 0){
							if(seqan::isNotEqual(getValue(source,0), '@')){
								std::stringstream error;
								error << "Incorrect FASTQ entry, missing new @ line. Input: " << source << std::endl;
								throw std::runtime_error(error.str());

							}
							else sequence_tag = seqan::suffix(source,1);

							if(seqan::length(sequence_tag)==0){
								std::stringstream error;
								error << "Incorrect FASTQ entry, missing readname after @. Input: " << source << std::endl;
								throw std::runtime_error(error.str());
							}
						}
						else return NULL;

						source = readLine();
						if(seqan::length(source) < 1){
							std::stringstream error;
							error << "Warning, found sequence tag without read! Tag: " << sequence_tag << std::endl;
							throw std::runtime_error(error.str());

						}

						dummy = readLine();
						if((seqan::length(dummy) == 0)||(seqan::isNotEqual(getValue(dummy,0), '+'))){
								std::stringstream error;
								error << "Incorrect FASTQ entry, missing + line. Readname: " << sequence_tag << std::endl;
								throw std::runtime_error(error.str());

						}

						quality = readLine();

						//In case CSFASTQ format has same quality and readlength it will be trimmed accordingly
						if(m_format==FlexAR::CSFASTQ){
							if(length(quality) == length(source)){
								quality = seqan::suffix(quality,1);
							}
						}

						if(seqan::length(source) < 1){
							std::stringstream error;
							error << "Warning, found sequence without quality values! Tag: " << sequence_tag << std::endl;
							throw std::runtime_error(error.str());

						}

						if(m_allowed_uncalled!=-1){
							if(!isUncalledSequence(source))search = false;
							else ++m_count_uncalled;
						}
						else search=false;

						//...and/or filter due to phred quality
						if(this->m_fixedPrePhredTrim!=-1){
							typename seqan::Iterator<TString >::Type it = seqan::begin(quality);
							typename seqan::Iterator<TString >::Type itEnd = seqan::end(quality);
							--itEnd;
							unsigned int n = length(quality);
							//std::cout << "QSTRING: " << quality << std::endl;
							//std::cout << m_fixedPrePhredTrim << std::endl;
							while (itEnd != it) {
								 if(static_cast<unsigned int>(*itEnd)>=this->m_fixedPrePhredTrim)break;
								 //std::cout << "QUAL: " << static_cast<unsigned int>(*itEnd) << " CHAR:" << (*itEnd) <<std::endl;
								 --n;
								 --itEnd;
							}

							if(n>this->m_minLength){
								source = seqan::prefix(source,n);
								search = false;
								if(m_format == FlexAR::CSFASTQ){
									quality = seqan::prefix(quality,n - 1);
								}
								else quality = seqan::prefix(quality,n);
	//std::cout << "TRIMMED READ " << sequence_tag << "-" << n << "," << m_minLength << std::endl;
							}
							else {
								search = true;
								++m_count_lowphred;
							}

						}
					}

					//now pre-cut read

					if((this->m_fixedPreTrim > 1)&&(seqan::length(source) > m_fixedPreTrim)){
						source = seqan::prefix(source,m_fixedPreTrim);

						if(m_format == FlexAR::CSFASTQ){
							quality = seqan::prefix(quality,m_fixedPreTrim - 1);
						}
						else quality = seqan::prefix(quality,m_fixedPreTrim);
					}



					myRead = new SequencingRead<TString,TIDString>(source,sequence_tag,quality,m_format);

					++m_nrReads;
				}

				return myRead;

			}
			catch(std::ios_base::failure &failure){
				std::cout << failure.what() << std::endl;
				fstrm.close();
				return NULL;
			}

		}
		else {
			return NULL;
		}
	}

	virtual ~SequenceInputFilter(){};


};

#endif /* SEQUENCEINPUTFILTER_H_ */
