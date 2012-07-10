/*
 * SequenceOutputFilter.h
 *
 *  Created on: Apr 12, 2010
 *      Author: mat
 */

#ifndef SEQUENCEOUTPUTFILTER_H_
#define SEQUENCEOUTPUTFILTER_H_
#include "SequencingRead.h"
#include <fstream>
#include <tbb/pipeline.h>
#include "Enums.h"
#include <map>

/**This class will write a SequencingRead in FASTA/FASTQ format (depending on the FlexAR::FileFormat) to a file specified in the constructor. */
template <typename TString,typename TIDString>
class SequenceOutputFilter {
private:
	std::map<unsigned int,unsigned int> m_lengthDist;
	std::fstream m_targetStream;
	std::string m_filePath;
	unsigned int m_fixedPostTrim;
	unsigned long m_countGood;
	FlexAR::FileFormat m_format;

public:
	void writeLengthDist(){
		std::string fname = m_filePath + ".lengthdist";
		std::fstream lstream;
		lstream.open(fname.c_str(), std::ios_base::out | std::ios_base::binary);
		if(!lstream.is_open())
		{
			std::cerr << "Error opening File: " << fname << std::endl;
		}
		else {
			std::map<unsigned int,unsigned int>::iterator iter;
			lstream << "Readlength\tCount" << std::endl;
			for( iter = m_lengthDist.begin(); iter != m_lengthDist.end(); iter++ ) {
				lstream << iter->first << "\t" << iter->second << std::endl;
			}
			lstream.close();
		}

	}

	std::string getFastString(SequencingRead<TString,TIDString> *myRead){

		std::stringstream ss;
		switch(m_format){
			case FlexAR::FASTQ:ss << "@" << myRead->getSequenceTag() << std::endl << myRead->getSequence() << std::endl << "+" << std::endl << myRead->getQuality() << std::endl;
				break;
			case FlexAR::FASTA:ss <<  ">" << myRead->getSequenceTag() << std::endl << myRead->getSequence() << std::endl;
				break;
			case FlexAR::CSFASTQ:ss << "@" << myRead->getSequenceTag() << std::endl <<  myRead->getSequence() << std::endl << "+" << std::endl << myRead->getQuality() << std::endl;
				break;
			case FlexAR::CSFASTA:ss <<  ">" << myRead->getSequenceTag() << std::endl <<  myRead->getSequence() << std::endl;
				break;
		}

		return ss.str();
	}

	void setFixedPostTrim(int nr){
		if(nr > 0)this->m_fixedPostTrim = nr;
	}

	unsigned long getNrGoodReads(){
		return m_countGood;
	}

	SequenceOutputFilter(std::string filePath, FlexAR::FileFormat format)// : filter(serial_in_order)
	{
		m_format = format;
		m_countGood = 0;
		m_fixedPostTrim = 0;
		m_filePath = filePath;
		
		m_targetStream.open(m_filePath.c_str(), std::ios_base::out | std::ios_base::binary);
		if(!m_targetStream.is_open())
		{
			std::cerr << "Error opening File: " << m_filePath << std::endl;
		}

	};

	std::string getFileName(){
		return m_filePath;
	}

	virtual ~SequenceOutputFilter(){
		m_targetStream.close();
	};

	void *writeRead(void* item){
		if(item)
		{
			SequencingRead<TString,TIDString> *myRead = static_cast< SequencingRead<TString,TIDString>* >(item);

			unsigned int readLength = length(myRead->getSequence());
			
			if(m_targetStream.is_open() && m_targetStream.good())
			{
				++m_countGood;
				if(m_fixedPostTrim > 0){
					if(readLength>m_fixedPostTrim){
						myRead->setSequence(prefix(myRead->getSequence(),m_fixedPostTrim));
						if(myRead->getFormat()==FlexAR::CSFASTQ){
							myRead->setQuality(prefix(myRead->getQuality(),m_fixedPostTrim - 1));
						}
						if(myRead->getFormat()==FlexAR::FASTQ){
							myRead->setQuality(prefix(myRead->getQuality(),m_fixedPostTrim));
						}

						readLength = m_fixedPostTrim;
					}
				}

				//store read length distribution
				std::map<unsigned int, unsigned int>::iterator it;
				it = m_lengthDist.find(readLength);

				if(it!=m_lengthDist.end()){
					it->second = it->second + 1;
				}
				else {
					m_lengthDist[readLength] = 1;
				}
				m_targetStream << getFastString(myRead);
				m_targetStream.flush();
				return NULL;
			}
			else {
				std::cerr << "Error writing target file " << m_filePath << std::endl;
				exit(1);
			}

		}
		return NULL;
	}
};

#endif /* SEQUENCEOUTPUTFILTER_H_ */
