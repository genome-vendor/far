/*
 * Read.h
 *
 *  Created on: Apr 12, 2010
 *      Author: mat
 */

#ifndef READ_H_
#define READ_H_
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include "Enums.h"
template <class TString,class TIDString>
class SequencingRead {
private:
	TString m_source;
	TIDString m_sequence_tag;
	TIDString m_qual;
	FlexAR::FileFormat m_format;
	bool m_discard;
public:
	SequencingRead(){
		m_source = "";
		m_sequence_tag = "";
		m_discard = false;
	}

	SequencingRead(TString source, TIDString sequence_tag, FlexAR::FileFormat format){
		m_sequence_tag = sequence_tag;
		m_format = format;
		m_discard = false;
		m_source = source;
	}

	SequencingRead(TString source, TIDString sequence_tag, TIDString qual, FlexAR::FileFormat format){
		m_sequence_tag = sequence_tag;
		m_format = format;
		m_discard = false;
		m_source = source;
		m_qual = qual;
	}

	TString getSequence(){
		return m_source;
	}

	void setSequence(TString seq){
		m_source = seq;
	}

	void setSequenceTag(TString tag){
		m_sequence_tag = tag;
	}


	void setQuality(TString qual){
		m_qual = qual;
	}

	TIDString getQuality(){
		return m_qual;
	}

	void setDiscard(bool discard){
		m_discard = discard;
	}


	bool getDiscard(){
		return m_discard;
	}

	TIDString getSequenceTag(){
		return m_sequence_tag;
	}

	FlexAR::FileFormat getFormat(){
		return m_format;
	}

	std::string getFastString(){
		std::stringstream ss;
			switch(m_format){
				case FlexAR::FASTQ:ss << "@" << this->getSequenceTag() << std::endl << this->getSequence() << std::endl << "+" << std::endl << this->getQuality() << std::endl;
					break;
				case FlexAR::FASTA:ss <<  ">" << this->getSequenceTag() << std::endl << this->getSequence() << std::endl;
					break;
				case FlexAR::CSFASTQ:ss << "@" << this->getSequenceTag() << std::endl <<  this->getSequence() << std::endl << "+" << std::endl << this->getQuality() << std::endl;
					break;
				case FlexAR::CSFASTA:ss <<  ">" << this->getSequenceTag() << std::endl <<  this->getSequence() << std::endl;
					break;
			}

			return ss.str();
	}

	virtual ~SequencingRead(){};
};

#endif /* READ_H_ */
