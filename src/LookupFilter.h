/*
 * m_adaptorAlignmentFilter.h
 *
 *  Created on: Apr 13, 2010
 *      Author: mat
 */

#ifndef LookupFilter_H_
#define LookupFilter_H_
#include <sstream>
#include "SequencingRead.h"
#include "tbb/pipeline.h"
#include "HashmapFilter.h"

/**This class is used by the PairedreadFinder. It will check for each processed read if it was stored by the HashmapFilter (if ID's match). */
template <class TString, class TIDString>
class LookupFilter : public tbb::filter{
private:
	std::string m_target1;
	std::string m_target2;
	std::fstream m_t1;
	std::fstream m_t2;
	std::fstream m_omitted2;
	unsigned long m_paired, m_single2;
	unsigned int m_suffixIgnore;
	unsigned int m_prefixIgnore;

	HashmapFilter<TString, TIDString> *m_filter;
public:
	LookupFilter(HashmapFilter<TString, TIDString> *filter,std::string target1, std::string target2) : tbb::filter(parallel){
		m_filter = filter;
		m_target1 = target1;
		m_target2 = target2;
		m_paired = 0;
		m_single2 = 0;
		m_suffixIgnore = 0;
		m_prefixIgnore = 0;

		m_t1.open(m_target1.c_str(), std::ios_base::out | std::ios_base::binary);
		if(!m_t1.is_open())
		{
			std::cout << "Error opening File: " << m_target1 << std::endl;
		}

		m_t2.open(m_target2.c_str(), std::ios_base::out | std::ios_base::binary);
		if(!m_t2.is_open())
		{
			std::cout << "Error opening File: " << m_target2 << std::endl;
		}

		std::string ofilename = m_target2 + ".single";
		m_omitted2.open(ofilename.c_str(), std::ios_base::out | std::ios_base::binary);
		if(!m_omitted2.is_open())
		{
			std::cout << "Error opening File: " << ofilename << std::endl;
		}

	}

	unsigned long getNrPairs(){
		return m_paired;
	}

	unsigned long getNrSingles(){
		return m_single2;
	}

	virtual ~LookupFilter(){

	}

	void setSuffixIgnore(unsigned int nr){
			this->m_suffixIgnore = nr;
	}

	void setPrefixIgnore(unsigned int nr){
		this->m_prefixIgnore = nr;
	}

	void* operator()( void* item ) {

		std::stringstream ss;
		SequencingRead<TString,TIDString> *myRead = static_cast< SequencingRead<TString,TIDString>* >(item), *hashRead;

		if(myRead==NULL){
			return NULL;
		}

		TString id = myRead->getSequenceTag();
		if(m_suffixIgnore>0){
			id = seqan::prefix(id,length(id)-m_suffixIgnore);
		}

		if(m_prefixIgnore>0){
			id = seqan::suffix(id,m_prefixIgnore);
		}

		//std::cout << "checking " << id << std::endl;
		hashRead = m_filter->hasHashKey(id);
		//std::cout << "done " << std::endl;
		if(hashRead!=NULL){
			//id was found
			//std::cout << "Read " << hashRead->getSequenceTag() << " " << myRead->getDiscard() << std::endl;
			m_t1 << hashRead->getFastString();
			m_t2 << myRead->getFastString();
			++m_paired;

			//reads which were referenced are marked (default false)
			hashRead->setDiscard(true);
		}
		else{
			//id was not found in hash table
			m_omitted2 << myRead->getFastString();
			++m_single2;
		}

		delete myRead;
		//delete hashRead;


		return NULL;
	}

};

#endif
