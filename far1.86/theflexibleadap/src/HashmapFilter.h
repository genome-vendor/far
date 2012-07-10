/*
 * m_adaptorAlignmentFilter.h
 *
 *  Created on: Apr 13, 2010
 *      Author: mat
 */

#ifndef HashmapFilter_H_
#define HashmapFilter_H_
#include <sstream>
#include "SequencingRead.h"
#include "tbb/pipeline.h"
//#include <tbb/pipeline.h>
#include "tbb/concurrent_hash_map.h"

template <class TString, class TIDString>
class HashmapFilter : public tbb::filter{
private:
	std::string m_target1;
	unsigned long m_singles;
	unsigned int m_suffixIgnore;
	unsigned int m_prefixIgnore;

	// Structure that defines hashing and comparison operations for user's type.
	struct MyHashCompare {
	static size_t hash( const TIDString& x ) {
	//size_t h = 0;

	//Hash function
	size_t h = 0;
	for (unsigned int i = 0; i < seqan::length(x); i++) {
	    h = 31*h + x[i];
	}

	return h;

	}

	//! True if strings are equal
	static bool equal( const TIDString& x, const TIDString& y ) {
		return x==y;
	}

	};
	typedef tbb::concurrent_hash_map<TIDString,SequencingRead<TString,TIDString>*,MyHashCompare> StringTable;
	StringTable m_ids;
public:
	HashmapFilter(std::string target1) : tbb::filter(parallel){
		m_target1 = target1;
		m_singles = 0;
		m_suffixIgnore = 0;
		m_prefixIgnore = 0;
	};

	virtual ~HashmapFilter(){

	};

	void* operator()( void* item ) {

		std::stringstream ss;
		SequencingRead<TString,TIDString> *myRead = static_cast< SequencingRead<TString,TIDString>* >(item);

		if(myRead==NULL)return NULL;

		TString id = myRead->getSequenceTag();
		if(m_suffixIgnore>0){
			id = seqan::prefix(id,length(id)-m_suffixIgnore);
		}

		if(m_prefixIgnore>0){
			id = seqan::suffix(id,m_prefixIgnore);
		}
		//std::cout << "ID: " << id << std::endl;

		typename StringTable::accessor a;
		m_ids.insert( a, id );
		a->second = myRead;
		//ids[id] = myRead;

		return myRead;
	}

	SequencingRead<TString,TIDString>* hasHashKey(TIDString key){
		typename StringTable::accessor a;
		if(m_ids.find(a, key)){
			return a->second;
		}
		return NULL;
	}

	void writeOmittedReads(){
		std::fstream m_omitted1;
		std::string ofilename = m_target1 + ".single";
		m_omitted1.open(ofilename.c_str(), std::ios_base::out | std::ios_base::binary);
		if(!m_omitted1.is_open())
		{
			std::cout << "Error opening File: " << ofilename << std::endl;
		}

		typename StringTable::iterator i;
		for( i = m_ids.begin(); i!=m_ids.end(); ++i ){
			//std::cout << i->second->getSequenceTag() << "bbb" <<i->second->getDiscard()<<std::endl;
			if(i->second->getDiscard()==false){
				++m_singles;
				m_omitted1 << i->second->getFastString();
			}
		}
	}

	unsigned long getNrSingles(){
		return m_singles;
	}

	void setSuffixIgnore(unsigned int nr){
		this->m_suffixIgnore = nr;
	}
	void setPrefixIgnore(unsigned int nr){
		this->m_prefixIgnore = nr;
	}
};

#endif
