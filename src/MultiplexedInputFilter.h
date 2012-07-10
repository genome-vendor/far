/*
 * MultiplexInputFilter.h
 *
 *  Created on: May 31, 2011
 *      Author: mat
 */

#ifndef MULTIPLEXINPUTFILTER_H_
#define MULTIPLEXINPUTFILTER_H_

#include "MultiplexedRead.h"
#include "SequenceInputFilter.h"

/**This class will handle up to 3 file sources at the same time (paired read input plus barcode reads) and create a MultiplexedRead
depending on the run type (single-end, paired-end and/or barcoded).*/
template <typename TString, typename TIDString>
class MultiplexedInputFilter : public tbb::filter {
private:
	bool m_isPaired;
	bool m_isBarcoded;
	int m_allowedUncalledBases;
	int m_pre_cut_length;
	int m_phred_pre_cut_qual;
	int m_min_read_length;
	::FlexAR::QualityType m_qual;
	std::string m_filePath1;
	FlexAR::FileFormat m_format;
	SequenceInputFilter<TString,TString > *m_f1;
	SequenceInputFilter<TString,TString > *m_f2;
	SequenceInputFilter<TString,TString > *m_b;
	long m_cnt_uncalled;
	long m_cnt_total;
	long m_cnt_invalid;
public:
	MultiplexedInputFilter(std::string filePath1, FlexAR::FileFormat format ,int allowedUncalledBases, int pre_cut_length, int min_read_length, int phred_pre_cut_qual,::FlexAR::QualityType qual) : filter(serial_in_order){
		m_filePath1 = filePath1;
		m_cnt_uncalled = 0;
		m_cnt_total = 0;
		m_cnt_invalid = 0;
		m_format = format;
		m_allowedUncalledBases = allowedUncalledBases;
		m_pre_cut_length = pre_cut_length;
		m_phred_pre_cut_qual = phred_pre_cut_qual;
		m_min_read_length = min_read_length;
		m_f1 = new SequenceInputFilter<TString, TString>( filePath1, m_format );
		m_f1->setFixedPreTrim(pre_cut_length);
		m_f1->SetMinReadLength(min_read_length)	;
		m_f1->setPrePhredTrim(phred_pre_cut_qual,qual);
		m_qual = qual;
		

		m_isPaired = false;
		m_isBarcoded = false;

		m_f2 = NULL;
		m_b = NULL;
	}

	virtual ~MultiplexedInputFilter(){

	}

	void setPairedFile(std::string filePath2){
		m_f2 = new SequenceInputFilter<TString, TString>( filePath2, m_format );
		m_f2->setFixedPreTrim(m_pre_cut_length);
		m_f2->SetMinReadLength(m_min_read_length);
		m_f2->setPrePhredTrim(m_phred_pre_cut_qual,m_qual);
		m_isPaired = true;

	}

	void setBarcodeReadsFile(std::string bfileName){
		m_b = new SequenceInputFilter<TString, TString>( bfileName, m_format );
		m_b->setFixedPreTrim(0);
		m_isBarcoded = true;
	}

	void* operator()(void*) {
		SequencingRead<TString, TIDString> *myRead1=NULL,*myRead2=NULL,*myBarcode=NULL;
		MultiplexedRead<TString, TIDString> *mRead = NULL;

		bool found = false;
		//Skip reads that have uncalled bases in r1 or r2
		do{
			//look for a new read
			//single read input
			if(!m_isPaired)
			{
				//this will skip all reads that have not an adequate readlength
				myRead1 = getValidRead(m_f1);
				if(myRead1==NULL)return NULL;
			}
			else {
			//paired read input
				//both reads have to be valid (>minlength). if not we will process a single read
				bool v1=false,v2=false;

				//find at least one valid read
				while(v1!=true || v2!=true)
				{
					myRead1 = static_cast< SequencingRead<TString,TIDString>* >(m_f1->getRead(v1));
					myRead2 = static_cast< SequencingRead<TString,TIDString>* >(m_f2->getRead(v2));
					//if both reads are NULL we reached EOF
					if(myRead1 == NULL && myRead2 == NULL)return NULL;
				}

				if(!v1)myRead1 = NULL;
				if(!v2)myRead2 = NULL;				
			}
			
			///just read barcode and don't worry if it's valid...
			bool valid = false;
			if(m_isBarcoded)myBarcode = static_cast< SequencingRead<TString,TIDString>* >(m_b->getRead(valid));

			//paired end barcoded
			if((myRead1!=NULL)&&(myRead2!=NULL)&&(myBarcode!=NULL)){
				++m_cnt_total;
				++m_cnt_total;
				if(!((myRead1->isUncalledSequence(m_allowedUncalledBases))||(myRead2->isUncalledSequence(m_allowedUncalledBases)))){
					mRead = new MultiplexedRead<TString, TIDString>(myRead1,myRead2,myBarcode);
					return mRead;
				}
				else found=false;
			}

			//paired end
			if((myRead1!=NULL)&&(myRead2!=NULL)&&(myBarcode==NULL)){
				++m_cnt_total;
				++m_cnt_total;
				if(!((myRead1->isUncalledSequence(m_allowedUncalledBases))||(myRead2->isUncalledSequence(m_allowedUncalledBases)))){
					mRead = new MultiplexedRead<TString, TIDString>(myRead1,myRead2);
					return mRead;
				}
				else found=false;
			}

			//single end
			if((myRead1!=NULL)&&(myRead2==NULL)&&(myBarcode==NULL)){
				++m_cnt_total;
				if(!((myRead1->isUncalledSequence(m_allowedUncalledBases)))){
					mRead = new MultiplexedRead<TString, TIDString>(myRead1);
					return mRead;
				}
				else found=false;
			}

			//single end barcoded
			if((myRead1!=NULL)&&(myRead2==NULL)&&(myBarcode!=NULL)){
				++m_cnt_total;
				myRead2 = NULL;
				if(!((myRead1->isUncalledSequence(m_allowedUncalledBases)))){
					//mRead = new MultiplexedRead<TString, TIDString>(myRead1,NULL,myBarcode);
					mRead = new MultiplexedRead<TString, TIDString>(myRead1,myRead2,myBarcode);
					return mRead;
				}
				else found=false;
			}

			//we didn't return any read, so look for the next valid one

			++m_cnt_uncalled;


		}
		while(!found);

		return mRead;

	}

	SequencingRead<TString,TIDString>* getValidRead(SequenceInputFilter<TString,TString > *input)
	{
		SequencingRead<TString,TIDString>* myRead1 =  NULL;
		bool valid = false;
		--m_cnt_invalid;
		while(!valid){
			++m_cnt_invalid;
			myRead1 = static_cast< SequencingRead<TString,TIDString>* >(m_f1->getRead(valid));
			if(myRead1==NULL)break;
		}


		return myRead1;
	}

	long getNrTotalReads(){
		return m_cnt_total;
	}

	long getNrUncalledReads(){
		return m_cnt_uncalled;
	}


	long getNrProcessedReads(){
		if(m_isPaired)return m_f1->getNrProcessedReads() + m_f2->getNrProcessedReads();
		else return m_f1->getNrProcessedReads();
	}

	long getNrLowPhredDiscarded(){
		if(m_isPaired)return m_f1->getNrLowPhredDiscarded() + m_f2->getNrLowPhredDiscarded();
		else return m_f1->getNrLowPhredDiscarded();
	}

};

#endif /* MULTIPLEXINPUTFILTER_H_ */
