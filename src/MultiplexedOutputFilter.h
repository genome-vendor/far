/*
 * MultiplexInputFilter.h
 *
 *  Created on: May 31, 2011
 *      Author: mat
 */

#ifndef MULTIPLEXOUTPUTFILTER_H_
#define MULTIPLEXOUTPUTFILTER_H_

#include "MultiplexedRead.h"
#include "SequenceOutputFilter.h"
#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>
#include "Enums.h"
#include "OutputFileStruct.h"


/**This class will process a MultiplexedRead and write it to a file depending on the runtype (single-end, paired-end and/or barcoded). */
template <typename TString, typename TIDString>
class MultiplexedOutputFilter : public tbb::filter {
private:
	std::string m_filePath;
	FlexAR::FileFormat m_format;
	FlexAR::RunType m_runType;
	typedef ::std::pair<SequencingRead<seqan::CharString,seqan::CharString>*,unsigned int> TAdapter;
	tbb::concurrent_vector<TAdapter> *m_barcodes;
	typedef SequenceOutputFilter<TString,TIDString> TOutputFilter;
	unsigned int m_minLength;

	typedef OutputFileStruct<TString,TIDString> filters;

	typedef std::map<int, filters > TMap;
	filters *m_outputMap;
	unsigned int m_mapsize;


	long m_cnt_total;
public:
	MultiplexedOutputFilter(std::string filePath, tbb::concurrent_vector<TAdapter> *barcodes, FlexAR::FileFormat format, int min_read_length, std::string omittedFilename, int post_cut_length, ::FlexAR::RunType runType) : filter(serial_in_order){
		m_filePath = filePath;
		m_format = format;
		m_barcodes = barcodes;
		m_runType = runType;
		m_minLength = min_read_length;
		m_mapsize = 0;

		std::stringstream ss;
		switch(m_runType){
			case FlexAR::PAIRED_BARCODED:{
				m_mapsize = barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				for(unsigned int i=0;i<barcodes->size();++i)
				{
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_1" << toFormatString(m_format);
					TOutputFilter *filter1 = new TOutputFilter(ss.str(),m_format);

					ss.str("");
					ss.clear();

					filter1->setFixedPostTrim(post_cut_length);

					//single
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_1_single" << toFormatString(m_format);
					TOutputFilter *single1 = new TOutputFilter(ss.str(),m_format);

					ss.str("");
					ss.clear();

					single1->setFixedPostTrim(post_cut_length);

					ss.str("");
					ss.clear();

					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_2"<< toFormatString(m_format);
					TOutputFilter *filter2 = new TOutputFilter(ss.str(),m_format);

					ss.str("");
					ss.clear();

					filter2->setFixedPostTrim(post_cut_length);

					//single reads
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_2_single"<< toFormatString(m_format);
					TOutputFilter *single2 = new TOutputFilter(ss.str(),m_format);

					ss.str("");
					ss.clear();

					single2->setFixedPostTrim(post_cut_length);
					ss.str("");
					ss.clear();

					filters f;
					f.f1 = filter1;
					f.f2 = filter2;
					f.single1 = single1;
					f.single2 = single2;

					m_outputMap[i + 1] = f;
				}
				ss << filePath << "_barcode_unassigned_1"<< toFormatString(m_format);
				TOutputFilter *filter3 = new SequenceOutputFilter<TString,TIDString>(ss.str(),m_format);

				ss.str("");
				ss.clear();

				ss << omittedFilename << "_barcode_unassigned_1"<< toFormatString(m_format);
				filter3->setFixedPostTrim(post_cut_length);

				ss.str("");
				ss.clear();

				//single reads
				ss << filePath << "_barcode_unassigned_1_single"<< toFormatString(m_format);
				TOutputFilter *single3 = new SequenceOutputFilter<TString,TIDString>(ss.str(),m_format);

				ss.str("");
				ss.clear();

				ss << omittedFilename << "_barcode_unassigned_1_single"<< toFormatString(m_format);
				single3->setFixedPostTrim(post_cut_length);

				ss.str("");
				ss.clear();


				ss << filePath << "_barcode_unassigned_2"<< toFormatString(m_format);
				TOutputFilter *filter4 = new SequenceOutputFilter<TString,TIDString>(ss.str(),m_format);

				ss.str("");
				ss.clear();

				filter4->setFixedPostTrim(post_cut_length);

				//single reads
				ss << filePath << "_barcode_unassigned_2_single"<< toFormatString(m_format);
				TOutputFilter *single4 = new SequenceOutputFilter<TString,TIDString>(ss.str(),m_format);

				ss.str("");
				ss.clear();

				single4->setFixedPostTrim(post_cut_length);


				filters f2;
				f2.f1 = filter3;
				f2.f2 = filter4;
				f2.single1 = single3;
				f2.single2 = single4;
				m_outputMap[0] = f2;

				break;
			}
			case FlexAR::PAIRED:{//if no barcodes have been specified just write paired output...
				m_mapsize = 1;
				m_outputMap = new filters[m_mapsize];
				ss << filePath << "_1"<< toFormatString(m_format);
				TOutputFilter *filter5 = new TOutputFilter(ss.str(),m_format);

				filter5->setFixedPostTrim(post_cut_length);

				ss.str("");
				ss.clear();

				//single
				ss << filePath << "_1_single" << toFormatString(m_format);
				TOutputFilter *single5 = new TOutputFilter(ss.str(),m_format);

				ss.str("");
				ss.clear();

				ss << omittedFilename << "_1"<< toFormatString(m_format);
				single5->setFixedPostTrim(post_cut_length);

				ss.str("");
				ss.clear();

				ss << filePath << "_2"<< toFormatString(m_format);
				TOutputFilter *filter6 = new TOutputFilter(ss.str(),m_format);

				ss.str("");
				ss.clear();

				filter6->setFixedPostTrim(post_cut_length);

				ss.str("");
				ss.clear();

				//single
				ss << filePath << "_2_single"<< toFormatString(m_format);
				TOutputFilter *single6 = new TOutputFilter(ss.str(),m_format);

				ss.str("");
				ss.clear();

				single6->setFixedPostTrim(post_cut_length);

				ss.str("");
				ss.clear();

				filters f3;
				f3.f1 = filter5;
				f3.f2 = filter6;

				f3.single1 = single5;
				f3.single2 = single6;

				m_outputMap[0] = f3;
				break;
			}
			case FlexAR::SINGLE:{
				m_mapsize = 1;
				m_outputMap = new filters[m_mapsize];
				std::stringstream ss;
				ss << filePath << toFormatString(m_format);
				TOutputFilter *filter7 = new SequenceOutputFilter<TString,TIDString>(ss.str(),m_format);

				ss.str("");
				ss.clear();

				ss << omittedFilename;
				filter7->setFixedPostTrim(post_cut_length);

				filters f4;
				f4.f1 = filter7;
				f4.f2 = NULL;
				m_outputMap[0] = f4;
				break;
			}
			case FlexAR::SINGLE_BARCODED:
			case FlexAR::SINGLE_BARCODED_WITHIN_READ:
				{
				m_mapsize = m_mapsize = barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				std::stringstream ss;
				for(unsigned int i=0;i<barcodes->size();++i)
				{
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << toFormatString(m_format);
					TOutputFilter *filter8 = new TOutputFilter(ss.str(),m_format);

					ss.str("");
					ss.clear();

					filter8->setFixedPostTrim(post_cut_length);

					filters f5;
					f5.f1 = filter8;
					f5.f2 = NULL;
					f5.single1 = NULL;
					f5.single2 = NULL;

					m_outputMap[i + 1] = f5;
				}

				ss.str("");
				ss.clear();

				ss << filePath << "_barcode_unassigned"<< toFormatString(m_format);
				TOutputFilter *filter9 = new SequenceOutputFilter<TString,TIDString>(ss.str(),m_format);

				ss.str("");
				ss.clear();

				filter9->setFixedPostTrim(post_cut_length);

				filters f6;
				f6.f1 = filter9;
				f6.f2 = NULL;
				m_outputMap[0] = f6;
				break;
			}
		}
	}
	virtual ~MultiplexedOutputFilter(){
		//delete ...
		/*TMapIterator it;
		for(it=m_outputMap.begin();it!=m_outputMap.end();++it){
			delete (*it).second.f1;
			delete (*it).second.f2;
		}*/



	}

	void setMinReadlength(int minLength){
		m_minLength = minLength;
	}

	unsigned long getNrGoodReads(){//omitted reads
		long good = 0;
		for(unsigned int i = 0;i < m_mapsize; ++i){
			good += m_outputMap[i].f1->getNrGoodReads();
			if(m_outputMap[i].f2!=NULL)good += m_outputMap[i].f2->getNrGoodReads();
		}
		return good;
	}

	void* operator()(void* item) {
		MultiplexedRead<TString,TIDString> *myRead = static_cast< MultiplexedRead<TString,TIDString>* >(item);
		bool l1ok=false,l2ok=false;
		switch(m_runType){
			case FlexAR::PAIRED_BARCODED:
			case FlexAR::PAIRED:{
			if((myRead->m_r1!=NULL)&&(myRead->m_r2!=NULL)){
				//now check if both reads have minLength
				if(length(myRead->m_r1->getSequence())>=m_minLength)l1ok = true;
				if(length(myRead->m_r2->getSequence())>=m_minLength)l2ok = true;
				if(l1ok&&l2ok){
					m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
					m_outputMap[myRead->m_barcode_id].f2->writeRead(myRead->m_r2);
				}
				if(l1ok&&!l2ok){
					m_outputMap[myRead->m_barcode_id].single1->writeRead(myRead->m_r1);
				}
				if(!l1ok&&l2ok){
					m_outputMap[myRead->m_barcode_id].single2->writeRead(myRead->m_r2);
				}

				if(!l1ok)m_outputMap[myRead->m_barcode_id].m_cnt_short_1+=1;
				if(!l2ok)m_outputMap[myRead->m_barcode_id].m_cnt_short_2+=1;

			}
			break;}
			case FlexAR::SINGLE_BARCODED:
			case FlexAR::SINGLE_BARCODED_WITHIN_READ:
			case FlexAR::SINGLE:{
			if(myRead->m_r1!=NULL){
				if(length(myRead->m_r1->getSequence())>=m_minLength){
					m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
					//m_outputMap[myRead->m_barcode_id].m_cnt_good_1+=1;
				}
				else m_outputMap[myRead->m_barcode_id].m_cnt_short_1+=1;
			}
			}

		}

		delete myRead;
		return NULL;

	}

	void writeLengthDist(){
		for(unsigned int i = 0;i < m_mapsize; ++i){
			m_outputMap[i].f1->writeLengthDist();
			if(m_outputMap[i].f2!=NULL)m_outputMap[i].f2->writeLengthDist();
		}
	}

	void printFilesSummary(){
		for(unsigned int i = 0;i < m_mapsize; ++i){
			std::cout << "File: " << m_outputMap[i].f1->getFileName() << std::endl;
			std::cout << "Nr. of reads dropped due to being shorter than minLength: " << m_outputMap[i].m_cnt_short_1 << std::endl;
			std::cout << "Nr. of reads written to the file: " << m_outputMap[i].f1->getNrGoodReads() << std::endl;

			if((m_runType==FlexAR::PAIRED_BARCODED)||(m_runType==FlexAR::PAIRED)){
				std::cout << "File: " << m_outputMap[i].f2->getFileName() << std::endl;
				std::cout << "Nr. of reads dropped due to being shorter than minLength: " << m_outputMap[i].m_cnt_short_2 << std::endl;
				std::cout << "Nr. of reads written to the file: " << m_outputMap[i].f2->getNrGoodReads() << std::endl;
			}

		}
	}


};

#endif /* MULTIPLEXINPUTFILTER_H_ */
