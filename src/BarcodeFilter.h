/*
 * BarcodeFilter.h
 *
 *  Created on: May 31, 2011
 *      Author: mat
 */

#ifndef BARCODEFILTER_H_
#define BARCODEFILTER_H_

template <class TString, class TIDString, class TAlgorithm >
class BarcodeFilter {//: public tbb::filter{
private:
	FlexAR::LogLevel m_verbose;
	int m_scoreMatch;
	int m_scoreMismatch;
	int m_scoreGapInsert;
	int m_scoreGapExtend;
	int m_mis_thresh;
	int m_min_overlap;
	TAlgorithm *m_alg;
	FlexAR::TrimEnd m_trimEnd;

	typedef ::std::pair<SequencingRead<seqan::CharString,seqan::CharString>*,unsigned int> TAdapter;
	tbb::concurrent_vector<TAdapter> *m_adapters;



public:
	BarcodeFilter(tbb::concurrent_vector<TAdapter> *adapters,int thresh, int min_overlap, FlexAR::TrimEnd end,FlexAR::LogLevel log = FlexAR::NONE) {//: tbb::filter(parallel) {
		m_verbose = log;
		m_adapters = adapters;
		m_scoreMatch = 3;
		m_scoreMismatch = -3;
		m_scoreGapExtend = -1;
		m_scoreGapInsert = -5;
		m_trimEnd = end;
		m_mis_thresh = thresh;
		m_min_overlap = min_overlap;
		//if(m_verbose==FlexAR::TAB)std::cout << "read-ID\tadapter-start\tadapter-end\toverlap-length\tmismatches\tindels\tallowed-errors" << std::endl;
	};
	//void* operator()( void* item ){
	
	/*! This method will align all barcode sequences to the passed read and return the ID of the read that aligns best.*/
	int classify( void* item , FlexAR::TrimEnd bend){
			std::stringstream ss;

			using namespace seqan;
			using namespace std;
			if(item!=NULL)
			{
				MultiplexedRead<TString,TIDString> *myMultiRead = static_cast< MultiplexedRead<TString,TIDString>* >(item);
			
				SequencingRead<TString,TIDString> myRead;
				
				if(m_bend!=FlexAR::BNONE){
					//the barcode read has to be aligned with the barcode sequences
					 &myRead = *static_cast< SequencingRead<TString,TIDString>* >(myMultiRead->m_b);
				}
				else {
					//the first read has to be aligned with the barcode sequences
					&myRead = *static_cast< SequencingRead<TString,TIDString>* >(myMultiRead->m_r1);
				}

				int readLength=0,adapterLength=0,overlapLength=0;

				TString bpRead="";
				TString sequence="",squality = "";

				switch(myRead.getFormat()){
					case FlexAR::CSFASTQ:
						sequence = suffix<TString>(myRead.getSequence(),2);
						squality = suffix<TString>(myRead.getQuality(),1);
						break;
					case FlexAR::FASTQ:
						sequence = myRead.getSequence();
						squality = myRead.getQuality();
						break;
					case FlexAR::CSFASTA:
						sequence = suffix<TString>(myRead.getSequence(),2);
						squality = "";
						break;
					case FlexAR::FASTA:
						sequence = myRead.getSequence();
						squality = "";
						break;
				}

				//common calculations
				readLength = length(sequence);

				TString current_adapter;
				TAlgorithm alg(m_scoreMatch,m_scoreMismatch,m_scoreGapInsert,m_scoreGapExtend);

				//Align all barcodes and remember the alignment scores
				//best one will be used
				int scoreIndex=-1,max_mismatches,max_overlap;
				float scoreMax=-10000;

				//Iterate over all barcodes & align
				for(unsigned int i=0;i<m_adapters->size();++i)
				{
					current_adapter = m_adapters->at(i).first->getSequence();
					adapterLength = length(current_adapter);
					//In tail mode trim read before alignment
					int startPos=0,endPos=0,startPosA=0,endPosA=0,startPosS=0,endPosS=0,alignmentScore=0,mismatches=0,gaps=0;


					std::stringstream alString;

					alg.align(current_adapter,sequence, squality, gaps, mismatches, startPos, endPos, startPosA, endPosA, startPosS, endPosS, alignmentScore, alString);

					//calculate overlap
					overlapLength = endPos - startPos;

					if(m_verbose==FlexAR::ALL){
						ss << endl;
						ss << "read-tag     : " << myRead.getSequenceTag() << endl;
						ss << "barcode-tag	: " << m_adapters->at(i).first->getSequenceTag()<< endl;
						ss << "read        : " << sequence << endl;
						ss << "quality     : " << squality << endl;

						ss << "beginR       : " << startPosS << endl;
						ss << "endR         : " << endPosS << endl;

						ss << "beginA       : " << startPosA << endl;
						ss << "endA         : " << endPosA << endl;

						ss << alString.str() << endl;
						ss << "overlap(bp) : " << overlapLength << endl;
						ss << "score       : " << alignmentScore << endl;

						ss << "indels        : " << gaps << endl;
						ss << "mismatches    : " << mismatches << endl;
						ss << "allowed errors: " << m_mis_thresh << endl;
						//ss << "cutOff        : " << m_cutOff << endl;
						ss << endl;
					}


					//remember min/max alignment score
					if(static_cast<float>(alignmentScore)>scoreMax){
						scoreMax = static_cast<float>(alignmentScore);
						scoreIndex = i;
						max_mismatches = mismatches;
						max_overlap = overlapLength;
					}


				}

				//now return the fragment of the read that remains (if barcode was within read)...
				if(m_bend!=FlexAR::BNONE){
					//cut the read depending on b_end
				
				}
				

				//Bundeled output for multithreading
				if(m_verbose!=FlexAR::NONE)cout << ss.str();

				//Recheck if overlap was sufficient
				if((max_mismatches>m_mis_thresh)||(max_overlap<m_min_overlap)){
					if(m_verbose==FlexAR::ALL)std::cout << "No barcode identified." << std::endl;
					return 0;
				}
				else {
					if(m_verbose==FlexAR::ALL)std::cout << "barcode " << m_adapters->at(scoreIndex).first->getSequenceTag() << " identified." << std::endl;
					return scoreIndex + 1;
				}
				//return &myRead;
			}
			else return 0;
	}
	virtual ~BarcodeFilter(){};
};

#endif /* BARCODEFILTER_H_ */
