/*
 * NeedlemanWunschQualityAlignmentAlgorithm.h
 *
 *  Created on: Jul 15, 2010
 *      Author: mat
 */

#ifndef NeedlemanWunschQualityAlignmentAlgorithm_H_
#define NeedlemanWunschQualityAlignmentAlgorithm_H_
#include <iostream>
template <typename TString>
class NeedlemanWunschQualityAlignmentAlgorithm {

private:
	int m_scoreMatch;
	int m_scoreMismatch;
	int m_scoreGapInsert;
	int m_scoreGapExtend;
	seqan::Score<int> *m_score;
public:
	NeedlemanWunschQualityAlignmentAlgorithm(int scoreMatch, int scoreMismatch, int scoreGapInsert, int scoreGapExtend){
		m_scoreMatch = scoreMatch;
		m_scoreMismatch = scoreMismatch;
		m_scoreGapInsert = scoreGapInsert;
		m_scoreGapExtend = scoreGapExtend;

		m_score = new seqan::Score<int>(scoreMatch,scoreMismatch,scoreGapInsert);
	};

	virtual ~NeedlemanWunschQualityAlignmentAlgorithm(){
		delete m_score;
	};

	void align(TString &adapter, TString &sequence, TString &squality, int &gaps, int &mismatches, int &startPos, int &endPos, int &startPosA, int &endPosA, int &startPosS, int &endPosS, int &alignmentScore, std::stringstream &alString){
		using namespace seqan;
		typedef Align<TString, seqan::ArrayGaps> AlignType;
		AlignType align;
		resize(rows(align), 3);

		AlignConfig<true,true,true,true> ac;

		assignSource(row(align, 0), sequence);
		assignSource(row(align, 1), adapter);
		assignSource(row(align, 2), squality);
		alignmentScore = globalAlignment(align, *m_score, ac, NeedlemanWunschQuality() );


		alString << align;

		//int startPos=0,endPos=0;
		startPosS = beginPosition(row(align,0));
		endPosS = endPosition(row(align,0));
		startPosA = beginPosition(row(align,1));
		endPosA = endPosition(row(align,1));

		//calculate overlap
		if(endPosA > endPosS)
		{
			endPos = endPosS;
		}
		else
		{
			endPos = endPosA;
		}

		if(startPosA > startPosS)
		{
			startPos = startPosA;
		}
		else
		{
			startPos = startPosS;
		}

		//Compute number of mismatches and gaps
		typedef Gaps<TString,ArrayGaps> GapType;
		GapType row1,row2;

		row1 = row(align,0);
		row2 = row(align,1);

		for(register int i=startPos;i<endPos;++i)
		{
			if(isGap(row1,i)||isGap(row2,i)){
				++gaps;
			}
			else{
				if(getValue(row1,i) != getValue(row2,i)){
					++mismatches;
				}
			}
		}
	}

};

#endif /* NeedlemanWunschQualityAlignmentAlgorithm_H_ */
