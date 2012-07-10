/*
 * NeedlemanWunschAlignmentAlgorithm.h
 *
 *  Created on: Jul 15, 2010
 *      Author: mat
 */

#ifndef ExactStringMatchingAlgorithm_H_
#define ExactStringMatchingAlgorithm_H_
#include <iostream>
template <typename TString>
class ExactStringMatchingAlgorithm {

private:
	int m_scoreMatch;
	int m_scoreMismatch;
	int m_scoreGapInsert;
	int m_scoreGapExtend;
	seqan::Score<int> *m_score;
public:
	ExactStringMatchingAlgorithm(int scoreMatch, int scoreMismatch, int scoreGapInsert, int scoreGapExtend){
		m_scoreMatch = scoreMatch;
		m_scoreMismatch = scoreMismatch;
		m_scoreGapInsert = scoreGapInsert;
		m_scoreGapExtend = scoreGapExtend;

		m_score = new seqan::Score<int>(scoreMatch,scoreMismatch,scoreGapInsert);
	};

	virtual ~ExactStringMatchingAlgorithm(){
		delete m_score;
	};

	void align(TString &adapter, TString &sequence, TString &squality, int &gaps, int &mismatches, int &startPos, int &endPos, int &startPosA, int &endPosA, int &startPosS, int &endPosS, int &alignmentScore, std::stringstream &alString){
		int la = length(adapter);
		int ls = length(sequence);
		int matches = 0;

		if(la!=ls){
			mismatches = ls;
			matches = 0;
		}
		else {
			mismatches = 0;

			alignmentScore = 0;
			for(register int i=0;i<la;++i){
				if(adapter[i]!=sequence[i])++mismatches;
				else ++matches;
			}

			startPosS = 0;
			endPosS = ls;
			startPosA = 0;
			endPosA = la;
			startPos = 0;
			endPos = ls;
		}
		alignmentScore = matches * m_scoreMatch - (mismatches * m_scoreMismatch);

	}

};

#endif /* NeedlemanWunschAlignmentAlgorithm_H_ */
