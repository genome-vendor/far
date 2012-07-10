/*
 * SequenceConverter.h
 *
 *  Created on: Jul 6, 2010
 *      Author: mat
 */

#ifndef SEQUENCECONVERTER_H_
#define SEQUENCECONVERTER_H_
template <class TString>
class SequenceConverter {
private:
	static SequenceConverter<TString>* instance;

	SequenceConverter(){
	};

public:
	static SequenceConverter<TString>* getInstance(){
		if(instance==NULL)instance = new SequenceConverter();
		return instance;
	}

	TString colorSpaceToBasepairSpace(TString csSequence){
		TString result = csSequence[0];
		TString substr = "XX";
		for(unsigned int i=1;i<length(csSequence);++i){
			substr[0] = result[i-1];
			substr[1] = csSequence[i];
			if(substr=="T0") append(result , "T");
			if(substr=="T1") append(result , "G");
			if(substr=="T2") append(result , "C");
			if(substr=="T3") append(result , "A");
			if(substr=="C0") append(result , "C");
			if(substr=="C1") append(result , "A");
			if(substr=="C2") append(result , "T");
			if(substr=="C3") append(result , "G");
			if(substr=="G0") append(result , "G");
			if(substr=="G1") append(result , "T");
			if(substr=="G2") append(result , "A");
			if(substr=="G3") append(result , "C");
			if(substr=="A0") append(result , "A");
			if(substr=="A1") append(result , "C");
			if(substr=="A2") append(result , "G");
			if(substr=="A3") append(result , "T");
		}
		return result;
	}

	TString basepairSpaceToColorSpace(TString csSequence){
		TString result = csSequence[0];
		TString substr = "XX";
		for(int i=1;i<length(csSequence);++i){
			substr[0] = result[i-1];
			substr[1] = csSequence[i];
			if(substr=="TT") append(result ,  "0");
			if(substr=="TG") append(result ,  "1");
			if(substr=="TC") append(result ,  "2");
			if(substr=="TA") append(result ,  "3");
			if(substr=="CC") append(result ,  "0");
			if(substr=="CA") append(result ,  "1");
			if(substr=="CT") append(result ,  "2");
			if(substr=="CG") append(result ,  "3");
			if(substr=="GG") append(result ,  "0");
			if(substr=="GT") append(result ,  "1");
			if(substr=="GA") append(result ,  "2");
			if(substr=="GC") append(result ,  "3");
			if(substr=="AA") append(result ,  "0");
			if(substr=="AC") append(result ,  "1");
			if(substr=="AG") append(result ,  "2");
			if(substr=="AT") append(result ,  "3");
		}
		return result;
	}

	TString getColorcodeFromT(TString substr){
		TString result;
		if(substr=="T") append(result ,  "0");
		if(substr=="G") append(result ,  "1");
		if(substr=="C") append(result ,  "2");
		if(substr=="A") append(result ,  "3");
		return result;
	}

	virtual ~SequenceConverter(){
	};
};

template <typename TString> SequenceConverter<TString>* SequenceConverter<TString>::instance = 0;

#endif /* SEQUENCECONVERTER_H_ */


