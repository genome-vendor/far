/*
 * MultiplexedRead.h
 *
 *  Created on: May 31, 2011
 *      Author: mat
 */

#ifndef MULTIPLEXEDREAD_H_
#define MULTIPLEXEDREAD_H_

/** This class represents either a single read or a paired read. In both cases a barcode-read can be present also (that is used to align a given set of barcodes to this read and identify the corresponding experiment).*/
template <class TString,class TIDString>
class MultiplexedRead {
	public:
		typedef SequencingRead<TString, TIDString> TSequencingRead;
		TSequencingRead *m_r1;
		TSequencingRead *m_r2;
		TSequencingRead *m_b;
		int m_barcode_id;

		MultiplexedRead(TSequencingRead *r1,TSequencingRead *r2,TSequencingRead *b){
			m_r1 = r1;
			m_r2 = r2;
			m_b = b;
			m_barcode_id = 0;
		};

		MultiplexedRead(TSequencingRead *r1,TSequencingRead *r2){
			m_r1 = r1;
			m_r2 = r2;
			m_b = NULL;
			m_barcode_id = 0;
		};

		MultiplexedRead(TSequencingRead *r1){
			m_r1 = r1;
			m_r2 = NULL;
			m_b = NULL;
			m_barcode_id = 0;
		};
		virtual ~MultiplexedRead(){
			delete m_r1;
			delete m_r2;
			delete m_b;
			m_r1 = m_r2 = m_b = NULL;
			m_barcode_id = 0;
		};
};

#endif /* MULTIPLEXEDREAD_H_ */
