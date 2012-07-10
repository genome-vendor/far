/*
 * Enums.h
 *
 *  Created on: Jun 29, 2010
 *      Author: mat
 */

#ifndef ENUMS_H_
#define ENUMS_H_

/**These enums are used by almost every class. */
namespace FlexAR{
	const unsigned int MAX_READLENGTH=2048;
	enum LogLevel {
		NONE,
		ALL,
		TAB,
		CHANGED
	};

	enum TrimEnd {
		ANY,
		LEFT,
		RIGHT,
		LEFT_TAIL,
		RIGHT_TAIL,
		OFF
	};

	enum FileFormat {
		FASTA,
		FASTQ,
		CSFASTA,
		CSFASTQ
	};
	enum QualityType {
		SANGER,
		SOLEXA,
		ILLUMINA13,
		ILLUMINA15
	};
	enum RunType {
		PAIRED_BARCODED,
		SINGLE_BARCODED,
		PAIRED,
		SINGLE,
		SINGLE_BARCODED_WITHIN_READ

	};

	std::string toFormatString(FlexAR::FileFormat format){
		switch(format){
			case FlexAR::FASTA:return ".fasta";
			case FlexAR::FASTQ:return ".fastq";
			case FlexAR::CSFASTA:return ".csfasta";
			case FlexAR::CSFASTQ:return ".csfastq";
		}

		return ".unknown";

	}

}
#endif /* ENUMS_H_ */
