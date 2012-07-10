/*
 * Enums.h
 *
 *  Created on: Jun 29, 2010
 *      Author: mat
 */

#ifndef ENUMS_H_
#define ENUMS_H_
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
	RIGHT_TAIL
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
}
#endif /* ENUMS_H_ */
