//============================================================================
// Name        : seqan_trimmer.cpp
// Author      : Matthias Dodt
// Copyright   : GPL V3
// Description : Flexible adapter remover
//============================================================================

#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <string>

#include <tbb/pipeline.h>
#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>

#include "SequenceInputFilter.h"
#include "MultiplexedInputFilter.h"
#include "MultiplexedOutputFilter.h"
#include "AdapterLoader.h"
#include "MultiplexedAlignmentFilter.h"
#include "NeedlemanWunschAlignmentAlgorithm.h"
#include "NeedlemanWunschQualityAlignmentAlgorithm.h"
#include "Enums.h"



using namespace std;
using namespace seqan;

void showProgress(int percent){
	std::cout << "\r";
	for(int i=0;i<percent;++i)std::cout << "=";

	std::cout << " " << percent << "%" << flush;
}

unsigned long int getFileLength(std::fstream &file)
{
	unsigned long int pos=0,length=0;

	if(file.is_open())
	{
		  // get length of file:
		  pos = file.tellg();
		  file.seekg (0, ios::end);
		  length = file.tellg();
		  file.seekg (pos);
		  return length;
	}
	else return 0;

}


int main(int argc, const char* argv[] ) {
	typedef std::pair<SequencingRead<CharString,CharString>*,unsigned int> TAdapter;
	time_t t_start,t_end;
	double totalTime;

	time(&t_start);

	//namespace po = boost::program_options;
	string input,output;

	// Declare the supported options.
	std::string program = "Flexible Adapter Remover, version 2.17 ";

	std::stringstream ss;
	ss << program << " - Allowed options";

	CommandLineParser cmdParser(program.c_str());

	addOption(cmdParser,CommandLineOption("v","version","return program version",OptionType::Label));
	addOption(cmdParser,CommandLineOption("s","source", "input file", OptionType::String | OptionType::Mandatory));
	addOption(cmdParser,CommandLineOption("s2","source2", "second input file (if paired run)", OptionType::String));
	addOption(cmdParser,CommandLineOption("t","target", "output file", OptionType::String | OptionType::Mandatory));
	addOption(cmdParser,CommandLineOption("f","format", "input file format - fastq (=fastq-illumina15), fastq-illumina13, fastq-sanger, fastq-solexa, fasta,csfastq,csfasta default is fasta, output will be in the same format", OptionType::String));
	addOption(cmdParser,CommandLineOption("a","adapters", "FASTA file with adaptersequences to be removed (will try subsequently, but only remove one). ", OptionType::String));
	addOption(cmdParser,CommandLineOption("b","barcodes", "FASTA file with barcode sequences. ",OptionType::String));
	addOption(cmdParser,CommandLineOption("br","barcode-reads", "FASTA/Q file with barcode sequences. Leave empty if barcodes are within the reads. ",OptionType::String));
	addOption(cmdParser,CommandLineOption("bm","barcode-mismatch", "Allowed mismatches for a barcode ",(int)OptionType::Int,1));
	addOption(cmdParser,CommandLineOption("bmo","barcode-min-overlap", "Minimum required overlap between known barcode and barcode sequence (default: barcode length - 1 )",(int)OptionType::Int));
	addOption(cmdParser,CommandLineOption("be","barcode-trim-end", "If you specify this option the barcode is assumed to be within the read. Trim options are  to --trim-end: left/right/any ",OptionType::String));
	addOption(cmdParser,CommandLineOption("do","demultiplex-only", "yes/no (if set to yes only barcode classification will be done.)", OptionType::String,"no"));
	addOption(cmdParser,CommandLineOption("as","adapter", "String to be removed (alternative to --adapters option). ", OptionType::String));
	addOption(cmdParser,CommandLineOption("alg","algorithm", "needlemanQuality (quality based global alignment via needleman wunsch quality alignment) or needleman (complete global overlap alignment via needleman wunsch)",(int)OptionType::String ,"needleman"));
	addOption(cmdParser,CommandLineOption("pre","fixed-pre-trim", "if set, trims the reads to the number of specified bases before alignment",(int)OptionType::Int,-1));
	addOption(cmdParser,CommandLineOption("ao","adaptive-overlap", "If set to 'yes' the parameter '--min-overlap' will be treated as adaptive measure (see docs, default: no)",(int)OptionType::String));
	addOption(cmdParser,CommandLineOption("post","fixed-post-trim", "if set, trims the reads to the number of specified bases after alignment", (int)OptionType::Int,-1));
	addOption(cmdParser,CommandLineOption("phredpre","phred-pre-trim", "if set, trims the reads to the length were phred/sanger (scaling depends on --format setting) quality will be bigger than this value (searching from 3' end of reads)", (int)OptionType::Int,-1));
	addOption(cmdParser,CommandLineOption("end", "trim-end", "Decides which type of adapterremoval is performed. Modes: \n right - adapter is removed if adapter aligns >= read start position \n left - adapter is removed if adapter aligns <= read end position \n any - adapter is removed if adapter aligns, longer part of read remains \n left_tail - adapter is searched only in the first n bases of the read (n = adapterlength) \n right_tail - adapter is searched only in the last n bases of the read (n = adapterlength)",(int)OptionType::String,"right"));
	addOption(cmdParser,CommandLineOption("c","cut-off", "Max. nr of allowed mismatches + indels for alignment per 10 bases ",(int)OptionType::Double,2.0));
	addOption(cmdParser,CommandLineOption("o","min-overlap","minimum required overlap of adapter and sequence in basepairs ",(int)OptionType::Int,10));
	addOption(cmdParser,CommandLineOption("ml","min-readlength", "minimum readlength in basepairs after adapter removal - read will be discarded otherwise ",(int)OptionType::Int,18));
	addOption(cmdParser,CommandLineOption("of","omitted-reads-file", "Prefix (name) of file which contains reads, that were to short after adapter removal according to min-readlength. (default: input_file.omitted)",OptionType::String));
	addOption(cmdParser,CommandLineOption("u","max-uncalled", "nr of allowed uncalled bases in a read (might be a '.' or an 'N'). ",(int)OptionType::Int,0));

	addOption(cmdParser,CommandLineOption("sm","match", "match score ",(int)OptionType::String,"3"));
	addOption(cmdParser,CommandLineOption("sf","mismatch", "mismatch penalty ",(int)OptionType::String,"-3"));
	addOption(cmdParser,CommandLineOption("sg","gap-cost", "gap penalty ",(int)OptionType::String,"-20"));
	addOption(cmdParser,CommandLineOption("th", "nr-threads", "Number of threads to use ",(int)OptionType::Int,"1"));
	addOption(cmdParser,CommandLineOption("l","log-level","Report detailed Alignment information: \n ALL: Full alignment for each read \n CHANGED: Alignments of reads where adapters have been removed \n TAB: Same as changed, but TAB-wise output. \n This will slow down the processing! ",(int)OptionType::String,"NONE"));
	addOption(cmdParser,CommandLineOption("n","modified-tag", "All reads which are modified by FAR will have a colon followed by this tag in their readname.",(int)OptionType::String,"NONE"));
	addOption(cmdParser,CommandLineOption("d","write-lengthdist", "Writes a length distribution to the specified filename. By default it is: <targetfilename>.lengthdist (Options: yes/no)",(int)OptionType::String,"yes"));

	if(!parse(cmdParser,argc,argv,std::cerr)){
		cerr << "Error parsing command line. Check the mandatory parameters." << endl;
		help(cmdParser,std::cerr);
		return -1;
	}

	std::string fileName,fileName2, targetFileName, omittedFilename, trim_end, b_trim_end,log_level;
	float cutOff=2;
	bool isColorSpace = false;
	CharString format = "fasta";

	FlexAR::FileFormat fformat;
	fformat = FlexAR::FASTA;

	FlexAR::RunType runType = FlexAR::SINGLE;

	int pre_cut_length = 0;
	int post_cut_length = 0;
	int phred_pre_cut_qual = 0;
	int allowedUncalledBases = 0;
	FlexAR::TrimEnd end = FlexAR::RIGHT;
	///default looking for barcodes in reads is disabled (assumed to be in separate file)
	FlexAR::TrimEnd b_end = FlexAR::OFF;
	string alignMode = "", adapterFile;
	CharString adapterString="";
	int min_overlap_length = 10;
	int min_read_length = 18;
	FlexAR::LogLevel log = FlexAR::NONE;
	int nrThreads = 1;
	bool adaptiveOverlap = false;

	CharString source, dest, modifiedTag;

	int match=3;
	int mismatch=-3;
	int gapOpen=-20;
	//in all algorithms currently unused:
	int gapExt=0;
	tbb::concurrent_vector<TAdapter> adapters;

	if(getOptionValueLong(cmdParser,"version",fileName)){
		std::cout << program << std::endl;
		//exit(0);
		return -1;
	}

	if(getOptionValueLong(cmdParser,"source",fileName)){
		std::cout << "source file was set to " << fileName << ".\n";
	}
	else {
		cerr << "No input file specified!\n"<<endl;
		return -1;
	}

	getOptionValueLong(cmdParser,"source2",fileName2);
	if(length(fileName2)!=0){
		std::cout << "source file 2 was set to " << fileName2 << ".\n";
		runType = FlexAR::PAIRED;
	}
	else {
		std::cout << "No 2nd input file specified! (Run is not paired)\n"<<endl;
	}

	if(getOptionValueLong(cmdParser,"target",targetFileName)){
		cout << "target file was set to " << targetFileName << ".\n";
	}
	else {
		cerr << "No target file specified!\n" << endl;
		return -1;
	}

	::FlexAR::QualityType qual = ::FlexAR::ILLUMINA15;
	if(getOptionValueLong(cmdParser,"format",format)){
		if (format=="fasta")
		{
			cout << "File format was set to FASTA" << endl;
		}
		if ((format=="fastq")||(format=="fastq-illumina15")){
			fformat = FlexAR::FASTQ;
			cout << "File format was set to FASTQ-illumina15" << endl;
			qual = ::FlexAR::ILLUMINA15;
		}
		if (format=="fastq-illumina13"){
			fformat = FlexAR::FASTQ;
			qual = ::FlexAR::ILLUMINA13;
			cout << "File format was set to FASTQ-illumina13" << endl;
		}

		if (format=="fastq-sanger"){
			fformat = FlexAR::FASTQ;
			qual = ::FlexAR::SANGER;
			cout << "File format was set to FASTQ-sanger" << endl;
		}

		if (format=="fastq-solexa"){
			fformat = FlexAR::FASTQ;
			cout << "File format was set to FASTQ-solexa" << endl;
		}

		if (format=="csfasta")
		{
			fformat = FlexAR::CSFASTA;
			cout << "File format was set to colorspace FASTA" << endl;
			isColorSpace = true;
		}
		if (format=="csfastq"){
			fformat = FlexAR::CSFASTQ;
			isColorSpace = true;
			qual = ::FlexAR::SANGER;
			cout << "File format was set to colorspace FASTQ (assumed sanger quality scaling)" << endl;
		}

	}
	else {
		cerr << "Please specify an input-file format! "<< endl << endl;
		//help(cmdParser,std::cerr);
		return -1;
	}

	getOptionValueLong(cmdParser,"fixed-pre-trim",pre_cut_length);
	if(pre_cut_length > 0){
		std::cout << "Reads will be trimmed to " << pre_cut_length << " bases before alignment." << std::endl;
	}

	getOptionValueLong(cmdParser,"fixed-post-trim",post_cut_length);
	if(post_cut_length > 0){
		std::cout << "Reads will be trimmed to " << post_cut_length << " bases after alignment." << std::endl;
	}

	CharString phred_qual_str;
	getOptionValueLong(cmdParser,"phred-pre-trim",phred_pre_cut_qual);
	if(phred_pre_cut_qual > 0){
		std::cout << "Reads will be trimmed to a quality of " << phred_pre_cut_qual << " (from 3' end) before alignment." << std::endl;
	}

	bool demultiplexOnly = false;
	CharString demultiplex;
	getOptionValueLong(cmdParser,"demultiplex-only",demultiplex);
	if(demultiplex=="yes"){
		std::cout << "Will do demultiplexing only... " << std::endl;
		demultiplexOnly = true;
	}
	else {
		std::cout << "Will do demultiplexing and adapter removal... " << std::endl;
	
		getOptionValueLong(cmdParser,"adapters",adapterFile);
		if(length(adapterFile)==0)
		{
			getOptionValueLong(cmdParser,"adapter",adapterString);
			if(length(adapterString)==0)
			{
				cerr << "Please set an adapters.fa file or a adapter string sequence to be removed!" << endl << endl;
				return -1;
			}
			else {
				std::cout << "Adapter sequence: " << adapterString << std::endl<< std::endl;
			}

		}
		else {
			std::cout << "Adapter file: " << adapterFile << std::endl<< std::endl;
		}

		if(getOptionValueLong(cmdParser,"cut-off",cutOff)){
		std::cout << "Nr. of allowed indels + mismatches per 10 bases: " << cutOff << std::endl;
		}

		getOptionValueLong(cmdParser,"max-uncalled",allowedUncalledBases);
		cout << std::endl<<"Allowed number of uncalled bases was set to " << allowedUncalledBases << endl;

		CharString match_cost;
		getOptionValueLong(cmdParser,"match",match_cost);
		match = atoi(toCString(match_cost));

		CharString mismatch_cost;
		getOptionValueLong(cmdParser,"mismatch",mismatch_cost);
		mismatch = atoi(toCString(mismatch_cost));

		CharString gap_cost;
		getOptionValueLong(cmdParser,"gap-cost",gap_cost);
		gapOpen = atoi(toCString(gap_cost));

		std::cout << "Chosen scoring scheme: match = " << match << " ,mismatch = " << mismatch << " ,gap opening = " << gapOpen << std::endl;

		if(getOptionValueLong(cmdParser,"trim-end",trim_end)){
			if(trim_end=="left"){
				end = FlexAR::LEFT;
				std::cout << "Trimming from left end. " << std::endl;
			} else{
			if(trim_end=="right"){
				end = FlexAR::RIGHT;
				std::cout << "Trimming from right end. " << std::endl;
			} else {
			if(trim_end=="any"){
				end = FlexAR::ANY;
				std::cout << "Trimming: Global alignment, longest part of read remains. " << std::endl;
			} else {
			if(trim_end=="left_tail"){
				end = FlexAR::LEFT_TAIL;
				std::cout << "Trimming left tail only. " << std::endl;
			} else {
			if(trim_end=="right_tail"){
				end = FlexAR::RIGHT_TAIL;
				std::cout << "Trimming right tail only. " << std::endl;
			} else {
				std::cout << "trim-end specified is not a valid trim-end!" << std::endl;
				return -1;
			}
			}
			}
			}
			}

		}
		else {
			std::cerr << "Please specify a --trim-end ! (left/right/both)" << std::endl;
			return -1;
		}

		if (getOptionValueLong(cmdParser,"algorithm",alignMode)) {
			if((alignMode!="needlemanQuality")&&(alignMode!="needlemanquality")&&(alignMode!="needleman"))
			{
				cerr << "Align mode not specified or not supported!\n"<<endl;
				return -1;
			}
			else {
				if((alignMode=="quality")&&((fformat!=FlexAR::CSFASTQ)||(format!=FlexAR::FASTQ))){
					cerr << "Quality mode only works with quality files!" << std::endl;
					return -1;
				}
			}
		}

		getOptionValueLong(cmdParser,"min-overlap",min_overlap_length);
		std::cout << "Minimum required overlap: " << min_overlap_length << std::endl;

		CharString adaptive_overlap = "no";
		getOptionValueLong(cmdParser,"adaptive-overlap",adaptive_overlap);
		if(adaptive_overlap=="yes"){
			std::cout << "Overlap will be treated as adaptive-overlap. " << std::endl;
			adaptiveOverlap = true;
		}
		else {
			std::cout << "Overlap will be treated as standard overlap. " << std::endl;
		}

		if(length(adapterString)==0){
			// Start task scheduler
			tbb::task_scheduler_init init_serial( nrThreads );

			// Create the pipeline
			tbb::pipeline prepipeline;

			// Read adapters file
			SequenceInputFilter<CharString,CharString > adapter_filter( adapterFile, FlexAR::FASTA );
			//adapter_filter.setNrAllowedUncalledBases(0);
			//adapter_filter.setFixedPreTrim(0);
			prepipeline.add_filter( adapter_filter );

			AdapterLoader<CharString,CharString> adapterLoader;
			prepipeline.add_filter( adapterLoader );

			prepipeline.run(1);

			adapters = adapterLoader.getAdapters();
		}
		else {
			SequencingRead<CharString,CharString> *myRead = new SequencingRead<CharString,CharString>(adapterString,"adapter1",FlexAR::FASTA);
			TAdapter adap;
			adap.first = myRead;
			adap.second = 0;
			adapters.push_back(adap);
		}

		cout << "Adapter sequences are:" << endl;
		for (unsigned int i=0;i<adapters.size();++i){
			std::cout << adapters.at(i).first->getSequenceTag() << endl << adapters.at(i).first->getSequence() << endl;
		}
	}

	std::string barcodeFile,barcodeReadsFile,barcodeTrimEnd;
	//now do barcoding
	tbb::concurrent_vector<TAdapter> barcodes;
	int barcode_mismatch=1,barcode_min_overlap=6;
	getOptionValueLong(cmdParser,"barcodes",barcodeFile);
	if(length(barcodeFile)==0)
	{
		std::cout << "No barcodes file specified... " << std::endl<< std::endl;
	}
	else {
		///If previously two files were specified treat this as paired_barcoded, single_barcoded otherwise
		if(runType == FlexAR::PAIRED)runType = FlexAR::PAIRED_BARCODED;
		if(runType == FlexAR::SINGLE)runType = FlexAR::SINGLE_BARCODED;

		getOptionValueLong(cmdParser,"barcode-mismatch",barcode_mismatch);
		if(barcode_mismatch > 0){
			std::cout << "Allowed mismatches for barcodes: " << barcode_mismatch << std::endl;
		}

		std::cout << "Barcode file: " << barcodeFile << std::endl<< std::endl;

		// Start task scheduler
		tbb::task_scheduler_init init_serial( nrThreads );

		// Create the pipeline
		tbb::pipeline bpipeline;

		// Read adapters file
		SequenceInputFilter<CharString,CharString > adapter_filter( barcodeFile, FlexAR::FASTA );
		//adapter_filter.setNrAllowedUncalledBases(0);
		//adapter_filter.setFixedPreTrim(0);
		bpipeline.add_filter( adapter_filter );

		AdapterLoader<CharString,CharString> adapterLoader;
		bpipeline.add_filter( adapterLoader );

		bpipeline.run(1);
		barcodes = adapterLoader.getAdapters();

		cout << "Barcode sequences are:" << endl;
		for (unsigned int i=0;i<barcodes.size();++i){
			std::cout << barcodes.at(i).first->getSequenceTag() << endl << barcodes.at(i).first->getSequence() << endl;
		}

		getOptionValueLong(cmdParser,"barcode-trim-end",barcodeTrimEnd);
		if(length(barcodeTrimEnd)>0){
			runType = FlexAR::SINGLE_BARCODED_WITHIN_READ;
			std::cout << std::endl << "Barcodes are assumed to be within a read." << std::endl << std::endl;
			if(barcodeTrimEnd=="left"){
				b_end = FlexAR::LEFT;
				std::cout << "Trimming from left end. " << std::endl;
			} else{
			if(barcodeTrimEnd=="right"){
				b_end = FlexAR::RIGHT;
				std::cout << "Trimming from right end. " << std::endl;
			} else {
			if(barcodeTrimEnd=="any"){
				b_end = FlexAR::ANY;
				std::cout << "Trimming: Global alignment, longest part of read remains. " << std::endl;
			}}}
		}
		else {
			
			getOptionValueLong(cmdParser,"barcode-reads",barcodeReadsFile);
			if(length(barcodeReadsFile)==0)
			{
				std::cout << "No barcode reads file specified. " << std::endl<< std::endl;

				exit(0);
			}
			std::cout << std::endl << "Barcodes for each read are assumed to be in a separate file (" << barcodeReadsFile << ")" << std::endl << std::endl;
		}

		// Create the pipeline
		if(barcodes.size()==0){
			std::cerr << "No barcodes found in file. exiting." << std::endl;
			exit(0);
		}

		barcode_min_overlap = length(barcodes.at(0).first->getSequence()) - 1;
		getOptionValueLong(cmdParser,"barcode-min-overlap",barcode_min_overlap);
		if(barcode_min_overlap > 0){
			std::cout << "Required overlap for barcodes: " << barcode_min_overlap << std::endl;
		}
	}

	getOptionValueLong(cmdParser,"min-readlength",min_read_length);
	std::cout << "Minimum required readlength: " << min_read_length << std::endl;

	//in case of cs space
	if(isColorSpace)
	{
		++min_read_length;
	}

	if (getOptionValueLong(cmdParser,"log-level",log_level)){
		if(log_level=="ALL"){
			cout << "WARNING! Producing ALL verbose output - this will slow down the processing!" << endl;
			log = FlexAR::ALL;
		}
		if(log_level=="TAB"){
			cout << "WARNING! Producing TAB verbose output - this will slow down the processing!" << endl;
			log = FlexAR::TAB;
		}
		if(log_level=="CHANGED"){
			cout << "WARNING! Producing CHANGED verbose output - this will slow down the processing!" << endl;
			log = FlexAR::CHANGED;
		}
	}
	else
	{
		cout << "--log-level variable not set, using non verbose output." << endl;
	}

	getOptionValueLong(cmdParser,"nr-threads",nrThreads);
	cout << "Using " << nrThreads << " threads." << endl << endl;

	tbb::task_scheduler_init init_serial(nrThreads);

	getOptionValueLong(cmdParser,"omitted-reads-file",omittedFilename);
	if(omittedFilename=="")omittedFilename = fileName + ".omitted";

	getOptionValueLong(cmdParser,"modified-tag",modifiedTag);

	/////////////////////////////////////
	///	Finished command line parsing ///
	/////////////////////////////////////

	cout << "Common filename prefix for omitted read: " << omittedFilename << endl << endl;
	
	// Create the pipeline
	tbb::pipeline pipeline;

	// Create file-reading writing stage and add it to the pipeline
	MultiplexedInputFilter<UnicodeString,UnicodeString > input_filter( fileName, fformat ,allowedUncalledBases, pre_cut_length, min_read_length, phred_pre_cut_qual, qual);
	if((runType==FlexAR::PAIRED_BARCODED)||(runType==FlexAR::SINGLE_BARCODED))input_filter.setBarcodeReadsFile(barcodeReadsFile);
	if((runType==FlexAR::PAIRED)||(runType==FlexAR::PAIRED_BARCODED))input_filter.setPairedFile(fileName2);

	pipeline.add_filter( input_filter );

	std::cout << "Starting processing (algorithm: needleman-wunsch)..." << std::endl;
	if(log==FlexAR::TAB)std::cout << "read-ID\tadapter-ID\tadapter-start\tadapter-end\toverlap-length\tmismatches\tindels\tallowed-errors" << std::endl;
	MultiplexedAlignmentFilter<UnicodeString,UnicodeString,NeedlemanWunschAlignmentAlgorithm<UnicodeString> > filter(&adapters,&barcodes,barcode_mismatch,barcode_min_overlap,b_end,match,mismatch,gapOpen,gapExt, end, log);
	filter.setMinOverlap(min_overlap_length,adaptiveOverlap);
	filter.setCutOff(cutOff);
	if(demultiplexOnly)filter.setDemultiplexOnly();

	pipeline.add_filter( filter );

	// Create file-writing stage and add it to the pipeline
	MultiplexedOutputFilter<UnicodeString, UnicodeString > output_filter( targetFileName,&barcodes,fformat ,min_read_length,omittedFilename,post_cut_length,runType);
	pipeline.add_filter( output_filter );

	// Run the pipeline
	tbb::tick_count t0 = tbb::tick_count::now();
	pipeline.run( nrThreads );
	tbb::tick_count t1 = tbb::tick_count::now();
	std::cout << "Done." << std::endl  << std::endl;
	time(&t_end);
	totalTime = difftime(t_end,t_start);
	std::cout << "Calculation Time: " << div(static_cast<int>(totalTime),60).quot << " minutes "<< div(static_cast<int>(totalTime), 60).rem << " seconds. " << std::endl << std::endl;

	std::cout << "Step 1 - filtering input files: " << std::endl;
	std::cout << "=============================== " << std::endl;
	std::cout << "Processed  " << input_filter.getNrProcessedReads() << " reads from input file." << std::endl<< std::endl;
	std::cout << "Discarded in total " << input_filter.getNrUncalledReads() << " reads due to containing uncalled bases." << std::endl;
	std::cout << "Discarded in total " << input_filter.getNrLowPhredDiscarded() << " reads due to having low quality." << std::endl;
	
	float nr = static_cast<float>(input_filter.getNrUncalledReads() + input_filter.getNrProcessedReads());
	if(nr>0){
		std::cout << output_filter.getNrGoodReads() << " reads remaining ( " << static_cast<float>(output_filter.getNrGoodReads()) / nr * 100 << " % of input reads )" << std::endl << std::endl;
	}

	std::cout << "Statistics on each output file:" << std::endl;
	std::cout << "===============================" << std::endl;
	output_filter.printFilesSummary();

	ss.clear();
	ss.str("");

	CharString writeLengthdist;
	getOptionValueLong(cmdParser,"write-lengthdist",writeLengthdist);
	if(writeLengthdist=="yes"){
		std::cout << "Writing length distributions of reads (for each file) " << std::endl;
		output_filter.writeLengthDist();
	}

	//Now output statistics
	std::cout << std::endl << "Statistics on adapter removal (input files):" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Adapter\tremoval_count" << std::endl;
	for(unsigned int i=0;i<adapters.size();++i){
		std::cout << adapters.at(i).first->getSequenceTag() << "\t" << adapters.at(i).second << std::endl;
	}

	//this will output statistics on min/mean/max overlap between adapter & seequence
	filter.printAlignmentSummary();

	return 0;
}
