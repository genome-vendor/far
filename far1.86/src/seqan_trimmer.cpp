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
#include "SequenceOutputFilter.h"
#include "AdapterLoader.h"
#include "AlignmentFilter.h"
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

	time_t t_start,t_end;
	double totalTime;

	time(&t_start);

	//namespace po = boost::program_options;
	string input,output;

	// Declare the supported options.
	std::string program = "Flexible Adapter Remover, version 1.86 ";

	std::stringstream ss;
	ss << program << " - Allowed options";

	CommandLineParser cmdParser(program.c_str());

	addOption(cmdParser,CommandLineOption("v","version","return program version",OptionType::Label));
	addOption(cmdParser,CommandLineOption("s","source", "input file", OptionType::String | OptionType::Mandatory));
	addOption(cmdParser,CommandLineOption("t","target", "output file", OptionType::String | OptionType::Mandatory));
	addOption(cmdParser,CommandLineOption("f","format", "input file format - fastq (=fastq-illumina15), fastq-illumina13, fastq-sanger, fastq-solexa, fasta,csfastq,csfasta default is fasta, output will be in the same format", OptionType::String));
	addOption(cmdParser,CommandLineOption("a","adapters", "FASTA file with adaptersequences to be removed (will try subsequently, but only remove one). ", OptionType::String));
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
	addOption(cmdParser,CommandLineOption("of","omitted-reads-file", "File which contains reads, that were to short after adapter removal according to min-readlength. (default: input_file.omitted)",OptionType::String));
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

	std::string fileName, targetFileName, omittedFilename, trim_end, log_level;
	float cutOff=2;
	bool isColorSpace = false;
	CharString format = "fasta";
	FlexAR::FileFormat fformat;
	fformat = FlexAR::FASTA;
	int pre_cut_length = 0;
	int post_cut_length = 0;
	int phred_pre_cut_qual = 0;
	int allowedUncalledBases = 0;
	FlexAR::TrimEnd end;
	string alignMode = "", adapterFile;
	CharString adapterString="";
	int min_overlap_length = 10;
	int min_read_length = 18;
	FlexAR::LogLevel log = FlexAR::NONE;
	int nrThreads = 1;

	CharString source, dest, modifiedTag;

	int match;
	int mismatch;
	int gapOpen;
	//in all algorithms currently unused:
	int gapExt=0;

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
			qual = ::FlexAR::ILLUMINA13;
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
/*
		//solexa/sanger/illumina15/illumina13
		getOptionValueLong(cmdParser,"quality-scaling",phred_qual_str);
		if(phred_qual_str=="sanger"){
			qual = ::FlexAR::SANGER;
		}

		if(phred_qual_str=="solexa"){
			qual = ::FlexAR::SOLEXA;
		}

		if(phred_qual_str=="illumina13"){
			qual = ::FlexAR::ILLUMINA13;
		}
*/
	}



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
	bool overlapMode = false;
	getOptionValueLong(cmdParser,"adaptive-overlap",adaptive_overlap);
	if(adaptive_overlap=="yes"){
		std::cout << "Overlap will be treated as adaptive-overlap. " << std::endl;
		overlapMode = true;
	}
	else {
		std::cout << "Overlap will be treated as standard overlap. " << std::endl;
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

	getOptionValueLong(cmdParser,"omitted-reads-file",omittedFilename);
	if(omittedFilename=="")omittedFilename = fileName + ".omitted";

	getOptionValueLong(cmdParser,"modified-tag",modifiedTag);

	/////////////////////////////////////
	///	Finished command line parsing ///
	/////////////////////////////////////

	cout << "Writing omitted read ids to: " << omittedFilename << endl << endl;
	typedef std::pair<SequencingRead<CharString,CharString>*,unsigned int> TAdapter;
	tbb::concurrent_vector<TAdapter> adapters;

	if(length(adapterString)==0){
		// Start task scheduler
		tbb::task_scheduler_init init_serial( nrThreads );

		// Create the pipeline
		tbb::pipeline prepipeline;

		// Read adapters file
		SequenceInputFilter<CharString,CharString > adapter_filter( adapterFile, FlexAR::FASTA );
		adapter_filter.setNrAllowedUncalledBases(0);
		adapter_filter.setFixedPreTrim(0);
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

	/*
	 * The following code duplication is to ease compiler optimizations...
	 */

	// Create alignment filter
	if((alignMode == "needlemanQuality")||(alignMode == "needlemanquality"))
	{
		// Create the pipeline
		tbb::pipeline pipeline;

		// Create file-reading writing stage and add it to the pipeline
		SequenceInputFilter<CharString,CharString > input_filter( fileName, fformat );
		input_filter.setNrAllowedUncalledBases(allowedUncalledBases);
		input_filter.setFixedPreTrim(pre_cut_length);
		input_filter.SetMinReadLength(min_read_length);
		input_filter.setPrePhredTrim(phred_pre_cut_qual,qual);
		pipeline.add_filter( input_filter );

		std::cout << "Starting Trimming (algorithm: needlemanQuality)..." << std::endl;
		AlignmentFilter<CharString,CharString,NeedlemanWunschQualityAlignmentAlgorithm<CharString> > filter(&adapters,match,mismatch,gapOpen,gapExt, end, log);

		filter.setMinOverlap(min_overlap_length,overlapMode);
		filter.setCutOff(cutOff);
		if(modifiedTag!="NONE")filter.setModifiedTag(modifiedTag);


		pipeline.add_filter( filter );

		// Create file-writing stage and add it to the pipeline
		SequenceOutputFilter<CharString, CharString > output_filter( targetFileName,fformat );
		output_filter.SetMinReadLength(min_read_length);
		output_filter.SetOmittedFilename(omittedFilename);
		output_filter.setFixedPostTrim(post_cut_length);


		pipeline.add_filter( output_filter );

		// Run the pipeline
		tbb::tick_count t0 = tbb::tick_count::now();
		pipeline.run( nrThreads );
		tbb::tick_count t1 = tbb::tick_count::now();

		std::cout << "Done." << std::endl  << std::endl;

		ss.clear();
		ss.str("");
		std::cout << "Input file contained " << input_filter.getNrUncalledDiscarded() + input_filter.getNrProcessedReads() << " reads." << std::endl;
		std::cout << "Discarded " << input_filter.getNrUncalledDiscarded() << " reads due to containing uncalled bases." << std::endl;
		std::cout << "Discarded " << input_filter.getNrLowPhredDiscarded() << " reads due to having low quality." << std::endl;
		std::cout << "Used      " << input_filter.getNrProcessedReads() << " reads from input file." << std::endl;

		ss << std::endl << std::endl << "> Removed " << filter.getNrModifiedReads() << " adapters. " << std::endl;
		if(filter.getNrModifiedReads()>0){
		 ss << "Min-/Max-/Mean-/Median-overlap length: "<< filter.getMinOverlapLength()  << " / "<< filter.getMaxOverlapLength() << " / " << filter.getMeanOverlapLength() << " / " << filter.getMedianOverlapLength() << std::endl;
		}
		std::cout << ss.str();



		if(!isColorSpace)
		{
			std::cout << "Discarded " << output_filter.getNrOmittedReads() << " reads (shorter than " << min_read_length << " bp after adapter removal)." << std::endl;
		}
		else
		{
			std::cout << "Discarded " << output_filter.getNrOmittedReads() << " reads (shorter than " << min_read_length-1 << " bp after adapter removal)." << std::endl;
		}

		float nr = static_cast<float>(input_filter.getNrUncalledDiscarded() + input_filter.getNrProcessedReads());
		if(nr>0){
			std::cout << output_filter.getNrGoodReads() << " reads remaining ( " << static_cast<float>(output_filter.getNrGoodReads()) / nr * 100 << " % of input reads )" << std::endl << std::endl;
		}

		time(&t_end);
		totalTime = difftime(t_end,t_start);
		std::cout << "Calculation Time: " << div(static_cast<int>(totalTime),60).quot << " minutes "<< div(static_cast<int>(totalTime), 60).rem << " seconds. " << std::endl;

		CharString writeLengthdist;
		getOptionValueLong(cmdParser,"write-lengthdist",writeLengthdist);
		if(writeLengthdist=="yes")output_filter.writeLengthDist();
	}

	if(alignMode == "needleman")
	{
		// Create the pipeline
		tbb::pipeline pipeline;

		SequenceInputFilter<CharString,CharString > input_filter( fileName, fformat );
		input_filter.setNrAllowedUncalledBases(allowedUncalledBases);
		input_filter.setFixedPreTrim(pre_cut_length);
		input_filter.SetMinReadLength(min_read_length);
		input_filter.setPrePhredTrim(phred_pre_cut_qual,qual);
		pipeline.add_filter( input_filter );


		std::cout << "Starting Trimming (algorithm: needleman)..." << std::endl;
		AlignmentFilter<CharString,CharString,NeedlemanWunschAlignmentAlgorithm<CharString> > filter(&adapters,match,mismatch,gapOpen,gapExt, end, log);

		CharString tag;
		getOptionValueLong(cmdParser,"modified-tag",tag);
		if(tag!="NONE")filter.setModifiedTag(tag);

		filter.setMinOverlap(min_overlap_length,overlapMode);
		filter.setCutOff(cutOff);

		pipeline.add_filter( filter );
		SequenceOutputFilter<CharString, CharString > output_filter( targetFileName,fformat );
		output_filter.SetMinReadLength(min_read_length);
		output_filter.SetOmittedFilename(omittedFilename);
		output_filter.setFixedPostTrim(post_cut_length);


		pipeline.add_filter( output_filter );

		// Run the pipeline
		tbb::tick_count t0 = tbb::tick_count::now();
		pipeline.run( nrThreads );
		tbb::tick_count t1 = tbb::tick_count::now();

		std::cout << "Done." << std::endl  << std::endl;

		ss.clear();
		ss.str("");
		std::cout << "Input file contained " << input_filter.getNrUncalledDiscarded() + input_filter.getNrProcessedReads() << " reads." << std::endl;
		std::cout << "Discarded " << input_filter.getNrUncalledDiscarded() << " reads due to containing uncalled bases." << std::endl;
		std::cout << "Discarded " << input_filter.getNrLowPhredDiscarded() << " reads due to having low quality." << std::endl;
		std::cout << "Used      " << input_filter.getNrProcessedReads() << " reads from input file." << std::endl;

		ss << std::endl << std::endl << "> Removed " << filter.getNrModifiedReads() << " adapters. " << std::endl;
		if(filter.getNrModifiedReads()>0){
		 ss << "Min-/Max-/Mean-/Median-overlap length: "<< filter.getMinOverlapLength()  << " / "<< filter.getMaxOverlapLength() << " / " << filter.getMeanOverlapLength() << " / " << filter.getMedianOverlapLength() << std::endl;
		}
		std::cout << ss.str();


		if(!isColorSpace)
		{
			std::cout << "Discarded " << output_filter.getNrOmittedReads() << " reads (shorter than " << min_read_length << " bp after adapter removal)." << std::endl;
		}
		else
		{
			std::cout << "Discarded " << output_filter.getNrOmittedReads() << " reads (shorter than " << min_read_length-1 << " bp after adapter removal)." << std::endl;
		}

		float nr = static_cast<float>(input_filter.getNrUncalledDiscarded() + input_filter.getNrProcessedReads());
		if(nr>0){
			std::cout << output_filter.getNrGoodReads() << " reads remaining ( " << static_cast<float>(output_filter.getNrGoodReads()) / nr * 100 << " % of input reads )" << std::endl << std::endl;
		}

		time(&t_end);
		totalTime = difftime(t_end,t_start);
		std::cout << "Calculation Time: " << div(static_cast<int>(totalTime),60).quot << " minutes "<< div(static_cast<int>(totalTime), 60).rem << " seconds. " << std::endl;

		CharString writeLengthdist;
		getOptionValueLong(cmdParser,"write-lengthdist",writeLengthdist);
		if(writeLengthdist=="yes")output_filter.writeLengthDist();

		if(modifiedTag!="NONE")filter.setModifiedTag(modifiedTag);

	}

	//Now output statistics
	std::cout << "Adapter\tremoval_count" << std::endl;
	for(unsigned int i=0;i<adapters.size();++i){
		std::cout << adapters.at(i).first->getSequenceTag() << "\t" << adapters.at(i).second << std::endl;
	}

	return 0;
}
