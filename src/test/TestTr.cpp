/*
 * TestTr.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: Hani Zakaria Girgis, PhD and Joseph Valencia
 */

#include "../nonltr/ChromosomeOneDigit.h"
#include "../nonltr/Chromosome.h"
#include "../nonltr/ChromosomeRandom.h"
#include "../nonltr/ChromListMaker.h" 

#include "../tr/ScorerTr.h"
#include "../tr/FilterTr.h"
#include "../tr/Tr.h"
#include "../tr/ForwardTr.h"
#include "../tr/BackwardTr.h"
#include "../tr/TrCollector.h"
#include "../utility/Util.h"
#include "../utility/LCS.h"
#include "../utility/LCSLen.h"
#include "../utility/Location.h"
#include "../utility/ILocation.h"
#include "../utility/LCSubStr.h"
#include "../utility/TSD.h"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <experimental/filesystem>

using namespace std;
using namespace nonltr;
using namespace utility;
using namespace tr;


namespace fs = std::experimental::filesystem;

int getOptionalArg(vector<string> * args, string option){
    string answer = "-1";
	auto it = std::find(args->begin(),args->end(),option);

	if(it!=args->end() && ++it !=args->end()){
		answer = * it;
	}
	return stoi(answer);
}

int main(int argc, char * argv[]) {
	int minDist = 400;
	int maxDist = 22000;
	int minLenLTR = 100;
	int maxLenLTR = 6000;
	int identity = 85;
	int k = 13;
	int minPlateauLength = 10;
	int gapTol = 200;
	int threads = 1;
	string chromDir ="";
	string seqFile = "";
	string outputDir = "";
	string outputPrefix = "";
	bool seqLevel = false;
	bool printRawScores = false;
	bool printCleanScores = false;
	bool displayHelp = false;
	bool bedOutput = false;
	bool displayNested = false;


	string helpMessage = 
	"| -arg     | Description | Default |\n"
	"| ---------------- | ----------- | ------- |\n"
	"| -seqFile | FASTA file to be scanned. May contain multiple sequences. | required unless -fasta or -chromDir is used |\n"
	"| -fasta | Directory containing files to be scanned. Files must have .fa extension. | required unless -seqFile is used |\n"
	"| -chromDir | Alias for -fasta. | required unless -seqFile is used |\n"
	"| -destDir | Output directory where the results are stored. | required |\n"
	"| -prefix | Prefix for output file names. A dot is inserted before the output suffix. | input file name |\n"
	"( Existing files with the same output names are overwritten; files with other prefixes are kept.)\n"
	"| -minLen | Minimum length of complete LTR-RTs. Constrains scoring system and filters. | 400 |\n"
	"| -maxLen |  Maximum length of complete LTR-RTs. Constrains scoring system and filters. | 22000 |\n"
	"| -minLenLTR | Minimum length of LTR direct repeat. Constrains filters. | 100 |\n"
	"| -maxLenLTR | Maximum length of LTR direct repeat. Constrains filters. | 6000 |\n"
	"( Note run time is highly dependent on this parameter, as it provides an upper bound for alignment length in the boundary adjustment step.)\n"               
	"| -id | Minimum identity [0-100] between 5' and 3' LTRs. | 85 |\n"
	"| -k  | Length of k-mers to adjust scoring system. Tradeoff between noise and resistance to mutation. | 13 |\n"
	"| -plateauSeed | Minimum length of plateaus to be initially considered 'Keep' in merging step. | 10 |\n"
	"| -nThreads | Number of cores to be used. | 1 |\n"
	"| -gapTol | Number of base pairs that two plateaus can differ by in height/distance. Affects both plateau merging and pairing steps. | 200 |\n"
	"|-seqLevel| Forces parallel execution on sequences within multi-FASTA file. Loads all sequences into memory | disabled |\n"
	"| -rawScores | prints the raw scores to a file called xxxxRawScores.csv under the output directory. | disabled |\n"
	"| -cleanedScores | prints the scores after merging to a file called xxxxCleanedScores.csv under the output directory. | disabled |\n"
	"| -nested | searches for nested elements. Results are stored in seperate files (marked as xxxxNestedDetector.bed) under the output directory | disabled |\n"
	"| -bedFormat | prints BED format without additional annotations (PPT and TSD). | disabled |\n"
	"| -help | prints this help message. | disabled |\n"
	;

	std::vector<string> * argList = new vector<string>();

	for(int i = 1;i<argc;i++){

		argList->push_back(string(argv[i]));
	}

	auto help = std::find(argList->begin(),argList->end(),"-help");

	if(help!=argList->end()){
		cout <<helpMessage<<endl;
		exit(1);
	}

	auto printRaw = std::find(argList->begin(),argList->end(),"-rawScores");

	if(printRaw!=argList->end()){

		printRawScores = true;
	}

	auto printClean = std::find(argList->begin(),argList->end(),"-cleanedScores");

		if(printClean!=argList->end()){

		printCleanScores = true;
	}


	auto bed = std::find(argList->begin(),argList->end(),"-bedFormat");

			if(bed!=argList->end()){

			bedOutput= true;
			cout <<"Bed format activated"<<endl;
		}

	auto nest = std::find(argList->begin(),argList->end(),"-nested");

			if(nest!=argList->end()){

			displayNested= true;
			
		}

	auto seq = std::find(argList->begin(),argList->end(),"-seqFile");
	auto src = std::find(argList->begin(),argList->end(),"-fasta");
	auto legacySrc = std::find(argList->begin(),argList->end(),"-chromDir");

	bool hasSeqFile = seq!=argList->end() && ++seq !=argList->end();
	bool hasFasta = src!=argList->end() && ++src !=argList->end();
	bool hasChromDir = legacySrc!=argList->end() && ++legacySrc !=argList->end();

	if((hasSeqFile ? 1 : 0) + (hasFasta ? 1 : 0) + (hasChromDir ? 1 : 0) != 1){
		cout<<" Specify exactly one of '-seqFile', '-fasta', or '-chromDir', followed by a path"<<endl;
		exit(1);
	}

	if(hasSeqFile){
		seqFile = *seq;
	}else if(hasFasta){
		chromDir = *src;
	}else{
		chromDir = *legacySrc;
	}

	auto dest = std::find(argList->begin(),argList->end(),"-destDir");

	if(dest!=argList->end() && ++dest != argList->end()){

		outputDir = *dest;
	}
	else{
		cout<<" \'-destDir\' is a required argument"<<endl;
		exit(1);
	}

	auto prefix = std::find(argList->begin(),argList->end(),"-prefix");
	if(prefix!=argList->end()){
		if(++prefix == argList->end()){
			cout<<" '-prefix' requires a value"<<endl;
			exit(1);
		}
		outputPrefix = *prefix;
	}

	int test = getOptionalArg(argList, "-minLen");

	if(test!=-1){
		minDist = test;
	}

	test = getOptionalArg(argList, "-maxLen");

	if(test!=-1){
		maxDist = test;
	}

	test = getOptionalArg(argList,"-minLenLTR");

	if(test!=-1){

		minLenLTR = test;
	}

	test = getOptionalArg(argList, "-maxLenLTR");

	if(test!=-1){
		maxLenLTR = test;
	}

	test = getOptionalArg(argList,"-id");

	if(test!= -1){

		identity = test;
	}

	test = getOptionalArg(argList, "-k");

	if(test!=-1){
		k = test;
	}

	test = getOptionalArg(argList, "-plateauSeed");

	if(test!=-1){
		minPlateauLength = test;
	}

	test = getOptionalArg(argList, "-gapTol");

	if(test!=-1){
		gapTol = test;
	}

	test = getOptionalArg(argList,"-nThreads");

	if(test!=-1){
		threads = test;
	}

	auto level = std::find(argList->begin(),argList->end(),"-seqLevel");

			if(level!=argList->end() && threads >1){

			seqLevel= true;
			
		}


	vector <string> * chromList = new vector<string>();

	if(!fs::exists(outputDir)){

		bool status = fs::create_directory(outputDir);

		if(!status){
			cerr<<"Output Directory could not be made"<<endl;
			throw std::exception();
		}
	}

	if(!seqFile.empty()){
		cout << "Reading sequence file ..." << endl;
		chromList->push_back(seqFile);
	}else{
		cout << "Reading chromosome directory ..." << endl;
		Util::readChromList(chromDir, chromList,"fa");
	}

	ChromosomeOneDigit *chrom;

	cout<<"Num threads: "<<threads<<endl;
	cout<<minDist<<endl;
	cout<<maxDist<<endl;
	cout<<minLenLTR<<endl;
	cout<<maxLenLTR<<endl;
	cout<<identity<<endl;
	cout<<k<<endl;
	cout<<minPlateauLength<<endl;
	cout<<gapTol<<endl;
	
	if(!seqLevel){
		#pragma omp parallel for schedule(dynamic) num_threads(threads)
		for (int i = 0; i < chromList->size(); i++)
		{

			string chromFile = chromList->at(i);

			int nameBegin = chromFile.find_last_of("/\\")+1;
			int nameEnd = chromFile.find_last_of(".") ;
			int len = nameEnd - nameBegin ;

			string name = chromFile.substr(nameBegin,len);

			string outFile = outputDir+"/"+(outputPrefix.empty() ? name : outputPrefix + ".");
			ChromListMaker chromListMaker(chromFile);
			const vector<ChromosomeOneDigit *> *  chromList = chromListMaker.makeChromOneDigitList();
			for(auto chrom : *chromList){
				TrCollector * collector = new TrCollector(chrom, outFile,name, minDist, maxDist, minLenLTR, maxLenLTR,identity,k, minPlateauLength, gapTol,printCleanScores,printRawScores,bedOutput,displayNested);
				
				delete collector;
			}
		}
	}else{
		for (int i = 0; i < chromList->size(); i++)
		{

			string chromFile = chromList->at(i);

			int nameBegin = chromFile.find_last_of("/\\")+1;
			int nameEnd = chromFile.find_last_of(".") ;
			int len = nameEnd - nameBegin ;

			string name = chromFile.substr(nameBegin,len);

			string outFile = outputDir+"/"+(outputPrefix.empty() ? name : outputPrefix + ".");
			ChromListMaker chromListMaker(chromFile);
			const vector<ChromosomeOneDigit *> *  chromList = chromListMaker.makeChromOneDigitList();
			#pragma omp parallel for schedule(dynamic) num_threads(threads)
			for(unsigned int j = 0; j < chromList->size(); j++){
				auto chrom = chromList->at(j);
				TrCollector * collector = new TrCollector(chrom, outFile,name, minDist, maxDist, minLenLTR, maxLenLTR,identity,k, minPlateauLength, gapTol,printCleanScores,printRawScores,bedOutput,displayNested);
				
				delete collector;
			}
	}

	chromList->clear();
	delete chromList;

	return 0;
	}
}

