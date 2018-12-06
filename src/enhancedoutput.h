#ifndef ENHANCEDOUTPUT_H
#define ENHANCEDOUTPUT_H

// Includes
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>

#ifdef _WIN32
#include <Windows.h>
#include <fcntl.h>
#include <io.h>
#endif
#if defined(_MSC_VER)
#include <BaseTsd.h>
#ifndef ssize_t
typedef SSIZE_T ssize_t;
#endif
#endif
extern "C" {
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/sam.h"
}

#include "common.h"
#include "KmerIndex.h"

#include "wincompat.h"

// Constants
#define MAX_BAM_ALIGN_SIZE 32768
//std::map<char, char> nucleotide_map {
//	{ '=', 0 },
//	{ 'A', 1 }, { 'a', 1 },
//	{ 'C', 2 }, { 'c', 2 },
//	{ 'M', 3 }, { 'm', 3 },
//	{ 'G', 4 }, { 'g', 4 },
//	{ 'R', 5 }, { 'r', 5 },
//	{ 'S', 6 }, { 's', 6 },
//	{ 'V', 7 }, { 'v', 7 },
//	{ 'T', 8 }, { 't', 8 },
//	{ 'W', 9 }, { 'w', 9 },
//	{ 'Y', 10 }, { 'y', 10 },
//	{ 'H', 11 }, { 'h', 11 },
//	{ 'K', 12 }, { 'k', 12 },
//	{ 'D', 13 }, { 'd', 13 },
//	{ 'B', 14 }, { 'b', 14 },
//	{ 'N', 15 }, { 'n', 15 }
//};

// Classes

class HighResTimer
{
public:
	HighResTimer();

#if defined(_MSC_VER)
	typedef std::chrono::duration<double, std::ratio<1, 1>> duration;
#else
	typedef std::chrono::duration<double> duration;
#endif

	void resetTimer();
	duration timeSinceReset();
	duration timeSincePrevious();

private:
#if defined(_MSC_VER)
	LARGE_INTEGER reset_time, previous_time, current_time, frequency;
#else
	typedef std::chrono::time_point<std::chrono::high_resolution_clock> timePoint;
	timePoint reset_time, previous_time, current_time;
#endif
};

class EnhancedOutput
{
public:
	EnhancedOutput(KmerIndex &index, const ProgramOptions& opt);
	~EnhancedOutput();

	// Gene coordinate map:
	bool enhancedoutput;
	static enum IntronFlag {intronNone, intronStart, intronEnd, intronFull};
	typedef std::vector<int> SegmentData;
	typedef std::vector<SegmentData> SegmentArray;
	typedef std::tuple<std::string, std::string, int, IntronFlag, SegmentArray> GeneData;
	std::unordered_map<std::string, GeneData> gene_map;

	// General SAM/BAM output data:
	bool pseudobam;
	bool sortedbam;
	bool outputunmapped;
	std::map<std::string, std::vector<int>> ref_seq_map;
	std::string sam_header;
	std::string sort_dir;
	char* outBamBuffer;
	BGZF* bam_stream;
	int num_threads;
	int current_sorting_index;
	std::mutex sorting_lock;
	std::vector<std::vector<std::fstream>> sorting_streams;
	std::vector<std::vector<uint>> num_alignments;

	// BED output data:
	bool outputbed;
	std::string bed_file;
	typedef std::tuple<std::string, int, int, std::string> JunctionKey;
	typedef std::map<JunctionKey, std::tuple<int, IntronFlag, char, int, int, int>> JunctionMap;
	std::vector<JunctionMap> junction_map;

	// Methods:
	void processAlignment(std::string trans_name, int flag1, int posread, int slen1, const char *name1, const char *seq1, const char *qual1, int flag2, int posmate, int slen2, const char *name2, const char *seq2, const char *qual2, int nmap, int id);
	void buildBAMCigar(std::vector<uint> &bam_cigar, bool prepend, uint op_len, uint cig_int);
	void buildSAMCigar(std::string &sam_cigar, bool prepend, uint op_len, const char cig_char);
	void mapJunction(int id, std::string chrom_name, std::string trans_name, IntronFlag intron_flag, bool negstrand, int junction_coord, int pair_coord, int size1 = 0, int size2 = 0);
	void outputJunction();
	void outputSortedBam();
	void fetchChromosome();
	void sortChromosome(int ref_ID);
	void removeSortingFiles(int ref_ID);

	static int reg2bin(int beg, int end);
	static void packseq(const char *seq_in, char *seq_out, uint seq_len);
	static char encodeNucleotide(char);
};

class DelimRow
{
// Adapted from a Stack Overflow answer from user Loki Astari: http://stackoverflow.com/a/1120224/52738
public:
	DelimRow(char delimiter = ',');

	std::string const& operator[](std::size_t index) const;
	std::size_t size() const;
	void readNextRow(std::istream& str);

private:
	char delim;
	std::vector<std::string> m_data;
};

std::istream& operator>>(std::istream& str, DelimRow& data);

// Templates
template <class arrayType, int arraySize>
inline int funCompareArrays(const void *a, const void *b)
{
	arrayType* va = (arrayType*)a;
	arrayType* vb = (arrayType*)b;

	for (int i = 0; i < arraySize; i++) {
		if (va[i] > vb[i]) {
			return 1;
		} else if (va[i] < vb[i]) {
			return -1;
		};
	};

	return 0;

};

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

#endif // ENHANCEDOUTPUT_H