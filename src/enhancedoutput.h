#ifndef ENHANCEDOUTPUT_H
#define ENHANCEDOUTPUT_H

// Includes
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#ifdef _WIN32
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

// Typedefs
typedef std::map<std::string, std::tuple<std::string, std::string, std::vector<std::vector<int>>>> ExonMap;

// Classes
class EnhancedOutput
{
public:
	EnhancedOutput(KmerIndex &index, const ProgramOptions& opt);
	EnhancedOutput();

	// Exon coordinate map
	ExonMap exon_map;

	// SAM/BAM output data
	bool sortedbam;
	bool outputunmapped;
	std::map<std::string, std::vector<int>> ref_seq_map;
	std::string sam_header;
	std::string sort_dir;
	std::map<std::string, std::fstream*> sort_file_map;

	// Temporary and buffer storage for mapping reads and BAM output
	std::set<std::string> gene_list;
	std::vector<uint> bam_cigar;
	uint cigar_len;
	char outBamBuffer[MAX_BAM_ALIGN_SIZE];

	bool getSamData(std::string &ref_name, char *cig, int& strand, bool mapped, int &posread, int &posmate, int slen1, int slen2);
	void outputBamAlignment(std::string ref_name, int posread, int flag, int slen, int posmate, int tlen, const char *n1, const char *seq, const char *qual, int nmap, int strand);
	void outputSortedBam();

	static int reg2bin(int beg, int end);
	static void packseq(const char *seq_in, char *seq_out, uint seq_len);
	static char encodeNucleotide(char);
};

class CSVRow
{
// Adapted from a Stack Overflow answer from user Loki Astari: http://stackoverflow.com/a/1120224/52738
public:
	std::string const& operator[](std::size_t index) const;
	std::size_t size() const;
	void readNextRow(std::istream& str);
private:
	std::vector<std::string> m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data);

// Templates
template <class arrayType, int arraySize>
inline int funCompareArrays(const void *a, const void *b)
{
	arrayType* va = (uint*)a;
	arrayType* vb = (uint*)b;

	for (int i = 0; i < arraySize; i++) {
		if (va[i] > vb[i]) {
			return 1;
		} else if (va[i] < vb[i]) {
			return -1;
		};
	};

	return 0;

};

#endif // ENHANCEDOUTPUT_H