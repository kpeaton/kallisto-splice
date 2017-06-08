#ifndef KALLISTO_PROCESSREADS_H
#define KALLISTO_PROCESSREADS_H

#include <zlib.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

#if defined(WIN32) || defined(WIN64)
#include "wincompat.h"
#endif
#include "csvread.h"
#include "exonmap.h"

#include "MinCollector.h"

#include "common.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

int ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc);
int ProcessBatchReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc, std::vector<std::vector<int>> &batchCounts);

class SequenceReader {
public:

  SequenceReader(const ProgramOptions& opt) :
  fp1(0),fp2(0),seq1(0),seq2(0),
  l1(0),l2(0),nl1(0),nl2(0),
  paired(!opt.single_end), files(opt.files),
  f_umi(new std::ifstream{}),
  current_file(0), state(false) {}
  SequenceReader() :
  fp1(0),fp2(0),seq1(0),seq2(0),
  l1(0),l2(0),nl1(0),nl2(0),
  paired(false), 
  f_umi(new std::ifstream{}),
  current_file(0), state(false) {}
  SequenceReader(SequenceReader&& o);
  
  bool empty();
  ~SequenceReader();

  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<std::string>& umis, 
                      bool full=false);

public:
  gzFile fp1 = 0, fp2 = 0;
  kseq_t *seq1 = 0, *seq2 = 0;
  int l1,l2,nl1,nl2;
  bool paired;
  std::vector<std::string> files;
  std::vector<std::string> umi_files;
  std::unique_ptr<std::ifstream> f_umi;
  int current_file;
  bool state; // is the file open
};

class MasterProcessor {
public:
  MasterProcessor (KmerIndex &index, const ProgramOptions& opt, MinCollector &tc)
    : tc(tc), index(index), opt(opt), SR(opt), numreads(0)
    ,nummapped(0), num_umi(0), tlencount(0), biasCount(0), maxBiasCount((opt.bias) ? 1000000 : 0) { 
      if (opt.batch_mode) {
        batchCounts.resize(opt.batch_ids.size(), {});
        
        for (auto &t : batchCounts) {
          t.resize(tc.counts.size(),0);
        }
        newBatchECcount.resize(opt.batch_ids.size());
        newBatchECumis.resize(opt.batch_ids.size());
        batchUmis.resize(opt.batch_ids.size());
      }

	  // Initialize the exon coordinate map
	  if ((opt.pseudobam) && (!opt.exon_coords_file.empty())) {

		  // Load the exon coordinate map data
		  std::ifstream       file(opt.exon_coords_file);
		  CSVRow              row;
		  std::string         last_key;
		  while (file >> row) {
			  if (last_key == row[0]){
				  std::get<2>(exon_map[last_key]).emplace_back(std::vector<int>({ std::stoi(row[1]), std::stoi(row[2]), std::stoi(row[3]) }));
			  } else {
				  last_key = row[0];
				  exon_map.emplace(last_key, std::make_tuple(row[4], row[5], std::vector<std::vector<int>>({ std::vector<int>({ std::stoi(row[1]), std::stoi(row[2]), std::stoi(row[3]) }) })));
			  }
		  }

		  // Collect information for the bam header
		  std::map<std::string, int> header_list;
		  for (int i = 0; i < index.num_trans; i++) {

			  const char * key = index.target_names_[i].c_str();
			  ExonMap::const_iterator map_entry = exon_map.find(key);
			  if (map_entry == exon_map.end()) {
				  std::cerr << "Transcript name could not be found in exon coordinate file: " << key << std::endl;
				  exit(1);
			  }

			  auto &entry = std::get<2>(map_entry->second);
			  int len;
			  if (entry[0][2] < 0) {
				  len = entry.front()[1] - entry.front()[0] - entry.front()[2];
			  } else {
				  len = entry.back()[1] - entry.back()[0] + entry.back()[2];
			  }

			  key = std::get<1>(map_entry->second).c_str();
			  if (header_list.find(key) == header_list.end()) {
				  header_list.emplace(key, len);
			  } else {
				  header_list[key] = (std::max)(header_list[key], len);
			  }

		  }

		  // Write the bam header
		  std::cout << "@HD\tVN:1.0\n";
		  for (auto const &entry : header_list) {
			  std::cout << "@SQ\tSN:" << entry.first << "\tLN:" << entry.second << "\n";
		  }
		  std::cout << "@PG\tID:kallisto\tPN:kallisto\tVN:" << KALLISTO_VERSION << "\n";
		  std::cout.flush();
	  }
    }

  std::mutex reader_lock;
  std::mutex writer_lock;

  SequenceReader SR;
  MinCollector& tc;
  KmerIndex& index;
  const ProgramOptions& opt;
  int numreads;
  int nummapped;
  int num_umi;
  std::atomic<int> tlencount;
  std::atomic<int> biasCount;
  std::vector<std::vector<int>> batchCounts;
  const int maxBiasCount;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> newECcount;
  std::vector<std::unordered_map<std::vector<int>, int, SortedVectorHasher>> newBatchECcount;
  std::vector<std::vector<std::pair<int, std::string>>> batchUmis;
  std::vector<std::vector<std::pair<std::vector<int>, std::string>>> newBatchECumis;

  ExonMap exon_map;

  void processReads();

  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& newEcs, std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi, int n, std::vector<int>& flens, std::vector<int> &bias, int id = -1);
};

class ReadProcessor {
public:
  ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int id = -1);
  ReadProcessor(ReadProcessor && o);
  ~ReadProcessor();
  char *buffer;
  size_t bufsize;
  bool paired;
  const MinCollector& tc;
  std::vector<std::pair<int, std::string>> ec_umi;
  std::vector<std::pair<std::vector<int>, std::string>> new_ec_umi;
  const KmerIndex& index;
  MasterProcessor& mp;
  SequenceReader batchSR;
  int numreads;
  int id;

  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<std::string> umis;
  std::vector<std::vector<int>> newEcs;
  std::vector<int> flens;
  std::vector<int> bias5;

  std::vector<int> counts;

  void operator()();
  void processBuffer();
  void clear();
};




#endif // KALLISTO_PROCESSREADS_H
