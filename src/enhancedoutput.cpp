#include "enhancedoutput.h"

HighResTimer::HighResTimer()
{
#if defined(_MSC_VER)
	QueryPerformanceFrequency(&frequency);
#endif
	resetTimer();
}

void HighResTimer::resetTimer()
{
#if defined(_MSC_VER)
	QueryPerformanceCounter(&reset_time);
#else
	reset_time = std::chrono::high_resolution_clock::now();
#endif
	previous_time = reset_time;
}

HighResTimer::duration HighResTimer::timeSinceReset()
{
#if defined(_MSC_VER)
	QueryPerformanceCounter(&current_time);
	duration elapsed = std::chrono::microseconds(((current_time.QuadPart - reset_time.QuadPart) * 1000000) / frequency.QuadPart);
#else
	current_time = std::chrono::high_resolution_clock::now();
	duration elapsed = current_time - reset_time;
#endif
	previous_time = current_time;
	return elapsed;
}

HighResTimer::duration HighResTimer::timeSincePrevious()
{
#if defined(_MSC_VER)
	QueryPerformanceCounter(&current_time);
	duration elapsed = std::chrono::microseconds(((current_time.QuadPart - previous_time.QuadPart) * 1000000) / frequency.QuadPart);
#else
	current_time = std::chrono::high_resolution_clock::now();
	duration elapsed = current_time - previous_time;
#endif
	previous_time = current_time;
	return elapsed;
}

EnhancedOutput::EnhancedOutput(KmerIndex &index, const ProgramOptions& opt)
	: enhancedoutput(opt.pseudobam && !opt.exon_coords_file.empty()), sortedbam(opt.sortedbam), outputunmapped(false),
	  get_sam_time(0), out_align_time(0), output_time(0), pre_sort_time(0), sort_time(0), post_sort_time(0)
{
	if (enhancedoutput) {

		// Load the exon coordinate map data
		std::ifstream file(opt.exon_coords_file);
		CSVRow row;
		std::string last_key;
		while (file >> row) {
			if (last_key == row[0]){
				std::get<2>(exon_map[last_key]).emplace_back(std::vector<int>({ std::stoi(row[1]), std::stoi(row[2]), std::stoi(row[3]) }));
			} else {
				last_key = row[0];
				exon_map.emplace(last_key, std::make_tuple(row[4], row[5], std::vector<std::vector<int>>({ std::vector<int>({ std::stoi(row[1]), std::stoi(row[2]), std::stoi(row[3]) }) })));
			}
		}
		file.close();

		// Collect information for the bam header
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
			if (ref_seq_map.find(key) == ref_seq_map.end()) {
				ref_seq_map.emplace(key, std::vector<int>({ len, -1, 0 }));
			} else {
				ref_seq_map[key][0] = (std::max)(ref_seq_map[key][0], len);
			}

		}

		if (sortedbam) {

			// Build header and initialize output file streams for sorting
			std::ostringstream header_stream;
			header_stream << "@HD\tVN:1.0\tSO:coordinate\n";

			sort_dir = opt.output + "/sorting";
			mkdir(sort_dir.c_str(), 0777);  // Move this to CheckOptions?

			int refid = 0;
			for (auto &entry : ref_seq_map) {

				header_stream << "@SQ\tSN:" << entry.first << "\tLN:" << entry.second[0] << "\n";
				sort_file_map[entry.first] = new std::fstream(sort_dir + "/" + entry.first, std::fstream::in | std::fstream::out | std::fstream::trunc | std::fstream::binary);
				entry.second[1] = refid++;

			}

			header_stream << "@PG\tID:kallisto\tPN:kallisto\tVN:" << KALLISTO_VERSION << "\n";
			sam_header = header_stream.str();

		} else {

			// Write the bam header
			std::cout << "@HD\tVN:1.0\n";
			for (auto const &entry : ref_seq_map) {
				std::cout << "@SQ\tSN:" << entry.first << "\tLN:" << entry.second[0] << "\n";
			}
			std::cout << "@PG\tID:kallisto\tPN:kallisto\tVN:" << KALLISTO_VERSION << "\n";
			std::cout.flush();

		}

	} else if (opt.pseudobam) {

		// Write the bam header
		index.writePseudoBamHeader(std::cout);

	}
}

EnhancedOutput::EnhancedOutput()
	: enhancedoutput(false), sortedbam(false), outputunmapped(false),
	  get_sam_time(0), out_align_time(0), output_time(0), pre_sort_time(0), sort_time(0), post_sort_time(0)
{
}

bool EnhancedOutput::getSamData(std::string &ref_name, char *cig, int& strand, bool mapped, int &posread, int &posmate, int slen1, int slen2)
{
	timer.resetTimer();
	bam_cigar.clear();
	align_len = 0;

	ExonMap::const_iterator map_entry = exon_map.find(ref_name);
	if (map_entry == exon_map.end()) {
		std::cerr << "Transcript name could not be found in exon coordinate file: " << ref_name << std::endl;
		exit(1);
	}

	const char *gene_name = std::get<0>(map_entry->second).c_str();
	if (gene_list.find(gene_name) == gene_list.end()) {
		gene_list.emplace(gene_name);
	} else {
		get_sam_time += timer.timeSinceReset();
		return 1;
	}

	ref_name = std::get<1>(map_entry->second);
	strand = (std::get<2>(map_entry->second)[0][2] < 0) ? -1 : 1;  // Apply elsewhere!

	std::string cig_string;
	uint op_len;
	int read_len = slen1;
	int mate_len = slen2;
	int read_offset = 0;
	int mate_offset = 0;
	for (auto &span : std::get<2>(map_entry->second)) {

		if (read_len != 0) {  // Not done mapping coordinates for read
			if (read_len != slen1) {  // In the process of mapping coordinates for read

				if (span[2] < 0) {  // Negative strand               Use trsense from findPosition instead??!!

					op_len = span[2] - span[1] + span[0] - 1 - read_offset;
					read_offset = span[1] - span[0] + 1;
					buildCigar(cig_string, 1, op_len, 'N', 3);
					if (read_len > read_offset) {
						op_len = read_offset;
						read_len -= read_offset;
						read_offset = span[2];
					} else {
						op_len = read_len;
						posread = read_offset - span[2] - read_len;
						read_len = 0;
					}
					buildCigar(cig_string, 1, op_len, 'M', 0);

				} else {  // Positive strand

					op_len = span[2] - read_offset - 1;
					read_offset = span[1] - span[0] + 1;
					buildCigar(cig_string, 0, op_len, 'N', 3);
					if (read_len > read_offset){
						op_len = read_offset;
						read_len -= read_offset;
						read_offset += span[2] - 1;
					} else {
						op_len = read_len;
						read_len = 0;
					}
					buildCigar(cig_string, 0, op_len, 'M', 0);

				}

			} else if (posread <= span[1]) {  // Begin mapping coordinates for read  // MOVE OUTSIDE LOOP?!

				if (posread < span[0]) {

					read_offset = span[0] - posread;
					op_len = read_offset;
					read_len -= read_offset;
					buildCigar(cig_string, 0, op_len, 'S', 4);
					
				}

				if (span[2] < 0) {  // Negative strand

					if ((posread + slen1 - 1) > span[1]) {
						read_offset = span[1] - posread - read_offset + 1;
						op_len = read_offset;
						read_len -= read_offset;
						read_offset = span[2];
					} else {
						op_len = read_len;
						read_len = 0;
						posread = span[1] - span[2] - posread - slen1 + 1;
					}
					buildCigar(cig_string, 1, op_len, 'M', 0);

				} else {  // Positive strand

					if ((posread + slen1 - 1) > span[1]) {
						read_offset = span[1] - posread - read_offset + 1;
						op_len = read_offset;
						read_len -= read_offset;
						read_offset = span[2] + span[1] - span[0];
					} else {
						op_len = read_len;
						read_len = 0;
					}
					posread += span[2] - span[0];
					buildCigar(cig_string, 0, op_len, 'M', 0);

				}
			}
		}

		if (mate_len != 0) {  // Not done mapping coordinates for mate
			if (mate_len != slen2) {  // In the process of mapping coordinates for mate

				if (span[2] < 0) {  // Negative strand

					mate_offset = span[1] - span[0] + 1;
					if (mate_len > mate_offset) {
						mate_len -= mate_offset;
						mate_offset = span[2];
					} else {
						posmate = mate_offset - span[2] - mate_len;
						mate_len = 0;
					}

				} else {  // Positive strand

					mate_offset = span[1] - span[0] + 1;
					if (mate_len > mate_offset){
						mate_len -= mate_offset;
						mate_offset += span[2] - 1;
					} else {
						mate_len = 0;
					}

				}

			} else if (posmate <= span[1]) {  // Begin mapping coordinates for mate

				if (posmate < span[0]) {
					mate_offset = span[0] - posmate;
					mate_len -= mate_offset;
				}

				if (span[2] < 0) {  // Negative strand

					if ((posmate + slen2 - 1) > span[1]) {
						mate_len -= span[1] - posmate - mate_offset + 1;
						mate_offset = span[2];
					} else {
						mate_len = 0;
						posmate = span[1] - span[2] - posmate - slen2 + 1;
					}

				} else {  // Positive strand

					if ((posmate + slen2 - 1) > span[1]) {
						mate_len -= span[1] - posmate - mate_offset + 1;
						mate_offset = span[2] + span[1] - span[0];
					} else {
						mate_len = 0;
					}
					posmate += span[2] - span[0];

				}

			}
		}

		if ((read_len == 0) && (mate_len == 0)) {  // Both have been mapped
			break;
		}

	}

	if (read_len != 0) {  // Account for read overhangs
		if (read_offset < 0) {  // Negative strand

			op_len = read_len;
			posread = -read_offset - read_len;
			buildCigar(cig_string, 1, op_len, 'S', 4);

		} else {  // Positive strand
			
			op_len = read_len;
			buildCigar(cig_string, 0, op_len, 'S', 4);

		}
	}

	if ((mate_len != 0) && (mate_offset < 0)) {  // Account for mate overhangs on negative strand
		posmate = -mate_offset - mate_len;
	}

	if (!sortedbam) {
		if (mapped) {
			sprintf(cig, "%s", cig_string.c_str());
		} else {
			sprintf(cig, "*");
		}
	}

	get_sam_time += timer.timeSinceReset();
	return 0;
}

void EnhancedOutput::buildCigar(std::string &cig_string, bool prepend, uint op_len, const char cig_char, uint cig_int)
{
	if (sortedbam) {

		if (prepend) {
			bam_cigar.insert(bam_cigar.begin(), ((op_len << 4) | cig_int));
		} else {
			bam_cigar.push_back(((op_len << 4) | cig_int));
		}
		if (cig_int != 4) {
			align_len += op_len;
		}

	} else {

		static char cig_[10];
		char *cig = &cig_[0];

		sprintf(cig, "%d%c", op_len, cig_char);
		if (prepend) {
			cig_string.insert(0, cig);
		} else {
			cig_string += cig;
		}
	}
}

void EnhancedOutput::outputBamAlignment(std::string ref_name, int posread, int flag, int slen, int posmate, int tlen, const char *n1, const char *seq, const char *qual, int nmap, int strand)
{
	timer.resetTimer();
	uint *buffer = (uint*)(outBamBuffer);
	uint n_bytes = 0;
	uint name_len;
	uint mapq = 255;
	uint n_cigar;
	uint cigar_bytes;

	// refID
	buffer[1] = ref_seq_map[ref_name][1];

	// pos
	buffer[2] = posread - 1;

	// bin_mq_nl
	name_len = strlen(n1) + 1;
	buffer[3] = ((reg2bin(posread - 1, posread + align_len - 1)<<16) | (mapq<<8) | name_len);

	// flag_nc
	n_cigar = bam_cigar.size();
	buffer[4] = ((flag << 16) | n_cigar);

	// l_seq
	buffer[5] = slen;

	// next_refID
	buffer[6] = ref_seq_map[ref_name][1];

	// next_pos
	buffer[7] = posmate - 1;

	// tlen
	buffer[8] = tlen;

	n_bytes = 9 * sizeof(uint);

	// read_name
	memcpy(outBamBuffer + n_bytes, n1, name_len);
	n_bytes += name_len;

	// cigar
	cigar_bytes = n_cigar*sizeof(uint);
	memcpy(outBamBuffer + n_bytes, bam_cigar.data(), cigar_bytes);
	n_bytes += cigar_bytes;

	// seq
	packseq(seq, outBamBuffer + n_bytes, slen);
	n_bytes += (slen + 1)/2;

	// qual
	for (uint i = 0; i < slen; i++) {
		(outBamBuffer + n_bytes)[i] = qual[i] - 33;
	};
	n_bytes += slen;

	// attributes
	memcpy(outBamBuffer + n_bytes, "NHi", 3);
	int value = 1;
	memcpy(outBamBuffer + n_bytes + 3, &value, sizeof(int));
	n_bytes += 3 + sizeof(int);
	if (strand < 0) {
		memcpy(outBamBuffer + n_bytes, "XSA-", 4);
	} else {
		memcpy(outBamBuffer + n_bytes, "XSA+", 4);
	}
	n_bytes += 4;

	// block_size
	buffer[0] = n_bytes - sizeof(uint);

	// Output to sorting file
	sort_file_map[ref_name]->write(outBamBuffer, n_bytes);
	ref_seq_map[ref_name][2]++; // Combine with sort_file_map?

	out_align_time += timer.timeSinceReset();
}

void EnhancedOutput::outputSortedBam()
{
	timer.resetTimer();
	std::fstream *sort_file;
	char *align_buffer;
	uint *sort_buffer;

	// Connect BAM output to stdout
	BGZF *bam_stream;
#ifdef _WIN32
	int result = _setmode(_fileno(stdout), _O_BINARY);
#endif
//	std::setvbuf(stdout, NULL, _IOFBF, 65536);
	bam_stream = bgzf_dopen(fileno(stdout), "w1");

	//// ABORT!!!
	//for (auto &entry : ref_seq_map) {
	//	sort_file_map[entry.first]->close();
	//	delete sort_file_map[entry.first];
	//	remove((sort_dir + "/" + entry.first).c_str());
	//}
	//bgzf_flush(bam_stream);
	//bgzf_close(bam_stream);
	//return;

	// Output header
	bgzf_write(bam_stream, "BAM\001", 4);
	int hlen = (int)sam_header.size();
	bgzf_write(bam_stream, (char*)&hlen, sizeof(hlen));
	bgzf_write(bam_stream, sam_header.c_str(), hlen);
//	std::cerr << hlen << std::endl << sam_header.c_str() << std::endl;
	int nchr = (int)ref_seq_map.size();
	bgzf_write(bam_stream, (char*)&nchr, sizeof(nchr));
	for (auto const &entry : ref_seq_map) {
		int namelen = (int)(entry.first.size() + 1);
		int seqlen = entry.second[0];
		//std::cerr << namelen << ": " << entry.first.data() << " " << seqlen << std::endl;
		bgzf_write(bam_stream, (char*)&namelen, sizeof(namelen));
		bgzf_write(bam_stream, entry.first.data(), namelen);
		bgzf_write(bam_stream, (char*)&seqlen, sizeof(seqlen));
	}
	bgzf_flush(bam_stream);

	// Create buffers to store each set of alignments and sorting data
	uint max_file_size = 0;
	uint max_n_align = 0;
	for (auto &entry : ref_seq_map) {
		max_file_size = (std::max)(max_file_size, (uint)sort_file_map[entry.first]->tellg());
		max_n_align = (std::max)(max_n_align, (uint)entry.second[2]);
	}
	align_buffer = new char[max_file_size+1];
	sort_buffer = new uint[max_n_align * 2];

	pre_sort_time += timer.timeSinceReset();
	// Sort alignments for each sorting file
	for (auto &entry : ref_seq_map) {

		timer.resetTimer();
		// Get file stream
		sort_file = sort_file_map[entry.first];
		uint n_align = entry.second[2];
		//std::cerr << entry.first << ": " << n_align << std::endl;
		if (n_align == 0) {
			sort_file->close();
			delete sort_file;
			remove((sort_dir + "/" + entry.first).c_str());
			continue;
		}

		// Load data from file
		uint n_bytes = sort_file->tellg();  // uint64_t??
		sort_file->seekg(std::fstream::beg);
		sort_file->read(align_buffer, n_bytes);
		sort_file->close();
		delete sort_file;
		remove((sort_dir + "/" + entry.first).c_str());

		// Collect sorting data
		for (uint position = 0, i = 0; i < n_align; i++) {
			uint *buffer = (uint*)(align_buffer + position);
			sort_buffer[i * 2] = buffer[2];
			sort_buffer[i * 2 + 1] = position;
			position += buffer[0] + sizeof(uint);
		}

		pre_sort_time += timer.timeSinceReset();
		// Sort by genomic position
		qsort((void*) sort_buffer, n_align, sizeof(uint) * 2, funCompareArrays<uint, 2>);

		sort_time += timer.timeSincePrevious();
		// Output sorted alignments
		for (uint i = 0; i < n_align; i++) {
			char *buffer = align_buffer + sort_buffer[i * 2 + 1];
			bgzf_write(bam_stream, buffer, *((uint*)buffer) + sizeof(uint));
		}
		post_sort_time += timer.timeSincePrevious();

	}

	// Flush and close output and delete buffers
	bgzf_flush(bam_stream);
	bgzf_close(bam_stream);
	delete[] align_buffer;
	delete[] sort_buffer;

}

// calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
int EnhancedOutput::reg2bin(int beg, int end)
{
	--end;
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
	return 0;
}

// pack nucleotides for BAM
void EnhancedOutput::packseq(const char *seq_in, char *seq_out, uint seq_len)
{
	for (uint i = 0; i < seq_len/2; i++) {
		seq_out[i] = ((encodeNucleotide(seq_in[2 * i])<<4) | encodeNucleotide(seq_in[2 * i + 1]));
	};
	if (seq_len%2 == 1) {
		seq_out[seq_len/2] = (encodeNucleotide(seq_in[seq_len - 1])<<4);
	};
}

// encode nucleotides for packing
char EnhancedOutput::encodeNucleotide(char cc) // Use a std::map instead?
{
	switch (cc) {
		case ('=') : cc = 0; break;
		case ('A') : case ('a') : cc = 1; break;
		case ('C') : case ('c') : cc = 2; break;
		case ('M') : case ('m') : cc = 3; break;
		case ('G') : case ('g') : cc = 4; break;
		case ('R') : case ('r') : cc = 5; break;
		case ('S') : case ('s') : cc = 6; break;
		case ('V') : case ('v') : cc = 7; break;
		case ('T') : case ('t') : cc = 8; break;
		case ('W') : case ('w') : cc = 9; break;
		case ('Y') : case ('y') : cc = 10; break;
		case ('H') : case ('h') : cc = 11; break;
		case ('K') : case ('k') : cc = 12; break;
		case ('D') : case ('d') : cc = 13; break;
		case ('B') : case ('b') : cc = 14; break;
		case ('N') : case ('n') : cc = 15; break;
		default: cc = 15;
	};
	return cc;
};

std::string const& CSVRow::operator[](std::size_t index) const
{
	return m_data[index];
}

std::size_t CSVRow::size() const
{
	return m_data.size();
}

void CSVRow::readNextRow(std::istream& str)
{
	std::string line;
	std::getline(str, line);

	std::stringstream lineStream(line);
	std::string cell;

	m_data.clear();
	while (std::getline(lineStream, cell, ',')) {
		m_data.push_back(cell);
	}
	// This checks for a trailing comma with no data after it.
	if (!lineStream && cell.empty()) {
		// If there was a trailing comma then add an empty element.
		m_data.push_back("");
	}
}

std::istream& operator>>(std::istream& str, CSVRow& data)
{
	data.readNextRow(str);
	return str;
}