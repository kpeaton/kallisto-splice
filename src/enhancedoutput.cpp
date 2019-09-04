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
	: enhancedoutput(!opt.gene_coords_file.empty()),
	  pseudobam(opt.pseudobam),
	  sortedbam(opt.sortedbam),
	  outputunmapped(false),
	  num_threads(opt.threads),
	  current_sorting_index(-1),
	  outputbed(opt.outputbed),
	  bed_file(opt.bed_file)
{
	if (enhancedoutput) {
		
		// Open the coordinate map file
		std::ifstream file(opt.gene_coords_file);
		DelimRow row('\t');
		std::string last_key;
		int segment_end;
		auto last_entry = gene_map.end();
		
		// Valid chromosome ID set
		//std::set<std::string> valid_chromosomes = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
		//											"15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT" };

		// Read header line
		file >> row;

		// Load the exon/intron coordinate map data
		while (file >> row) {

			int segment_start;
			int genome_position = std::stoi(row[3]);
			int segment_length = std::stoi(row[4]) - genome_position + 1;

			if (last_key != row[7]) {  // New exon or intron segment

				std::string chromID = row[1];
				//if (valid_chromosomes.count(chromID) == 0) {
				//	continue;
				//}

				segment_start = 1;
				segment_end = segment_length;
				IntronFlag intron_flag = intronNone;
				int strand = std::stoi(row[2]);
				int str_pos = row[7].find("::");
				
				if (str_pos != std::string::npos) {  // Intron segment
					
					std::string short_key = row[7].substr(0, str_pos);
					
					if (last_key == short_key) {  // Intron end segment that pairs with previous intron start segment

						// Modify intron start segment data
						std::get<1>(last_entry->second) = short_key;
						std::get<3>(last_entry->second) = intronStart;
						
						// Create common segment data for start and end map entries to reference
						GeneData &new_entry = gene_map.emplace(short_key, std::forward_as_tuple(row[0], chromID, strand, intronFull, SegmentArray{})).first->second;
						std::get<4>(new_entry).swap(std::get<4>(last_entry->second));
						std::get<4>(new_entry).emplace_back(std::initializer_list<int>{ segment_start, segment_end, genome_position });

						// Create intron end segment map entry
						last_entry = gene_map.emplace(row[7], std::forward_as_tuple(row[0], short_key, strand, intronEnd, SegmentArray{})).first;
						continue;

					} else {  // Assume full intron segment (for now)

						intron_flag = intronFull;
						last_key = short_key;

					}

				} else {  // Exon segment

					last_key = row[7];

				}

				last_entry = gene_map.emplace(row[7], std::forward_as_tuple(row[0], chromID, strand, intron_flag, SegmentArray{})).first;

			} else {

				segment_start = segment_end + 1;
				segment_end += segment_length;

			}

			std::get<4>(last_entry->second).emplace_back(std::initializer_list<int>{ segment_start, segment_end, genome_position });

		}
		file.close();
		
		// Collect information for the bam header
		for (int tr = 0; tr < index.num_trans; tr++) {

			std::string key = index.target_names_[tr];
			auto map_entry = gene_map.find(key);
			if (map_entry == gene_map.end()) {
				std::cerr << "Transcript name could not be found in exon coordinate file: " << key << std::endl;
				exit(1);
			}

			GeneData gene_data = map_entry->second;
			int strand = std::get<2>(gene_data);
			IntronFlag intron_flag = std::get<3>(gene_data);
			if ((intron_flag == intronStart) || (intron_flag == intronEnd)) {
				gene_data = gene_map.at(std::get<1>(gene_data));
				strand = 1;
			}
			const SegmentArray &segments = std::get<4>(gene_data);

			int len;
			if (strand < 0) {
				len = segments.front()[2] + segments.front()[1] - segments.front()[0];
			} else {
				len = segments.back()[2] + segments.back()[1] - segments.back()[0];
			}

			key = std::get<1>(gene_data);
			auto ref_entry = ref_seq_map.find(key);
			if (ref_entry == ref_seq_map.end()) {
				ref_seq_map.emplace(key, std::initializer_list<int>{ len, -1 });
			} else {
				ref_entry->second[0] = (std::max)(ref_entry->second[0], len);
			}

		}

		if (sortedbam) {  // Build header and initialize output file streams for sorting

			std::ostringstream header_stream;
			header_stream << "@HD\tVN:1.0\tSO:coordinate\n";
			sort_dir = opt.output + "/sorting";
			mkdir(sort_dir.c_str(), 0777);  // Move this to CheckOptions?

			int num_chromosomes = (int)ref_seq_map.size();
			int ref_ID = 0;
			std::fstream::openmode stream_flags = std::fstream::in | std::fstream::out | std::fstream::trunc | std::fstream::binary;
			sorting_streams.reserve(num_threads);
			num_alignments.reserve(num_threads);

			for (auto &entry : ref_seq_map) {

				header_stream << "@SQ\tSN:" << entry.first << "\tLN:" << entry.second[0] << "\n";
				std::string sort_file = sort_dir + "/" + std::to_string(ref_ID) + "_";

				for (int thread = 0; thread < num_threads; thread++) {
					if (ref_ID == 0) {
						sorting_streams.emplace_back();
						sorting_streams[thread].reserve(num_chromosomes);
						num_alignments.emplace_back();
						num_alignments[thread].reserve(num_chromosomes);
					}
					sorting_streams[thread].emplace_back(sort_file + std::to_string(thread), stream_flags);
					num_alignments[thread].emplace_back(0);
				}

				entry.second[1] = ref_ID++;

			}

			header_stream << "@PG\tID:kallisto\tPN:kallisto\tVN:" << KALLISTO_VERSION << "\n";
			sam_header = header_stream.str();

			outBamBuffer = new char[MAX_BAM_ALIGN_SIZE*num_threads];

		} else if (pseudobam) {  // Write the bam header

			std::cout << "@HD\tVN:1.0\n";
			for (auto const &entry : ref_seq_map) {
				std::cout << "@SQ\tSN:" << entry.first << "\tLN:" << entry.second[0] << "\n";
			}
			std::cout << "@PG\tID:kallisto\tPN:kallisto\tVN:" << KALLISTO_VERSION << "\n";
			std::cout.flush();

		}

		if (outputbed) {
			junction_map.reserve(num_threads);
			for (int thread = 0; thread < num_threads; thread++) {
				junction_map.emplace_back();
			}
		}

	} else if (pseudobam) {  // Write the bam header

		index.writePseudoBamHeader(std::cout);

	}
}

EnhancedOutput::~EnhancedOutput()
{
	if (sortedbam && (outBamBuffer != nullptr)) {
		delete[] outBamBuffer;
		outBamBuffer = nullptr;
	}
}

void EnhancedOutput::processAlignment(std::string trans_name, int flag1, int posread, int slen1, const char *name1, const char *seq1, const char *qual1, int flag2, int posmate, int slen2, const char *name2, const char *seq2, const char *qual2, int nmap, int id)
{
	GeneData gene_data = gene_map.at(trans_name);
	int strand = (std::get<2>(gene_data) < 0) ? -1 : 1;
	IntronFlag intron_flag = std::get<3>(gene_data);
	SegmentArray::const_iterator iter_start, iter_end, intron_pair;

	switch (intron_flag) {
		case intronFull:
			trans_name = trans_name.substr(0, trans_name.find("::"));
		case intronNone:
		{
			const SegmentArray &segments = std::get<4>(gene_data);
			iter_start = segments.begin();
			iter_end = segments.end();
			break;
		}
		case intronStart:
		{
			trans_name = std::get<1>(gene_data);
			gene_data = gene_map.at(trans_name);
			const SegmentArray &segments = std::get<4>(gene_data);
			iter_start = segments.begin();
			iter_end = ++segments.begin();
			intron_pair = iter_end;
			break;
		}
		case intronEnd:
		{
			trans_name = std::get<1>(gene_data);
			gene_data = gene_map.at(trans_name);
			const SegmentArray &segments = std::get<4>(gene_data);
			intron_pair = segments.begin();
			iter_start = ++segments.begin();
			iter_end = segments.end();
			break;
		}
	}

	std::string ref_name = std::get<1>(gene_data);
	bool negstrand = strand < 0;
	std::vector<uint> read_bam_cigar;
	std::vector<uint> mate_bam_cigar;
	std::string read_sam_cigar;
	std::string mate_sam_cigar;
	uint read_len = 0;
	uint mate_len = 0;
	uint op_len;

	int read_rem = slen1;
	int read_mapq = 255;
	if (flag1 & 0x04) {
		posread = 0;
		read_rem = 0;
		read_mapq = 0;
	}

	int mate_rem = slen2;
	int mate_mapq = 255;
	if (flag2 & 0x04) {
		posmate = 0;
		mate_rem = 0;
		mate_mapq = 0;
	}

	int read_offset = 0;
	int mate_offset = 0;
	int start_coord = 0;
	int end_coord = 0;

	// Map alignment to chromosome coordinates
	for (auto span = iter_start; span != iter_end; span++) {

		const int segment_start = (*span)[0];
		const int segment_end = (*span)[1];
		const int genome_position = (*span)[2];

		if (read_rem > 0) {  // Not done mapping read

			if (read_rem < slen1) {  // In the process of mapping read

				// Find span of skipped region:
				if (negstrand) {
					start_coord = genome_position + segment_end - segment_start;
					end_coord = read_offset;
				} else {
					start_coord = read_offset;
					end_coord = genome_position;
				}

				// Update CIGAR string:
				op_len = end_coord - start_coord - 1;
				read_len += op_len;
				if (sortedbam) {
					buildBAMCigar(read_bam_cigar, negstrand, op_len, 3);
				} else if (pseudobam) {
					buildSAMCigar(read_sam_cigar, negstrand, op_len, 'N');
				}

				// Store BED file information:
				if (outputbed) {
					mapJunction(id, ref_name, trans_name, intron_flag, negstrand, start_coord, end_coord);
				}

				// Find overlap with segment:
				read_offset = segment_end - segment_start + 1;
				if (read_rem > read_offset){
					op_len = read_offset;
					read_offset = (negstrand) ? genome_position : read_offset + genome_position - 1;
				} else {
					op_len = read_rem;
					if (negstrand) {
						posread = start_coord - read_rem + 1;
					}
				}
				read_rem -= op_len;

				// Update CIGAR string:
				read_len += op_len;
				if (sortedbam) {
					buildBAMCigar(read_bam_cigar, negstrand, op_len, 0);
				} else if (pseudobam) {
					buildSAMCigar(read_sam_cigar, negstrand, op_len, 'M');
				}

			} else if (posread <= segment_end) {  // Begin mapping read

				// Find soft clipping at beginning of read:
				if (posread < segment_start) {
					op_len = segment_start - posread;
					read_rem -= op_len;
					posread += op_len;
					if (sortedbam) {
						buildBAMCigar(read_bam_cigar, false, op_len, 4);
					} else if (pseudobam) {
						buildSAMCigar(read_sam_cigar, false, op_len, 'S');
					}
				}

				// Find overlap with segment:
				read_offset = posread + read_rem - segment_end - 1;
				if (read_offset > 0) {
					op_len = read_rem - read_offset;
					if (negstrand) {
						read_offset = genome_position;
					} else {
						read_offset = genome_position + segment_end - segment_start;
						posread += genome_position - segment_start;
					}
				} else {
					op_len = read_rem;
					posread = (negstrand) ? genome_position - read_offset : posread + genome_position - segment_start;
				}
				read_rem -= op_len;

				// Update CIGAR string:
				read_len += op_len;
				if (sortedbam) {
					buildBAMCigar(read_bam_cigar, negstrand, op_len, 0);
				} else if (pseudobam) {
					buildSAMCigar(read_sam_cigar, negstrand, op_len, 'M');
				}

			}

		}

		if (mate_rem > 0) {  // Not done mapping mate

			if (mate_rem < slen2) {  // In the process of mapping mate

				// Find span of skipped region:
				if (negstrand) {
					start_coord = genome_position + segment_end - segment_start;
					end_coord = mate_offset;
				} else {
					start_coord = mate_offset;
					end_coord = genome_position;
				}

				// Update CIGAR string:
				op_len = end_coord - start_coord - 1;
				mate_len += op_len;
				if (sortedbam) {
					buildBAMCigar(mate_bam_cigar, negstrand, op_len, 3);
				} else if (pseudobam) {
					buildSAMCigar(mate_sam_cigar, negstrand, op_len, 'N');
				}

				// Store BED file information:
				if (outputbed) {
					mapJunction(id, ref_name, trans_name, intron_flag, negstrand, start_coord, end_coord);
				}

				// Find overlap with segment:
				mate_offset = segment_end - segment_start + 1;
				if (mate_rem > mate_offset){
					op_len = mate_offset;
					mate_offset = (negstrand) ? genome_position : mate_offset + genome_position - 1;
				} else {
					op_len = mate_rem;
					if (negstrand) {
						posmate = start_coord - mate_rem + 1;
					}
				}
				mate_rem -= op_len;

				// Update CIGAR string:
				mate_len += op_len;
				if (sortedbam) {
					buildBAMCigar(mate_bam_cigar, negstrand, op_len, 0);
				} else if (pseudobam) {
					buildSAMCigar(mate_sam_cigar, negstrand, op_len, 'M');
				}

			} else if (posmate <= segment_end) {  // Begin mapping mate

				// Find soft clipping at beginning of mate:
				if (posmate < segment_start) {
					op_len = segment_start - posmate;
					mate_rem -= op_len;
					posmate += op_len;
					if (sortedbam) {
						buildBAMCigar(mate_bam_cigar, false, op_len, 4);
					} else if (pseudobam) {
						buildSAMCigar(mate_sam_cigar, false, op_len, 'S');
					}
				}

				// Find overlap with segment:
				mate_offset = posmate + mate_rem - segment_end - 1;
				if (mate_offset > 0) {
					op_len = mate_rem - mate_offset;
					if (negstrand) {
						mate_offset = genome_position;
					} else {
						mate_offset = genome_position + segment_end - segment_start;
						posmate += genome_position - segment_start;
					}
				} else {
					op_len = mate_rem;
					posmate = (negstrand) ? genome_position - mate_offset : posmate + genome_position - segment_start;
				}
				mate_rem -= op_len;

				// Update CIGAR string:
				mate_len += op_len;
				if (sortedbam) {
					buildBAMCigar(mate_bam_cigar, negstrand, op_len, 0);
				} else if (pseudobam) {
					buildSAMCigar(mate_sam_cigar, negstrand, op_len, 'M');
				}

			}

		}

		if ((read_rem <= 0) && (mate_rem <= 0)) {  // Both have been mapped
			break;
		}

	}

	if (read_rem > 0) {  // Account for read overhangs
		if (negstrand) {
			posread = read_offset;
		}
		if (sortedbam) {
			buildBAMCigar(read_bam_cigar, negstrand, read_rem, 4);
		} else if (pseudobam) {
			buildSAMCigar(read_sam_cigar, negstrand, read_rem, 'S');
		}
	}

	if (mate_rem > 0) {  // Account for mate overhangs
		if (negstrand) {
			posmate = mate_offset;
		}
		if (sortedbam) {
			buildBAMCigar(mate_bam_cigar, negstrand, mate_rem, 4);
		} else if (pseudobam) {
			buildSAMCigar(mate_sam_cigar, negstrand, mate_rem, 'S');
		}
	}

	// Calculate tlen
	int tlen = 0;
	if ((posread != 0) && (posmate != 0)) {
		tlen = (std::max)(posread + read_len, posmate + mate_len) - (std::min)(posread, posmate);
	}

	// Store junction information
	if (outputbed && (intron_flag != intronNone)) {

		if (intron_flag != intronEnd) {  // Map intron start junction

			start_coord = (*iter_start)[2];
			if (intron_flag == intronFull) {
				end_coord = start_coord + (*iter_start)[1] - (*iter_start)[0];
			} else {
				end_coord = (*intron_pair)[2] + (*intron_pair)[1] - (*intron_pair)[0];
			}

			if ((posread >= start_coord) && (posread < start_coord + 50) &&
				(posread + slen1 >= start_coord + 50) && (posread + slen1 < end_coord - 50) &&
				(posmate < end_coord)) {  // Check posmate versus start_coord too?
				mapJunction(id, ref_name, trans_name, intronStart, negstrand, start_coord + 50, end_coord - 50, 10, 10);
			}

			if ((posmate >= start_coord) && (posmate < start_coord + 50) &&
				(posmate + slen2 >= start_coord + 50) && (posmate + slen2 < end_coord - 50) &&
				(posread < end_coord)) {  // Check posread versus start_coord too?
				mapJunction(id, ref_name, trans_name, intronStart, negstrand, start_coord + 50, end_coord - 50, 10, 10);
			}

		}

		if (intron_flag != intronStart) {  // Map intron end junction

			end_coord = (*iter_start)[2] + (*iter_start)[1] - (*iter_start)[0];
			if (intron_flag == intronFull) {
				start_coord = (*iter_start)[2];
			} else {
				start_coord = (*intron_pair)[2];
			}

			if ((posread >= start_coord + 50) && (posread < end_coord - 50) &&
				(posread + slen1 >= end_coord - 50) && (posread + slen1 < end_coord) &&
				(posmate + slen2 >= start_coord)) {  // Check posmate versus end_coord too?
				mapJunction(id, ref_name, trans_name, intronEnd, negstrand, end_coord - 50, start_coord + 50, 10, 10);
			}

			if ((posmate >= start_coord + 50) && (posmate < end_coord - 50) &&
				(posmate + slen2 >= end_coord - 50) && (posmate + slen2 < end_coord) &&
				(posread + slen1 >= start_coord)) {  // Check posread versus end_coord too?
				mapJunction(id, ref_name, trans_name, intronEnd, negstrand, end_coord - 50, start_coord + 50, 10, 10);
			}

		}

	}

	// Output alignment
	if (sortedbam) {

		char *threadBuffer = outBamBuffer + id*MAX_BAM_ALIGN_SIZE;
		uint *buffer = (uint*)(threadBuffer);
		uint n_bytes = 0;
		uint ref_ID = ref_seq_map.at(ref_name)[1];

		uint name_len;
		uint n_cigar;
		uint cigar_bytes;

		if (posread != 0) {

			name_len = strlen(name1) + 1;
			n_cigar = read_bam_cigar.size();
			cigar_bytes = n_cigar * sizeof(uint);

			buffer[1] = ref_ID;  // refID
			buffer[2] = posread - 1;  // pos
			buffer[3] = ((reg2bin(posread - 1, posread + read_len - 1) << 16) | (read_mapq << 8) | name_len);  // bin_mq_nl
			buffer[4] = ((flag1 << 16) | n_cigar);  // flag_nc
			buffer[5] = slen1;  // l_seq
			buffer[6] = ref_ID;  // next_refID
			buffer[7] = posmate - 1;  // next_pos
			buffer[8] = tlen;  // tlen
			n_bytes = 9 * sizeof(uint);

			// read_name
			memcpy(threadBuffer + n_bytes, name1, name_len);
			n_bytes += name_len;

			// cigar
			memcpy(threadBuffer + n_bytes, read_bam_cigar.data(), cigar_bytes);
			n_bytes += cigar_bytes;

			// seq
			packseq(seq1, threadBuffer + n_bytes, slen1);
			n_bytes += (slen1 + 1) / 2;

			// qual
			for (uint i = 0; i < slen1; i++) {
				(threadBuffer + n_bytes)[i] = qual1[i] - 33;
			};
			n_bytes += slen1;

			// attributes
			memcpy(threadBuffer + n_bytes, "NHi", 3);
			memcpy(threadBuffer + n_bytes + 3, &nmap, sizeof(int));
			n_bytes += 3 + sizeof(int);
			memcpy(threadBuffer + n_bytes, (strand < 0) ? "XSA-" : "XSA+", 4);
			//memcpy(threadBuffer + n_bytes, bool(flag & 0x10) ? "XSA-" : "XSA+", 4);  // Which is better to use?!
			n_bytes += 4;

			// block_size
			buffer[0] = n_bytes - sizeof(uint);

			// Output to sorting file
			sorting_streams[id][ref_ID].write(threadBuffer, n_bytes);
			num_alignments[id][ref_ID]++;

		}

		if (posmate != 0) {

			name_len = strlen(name2) + 1;
			n_cigar = mate_bam_cigar.size();
			cigar_bytes = n_cigar * sizeof(uint);

			buffer[1] = ref_ID;  // refID
			buffer[2] = posmate - 1;  // pos
			buffer[3] = ((reg2bin(posmate - 1, posmate + mate_len - 1) << 16) | (mate_mapq << 8) | name_len);  // bin_mq_nl
			buffer[4] = ((flag2 << 16) | n_cigar);  // flag_nc
			buffer[5] = slen2;  // l_seq
			buffer[6] = ref_ID;  // next_refID
			buffer[7] = posread - 1;  // next_pos
			buffer[8] = -tlen;  // tlen
			n_bytes = 9 * sizeof(uint);

			// read_name
			memcpy(threadBuffer + n_bytes, name2, name_len);
			n_bytes += name_len;

			// cigar
			memcpy(threadBuffer + n_bytes, mate_bam_cigar.data(), cigar_bytes);
			n_bytes += cigar_bytes;

			// seq
			packseq(seq2, threadBuffer + n_bytes, slen2);
			n_bytes += (slen2 + 1) / 2;

			// qual
			for (uint i = 0; i < slen2; i++) {
				(threadBuffer + n_bytes)[i] = qual2[i] - 33;
			};
			n_bytes += slen2;

			// attributes
			memcpy(threadBuffer + n_bytes, "NHi", 3);
			memcpy(threadBuffer + n_bytes + 3, &nmap, sizeof(int));
			n_bytes += 3 + sizeof(int);
			memcpy(threadBuffer + n_bytes, (strand < 0) ? "XSA-" : "XSA+", 4);
			//memcpy(threadBuffer + n_bytes, bool(flag & 0x10) ? "XSA-" : "XSA+", 4);  // Which is better to use?!
			n_bytes += 4;

			// block_size
			buffer[0] = n_bytes - sizeof(uint);

			// Output to sorting file
			sorting_streams[id][ref_ID].write(threadBuffer, n_bytes);
			num_alignments[id][ref_ID]++;

		}

	} else if (pseudobam) {

		if (posread != 0) {

			printf("%s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d", name1, flag1, ref_name.c_str(), posread, read_mapq, read_sam_cigar.c_str(), posmate, tlen, seq1, qual1, nmap);
			if (negstrand) {
				printf("\tXS:A:-\n");
			} else {
				printf("\tXS:A:+\n");
			}

		}

		if (posmate != 0) {

			printf("%s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d", name2, flag2, ref_name.c_str(), posmate, mate_mapq, mate_sam_cigar.c_str(), posread, -tlen, seq2, qual2, nmap);
			if (negstrand) {
				printf("\tXS:A:-\n");
			} else {
				printf("\tXS:A:+\n");
			}

		}

	}

}

void EnhancedOutput::buildBAMCigar(std::vector<uint> &bam_cigar, bool prepend, uint op_len, uint cig_int)
{
	if (prepend) {
		bam_cigar.insert(bam_cigar.begin(), ((op_len << 4) | cig_int));
	} else {
		bam_cigar.push_back(((op_len << 4) | cig_int));
	}
}

void EnhancedOutput::buildSAMCigar(std::string &sam_cigar, bool prepend, uint op_len, const char cig_char)
{
	static char cig_[10];
	char *cig = &cig_[0];

	sprintf(cig, "%d%c", op_len, cig_char);  // Would to_string be faster?
	if (prepend) {
		sam_cigar.insert(0, cig);
	} else {
		sam_cigar += cig;
	}
}

void EnhancedOutput::mapJunction(int id, std::string chrom_name, std::string trans_name, IntronFlag intron_flag, bool negstrand, int junction_coord, int pair_coord, int size1, int size2)
{
	JunctionKey key;
	switch (intron_flag) {
		case intronNone:
			key = std::make_tuple(chrom_name, junction_coord, pair_coord, trans_name);
			break;
		case intronStart:
			key = std::make_tuple(chrom_name, junction_coord - 11, junction_coord + 9, trans_name);
			break;
		case intronEnd:
			key = std::make_tuple(chrom_name, junction_coord - 10, junction_coord + 10, trans_name);
			break;
	}

	auto map_entry = junction_map[id].find(key);
	if (map_entry == junction_map[id].end()) {
		junction_map[id].emplace(key, std::forward_as_tuple(1, intron_flag, (negstrand) ? '-' : '+', size1, size2, pair_coord));
	} else {
		std::get<0>(map_entry->second)++;
	}
}

void EnhancedOutput::outputJunction()
{
	// Merge junction information from threads
	for (int thread = 1; thread < num_threads; thread++) {
		for (const auto &entry : junction_map[thread]) {
			JunctionKey key = entry.first;
			auto map_entry = junction_map[0].find(key);
			if (map_entry == junction_map[0].end()) {
				junction_map[0].emplace(key, entry.second);
			} else {
				std::get<0>(map_entry->second) += std::get<0>(entry.second);
			}
		}
	}

	// Output junctions
	std::fstream bed_out(bed_file.c_str(), std::fstream::out | std::fstream::trunc);
	for (const auto &entry : junction_map[0]) {

		JunctionKey key = entry.first;
		int start_coord = std::get<1>(key);
		int end_coord = std::get<2>(key);
		std::string trans_name = std::get<3>(key);

		const auto &junction_data = entry.second;
		IntronFlag intron_flag = std::get<1>(junction_data);

		if (intron_flag != intronNone) {

			int pair_coord = std::get<5>(junction_data);

			if (intron_flag == intronStart) {
				std::get<1>(key) = pair_coord - 10;
				std::get<2>(key) = pair_coord + 10;
				trans_name += '-' + std::to_string(start_coord + 11);
			} else {
				std::get<1>(key) = pair_coord - 11;
				std::get<2>(key) = pair_coord + 9;
				trans_name += '-' + std::to_string(start_coord + 10);
			}

			if (junction_map[0].count(key) == 0) {
				continue;
			}

		}

		bed_out << std::get<0>(key) << "\t" << start_coord << "\t" << end_coord << "\t";
		bed_out << trans_name << "\t" << std::get<0>(junction_data) << "\t" << std::get<2>(junction_data) << "\t";
		bed_out << start_coord << "\t" << end_coord << "\t255,0,0\t2\t" << std::get<3>(junction_data) << "," << std::get<4>(junction_data) << "\t0,0\n";

	}
	bed_out.close();
}

void EnhancedOutput::outputSortedBam()
{
	int header_size = (int)sam_header.size();
	int num_chromosomes = (int)ref_seq_map.size();

	// Connect BAM output to stdout
#ifdef _WIN32
	int result = _setmode(_fileno(stdout), _O_BINARY);
	//	std::setvbuf(stdout, NULL, _IOFBF, 65536);
	bam_stream = bgzf_dopen(_fileno(stdout), "w1");
#else
	bam_stream = bgzf_dopen(fileno(stdout), "w1");
#endif

	// Output header
	bgzf_write(bam_stream, "BAM\001", 4);
	bgzf_write(bam_stream, (char*)&header_size, sizeof(header_size));
	bgzf_write(bam_stream, sam_header.c_str(), header_size);
	bgzf_write(bam_stream, (char*)&num_chromosomes, sizeof(num_chromosomes));
	for (auto const &entry : ref_seq_map) {
		int namelen = (int)(entry.first.size() + 1);
		int seqlen = entry.second[0];
		bgzf_write(bam_stream, (char*)&namelen, sizeof(namelen));
		bgzf_write(bam_stream, entry.first.data(), namelen);
		bgzf_write(bam_stream, (char*)&seqlen, sizeof(seqlen));
	}
	bgzf_flush(bam_stream);

	// Sort and output alignments
	if (num_threads == 1) {

		//Sort chromosomes and output to BAM file
		for (int ref_ID = 0; ref_ID < num_chromosomes; ref_ID++) {
			sortChromosome(ref_ID);
		}

	} else {

		// Sort chromosomes
		std::vector<std::thread> workers;
		for (int thread = 0; thread < num_threads; thread++) {
			workers.emplace_back(&EnhancedOutput::fetchChromosome, this);
		}
		for (int thread = 0; thread < num_threads; thread++) {
			workers[thread].join();
		}

		// Allocate buffer storage for alignments
		char *align_buffer;
		uint64_t max_align_bytes = 0;
		for (int ref_ID = 0; ref_ID < num_chromosomes; ref_ID++) {
			max_align_bytes = (std::max)(max_align_bytes, (uint64_t)sorting_streams[0][ref_ID].tellg());
		}
		align_buffer = new char[max_align_bytes + 1];

		// Output alignments for each chromosome
		for (int ref_ID = 0; ref_ID < num_chromosomes; ref_ID++) {

			// Read alignments for chromosome
			uint64_t stream_bytes = (uint64_t)sorting_streams[0][ref_ID].tellg();
			sorting_streams[0][ref_ID].seekg(std::fstream::beg);
			sorting_streams[0][ref_ID].read(align_buffer, stream_bytes);
			sorting_streams[0][ref_ID].close();
			remove((sort_dir + "/" + std::to_string(ref_ID) + "_0").c_str());

			// Output alignment to BAM file
			uint64_t position = 0;
			for (uint i = 0; i < num_alignments[0][ref_ID]; i++) {
				char *buffer = align_buffer + position;
				uint align_bytes = *((uint*)buffer) + sizeof(uint);
				bgzf_write(bam_stream, buffer, align_bytes);
				position += align_bytes;
			}
		}

		// Deallocate buffer space
		delete[] align_buffer;

	}

	// Flush and close output stream
	bgzf_flush(bam_stream);
	bgzf_close(bam_stream);
}

void EnhancedOutput::fetchChromosome()
{
	int ref_ID;
	int num_chromosomes = (int)ref_seq_map.size();

	while (true) {

		{
			std::lock_guard<std::mutex> lock(sorting_lock);
			ref_ID = ++current_sorting_index;
		}

		if (ref_ID >= num_chromosomes) {
			return;
		}

		sortChromosome(ref_ID);

	}
}

void EnhancedOutput::sortChromosome(int ref_ID)
{
	// Allocate buffer space to store the alignments and sorting data
	char *align_buffer;
	uint64_t *sort_buffer;
	uint64_t total_align_bytes = 0;
	uint num_align = 0;
	for (int thread = 0; thread < num_threads; thread++) {
		total_align_bytes += (uint64_t)sorting_streams[thread][ref_ID].tellg();
		num_align += num_alignments[thread][ref_ID];
	}
	if (num_align == 0) {
		removeSortingFiles(ref_ID);
		return;
	}
	align_buffer = new char[total_align_bytes + 1];
	sort_buffer = new uint64_t[num_align * 2];

	// Load data from file(s)
	uint64_t total_stream_bytes = 0;
	for (int thread = 0; thread < num_threads; thread++) {
		uint64_t stream_bytes = (uint64_t)sorting_streams[thread][ref_ID].tellg();
		sorting_streams[thread][ref_ID].seekg(std::fstream::beg);
		sorting_streams[thread][ref_ID].read(align_buffer + total_stream_bytes, stream_bytes);
		total_stream_bytes += stream_bytes;
	}
	removeSortingFiles(ref_ID);

	// Collect sorting data
	uint64_t position = 0;
	for (uint i = 0; i < num_align; i++) {
		uint *buffer = (uint*)(align_buffer + position);
		sort_buffer[i * 2] = (uint64_t)buffer[2];
		sort_buffer[i * 2 + 1] = position;
		position += buffer[0] + sizeof(uint);
	}

	// Sort by genomic position
	qsort((void*)sort_buffer, num_align, sizeof(uint64_t) * 2, funCompareArrays<uint64_t, 2>);

	// Output sorted alignments
	if (num_threads > 1) {
		sorting_streams[0][ref_ID].seekg(std::fstream::beg);
		num_alignments[0][ref_ID] = num_align;
	}
	for (uint i = 0; i < num_align; i++) {
		char *buffer = align_buffer + sort_buffer[i * 2 + 1];
		uint align_bytes = *((uint*)buffer) + sizeof(uint);
		if (num_threads == 1) {
			bgzf_write(bam_stream, buffer, align_bytes);
		} else {
			sorting_streams[0][ref_ID].write(buffer, align_bytes);
		}
	}

	// Deallocate buffer space
	delete[] align_buffer;
	delete[] sort_buffer;
}

void EnhancedOutput::removeSortingFiles(int ref_ID)
{
	std::string sort_file = sort_dir + "/" + std::to_string(ref_ID) + "_";
	for (int thread = (num_threads == 1) ? 0 : 1; thread < num_threads; thread++) {
		sorting_streams[thread][ref_ID].close();
		remove((sort_file + std::to_string(thread)).c_str());
	}
}

// Calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
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

// Pack nucleotides for BAM
void EnhancedOutput::packseq(const char *seq_in, char *seq_out, uint seq_len)
{
	for (uint i = 0; i < seq_len/2; i++) {
		seq_out[i] = ((encodeNucleotide(seq_in[2 * i])<<4) | encodeNucleotide(seq_in[2 * i + 1]));
	};
	if (seq_len%2 == 1) {
		seq_out[seq_len/2] = (encodeNucleotide(seq_in[seq_len - 1])<<4);
	};
}

// Encode nucleotides for packing
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

DelimRow::DelimRow(char delimiter)
	: delim(delimiter)
{
}

std::string const& DelimRow::operator[](std::size_t index) const
{
	return m_data[index];
}

std::size_t DelimRow::size() const
{
	return m_data.size();
}

void DelimRow::readNextRow(std::istream& str)
{
	std::string line;
	std::getline(str, line);

	std::stringstream lineStream(line);
	std::string cell;

	m_data.clear();
	while (std::getline(lineStream, cell, delim)) {
		m_data.push_back(cell);
	}
	// This checks for a trailing delimiter with no data after it
	if (!lineStream && cell.empty()) {
		// If there was a trailing delimiter then add an empty element
		m_data.push_back("");
	}
}

std::istream& operator>>(std::istream& str, DelimRow& data)
{
	data.readNextRow(str);
	return str;
}