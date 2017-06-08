#include "PseudoBam.h"

/** --- pseudobam functions -- **/

void outputPseudoBam(const KmerIndex &index, const std::vector<int> &u,
    const char *s1, const char *n1, const char *q1, int slen1, int nlen1, const std::vector<std::pair<KmerEntry,int>>& v1,
    const char *s2, const char *n2, const char *q2, int slen2, int nlen2, const std::vector<std::pair<KmerEntry,int>>& v2,
	bool paired, const ExonMap& exmap) {

  static char buf1[32768];
  static char buf2[32768];
  static char cig_[1000];
  char *cig = &cig_[0];


  if (nlen1 > 2 && n1[nlen1-2] == '/') {
    ((char*)n1)[nlen1-2] = 0;
  }

  if (paired && nlen2 > 2 && n2[nlen2-2] == '/') {
    ((char*)n2)[nlen2-2] = 0;
  }

  if (u.empty()) {
    // no mapping
    if (paired) {
      printf("%s\t77\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n1,s1,q1);
      printf("%s\t141\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n2,s2,q2);
      //o << seq1->name.s << "" << seq1->seq.s << "\t" << seq1->qual.s << "\n";
      //o << seq2->name.s << "\t141\t*\t0\t0\t*\t*\t0\t0\t" << seq2->seq.s << "\t" << seq2->qual.s << "\n";
    } else {
      printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n1,s1,q1);
    }
  } else {
    if (paired) {

      int flag1 = 0x01 + 0x40;
      int flag2 = 0x01 + 0x80;

      if (v1.empty()) {
        flag1 += 0x04; // read unmapped
        flag2 += 0x08; // mate unmapped
      }

      if (v2.empty()) {
        flag1 += 0x08; // mate unmapped
        flag2 += 0x04; // read unmapped
      }

      if (!v1.empty() && !v2.empty()) {
        flag1 += 0x02; // proper pair
        flag2 += 0x02; // proper pair
      }


      int p1 = -1, p2 = -1;
      KmerEntry val1, val2;
      int nmap = u.size();//index.ecmap[ec].size();
      Kmer km1, km2;

      if (!v1.empty()) {
        val1 = v1[0].first;
        p1 = v1[0].second;
        for (auto &x : v1) {
          if (x.second < p1) {
            val1 = x.first;
            p1 = x.second;
          }
        }
        km1 = Kmer((s1+p1));
      }

      if (!v2.empty()) {
        val2 = v2[0].first;
        p2 = v2[0].second;
        for (auto &x : v2) {
          if (x.second < p2) {
            val2 = x.first;
            p2 = x.second;
          }
        }
        km2 = Kmer((s2+p2));
      }

      bool revset = false;

      // output pseudoalignments for read 1
      bool firstTr = true;
	  std::set<std::string> gene_list;
      for (auto tr : u) {
        int f1 = flag1;
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = index.findPosition(tr, km1, val1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
          if (!x1.second) {
            f1 += 0x10; // read reverse
            if (!revset) {
              revseq(&buf1[0], &buf2[0], s1, q1, slen1);
              revset = true;
            }
          }
        }
        if (!firstTr) {
          f1 += 0x100; // secondary alignment
        }
        if (p2 != -1) {
          x2 = index.findPosition(tr, km2 , val2, p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
          if (!x2.second) {
            f1 += 0x20; // mate reverse
          }
        }
        firstTr = false;

        int tlen = x2.first - x1.first;
        if (tlen != 0) {
          tlen += (tlen>0) ? 1 : -1;
        }

		int posread = (f1 & 0x10) ? (x1.first - slen1 + 1) : x1.first;
		int posmate = (f1 & 0x20) ? (x2.first - slen2 + 1) : x2.first;
		int strand = 0;
		const char * trans_name = index.target_names_[tr].c_str();
		if (exmap.empty()) {  // Default calculations when exon coordinates map is empty

			getCIGARandSoftClip(cig, bool(f1 & 0x10), (f1 & 0x04) == 0, posread, posmate, slen1, index.target_lens_[tr]);

		} else {  // Convert to genome coordinates and build CIGAR string

			ExonMap::const_iterator map_entry = exmap.find(trans_name);
			if (map_entry == exmap.end()) {
				std::cerr << "Transcript name could not be found in exon coordinate file: " << trans_name << std::endl;
				exit(1);
			}

			const char * gene_name = std::get<0>(map_entry->second).c_str();
			if (gene_list.find(gene_name) == gene_list.end()) {
				gene_list.emplace(gene_name);
			} else {
				continue;
			}

			trans_name = std::get<1>(map_entry->second).c_str();
			strand = (std::get<2>(map_entry->second)[0][2] < 0) ? -1 : 1;  // Apply elsewhere!
			
			std::string cig_string;
			int read_len = slen1;
			int mate_len = slen2;
			int read_offset = 0;
			int mate_offset = 0;
			for (auto &span : std::get<2>(map_entry->second)) {

				if (read_len != 0) {  // Not done mapping coordinates for read
					if (read_len != slen1) {  // In the process of mapping coordinates for read

						if (span[2] < 0) {  // Negative strand               Use trsense from findPosition instead??!!

							sprintf(cig, "%dN", span[2] - span[1] + span[0] - 1 - read_offset);
							cig_string.insert(0, cig);
							read_offset = span[1] - span[0] + 1;
							if (read_len > read_offset) {
								sprintf(cig, "%dM", read_offset);
								cig_string.insert(0, cig);
								read_len -= read_offset;
								read_offset = span[2];
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string.insert(0, cig);
								posread = read_offset - span[2] - read_len;
								read_len = 0;
							}

						} else {  // Positive strand

							sprintf(cig, "%dN", span[2] - read_offset - 1);
							cig_string += cig;
							read_offset = span[1] - span[0] + 1;
							if (read_len > read_offset){
								sprintf(cig, "%dM", read_offset);
								cig_string += cig;
								read_len -= read_offset;
								read_offset += span[2] - 1;
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string += cig;
								read_len = 0;
							}

						}

					} else if (posread <= span[1]) {  // Begin mapping coordinates for read
						
						if (posread < span[0]) {
							read_offset = span[0] - posread;
							sprintf(cig, "%dS", read_offset);
							cig_string = cig;
							read_len -= read_offset;
						}

						if (span[2] < 0) {  // Negative strand

							if ((posread + slen1 - 1) > span[1]) {
								read_offset = span[1] - posread - read_offset + 1;
								sprintf(cig, "%dM", read_offset);
								cig_string.insert(0, cig);
								read_len -= read_offset;
								read_offset = span[2];
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string.insert(0, cig);
								read_len = 0;
								posread = span[1] - span[2] - posread - slen1 + 1;
							}

						} else {  // Positive strand

							if ((posread + slen1 - 1) > span[1]) {
								read_offset = span[1] - posread - read_offset + 1;
								sprintf(cig, "%dM", read_offset);
								cig_string += cig;
								read_len -= read_offset;
								read_offset = span[2] + span[1] - span[0];
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string += cig;
								read_len = 0;
							}
							posread += span[2] - span[0];

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
					sprintf(cig, "%dS", read_len);
					cig_string.insert(0, cig);
					posread = -read_offset - read_len;
				} else {  // Positive strand
					sprintf(cig, "%dS", read_len);
					cig_string += cig;
				}
			}

			if ((mate_len != 0) && (mate_offset < 0)) {  // Account for mate overhangs on negative strand
				posmate = -mate_offset - mate_len;
			}

			if ((f1 & 0x04) == 0) {
				sprintf(cig, "%s", cig_string.c_str());
			} else {
				sprintf(cig, "*");
			}

		}

        printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d", n1, f1 & 0xFFFF, trans_name, posread, cig, posmate, tlen, (f1 & 0x10) ? &buf1[0] : s1, (f1 & 0x10) ? &buf2[0] : q1, nmap);
		if (strand == 0) {
			printf("\n");
		} else if (strand < 0) {
			printf("\tXS:A:-\n");
		} else {
			printf("\tXS:A:+\n");
		}

      }

      revset = false;
      // output pseudoalignments for read 2
      firstTr = true;
	  gene_list.clear();
      for (auto tr : u) {
        int f2 = flag2;
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = index.findPosition(tr, km1, val1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
          if (!x1.second) {
            f2 += 0x20; // mate reverse
          }
        }
        if (p2 != -1) {
          x2 = index.findPosition(tr, km2, val2, p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
          if (!x2.second) {
            f2 += 0x10; // read reverse
            if (!revset) {
              revseq(&buf1[0], &buf2[0], s2, q2, slen2);
              revset = true;
            }

          }
        }
        if (!firstTr) {
          f2 += 0x100; // secondary alignment
        }

        firstTr = false;
        
        int tlen = x1.first - x2.first;
        if (tlen != 0) {
          tlen += (tlen > 0) ? 1 : -1;
        }

		int posread = (f2 & 0x10) ? (x2.first - slen2 + 1) : x2.first;
		int posmate = (f2 & 0x20) ? (x1.first - slen1 + 1) : x1.first;
		int strand = 0;
		const char * trans_name = index.target_names_[tr].c_str();
		if (exmap.empty()) {  // Default calculations when exon coordinates map is empty

			getCIGARandSoftClip(cig, bool(f2 & 0x10), (f2 & 0x04) == 0, posread, posmate, slen2, index.target_lens_[tr]);

		} else {  // Convert to genome coordinates and build CIGAR string

			ExonMap::const_iterator map_entry = exmap.find(trans_name);
			if (map_entry == exmap.end()) {
				std::cerr << "Transcript name could not be found in exon coordinate file: " << trans_name << std::endl;
				exit(1);
			}

			const char * gene_name = std::get<0>(map_entry->second).c_str();
			if (gene_list.find(gene_name) == gene_list.end()) {
				gene_list.emplace(gene_name);
			} else {
				continue;
			}

			trans_name = std::get<1>(map_entry->second).c_str();
			strand = (std::get<2>(map_entry->second)[0][2] < 0) ? -1 : 1;

			std::string cig_string;
			int read_len = slen2;
			int mate_len = slen1;
			int read_offset = 0;
			int mate_offset = 0;
			for (auto &span : std::get<2>(map_entry->second)) {

				if (read_len != 0) {  // Not done mapping coordinates for read
					if (read_len != slen2) {  // In the process of mapping coordinates for read

						if (span[2] < 0) {  // Negative strand

							sprintf(cig, "%dN", span[2] - span[1] + span[0] - 1 - read_offset);
							cig_string.insert(0, cig);
							read_offset = span[1] - span[0] + 1;
							if (read_len > read_offset) {
								sprintf(cig, "%dM", read_offset);
								cig_string.insert(0, cig);
								read_len -= read_offset;
								read_offset = span[2];
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string.insert(0, cig);
								posread = read_offset - span[2] - read_len;
								read_len = 0;
							}

						} else {  // Positive strand

							sprintf(cig, "%dN", span[2] - read_offset - 1);
							cig_string += cig;
							read_offset = span[1] - span[0] + 1;
							if (read_len > read_offset){
								sprintf(cig, "%dM", read_offset);
								cig_string += cig;
								read_len -= read_offset;
								read_offset += span[2] - 1;
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string += cig;
								read_len = 0;
							}

						}

					} else if (posread <= span[1]) {  // Begin mapping coordinates for read

						if (posread < span[0]) {
							read_offset = span[0] - posread;
							sprintf(cig, "%dS", read_offset);
							cig_string = cig;
							read_len -= read_offset;
						}

						if (span[2] < 0) {  // Negative strand

							if ((posread + slen2 - 1) > span[1]) {
								read_offset = span[1] - posread - read_offset + 1;
								sprintf(cig, "%dM", read_offset);
								cig_string.insert(0, cig);
								read_len -= read_offset;
								read_offset = span[2];
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string.insert(0, cig);
								read_len = 0;
								posread = span[1] - span[2] - posread - slen2 + 1;
							}

						} else {  // Positive strand

							if ((posread + slen2 - 1) > span[1]) {
								read_offset = span[1] - posread - read_offset + 1;
								sprintf(cig, "%dM", read_offset);
								cig_string += cig;
								read_len -= read_offset;
								read_offset = span[2] + span[1] - span[0];
							} else {
								sprintf(cig, "%dM", read_len);
								cig_string += cig;
								read_len = 0;
							}
							posread += span[2] - span[0];

						}
					}
				}

				if (mate_len != 0) {  // Not done mapping coordinates for mate
					if (mate_len != slen1) {  // In the process of mapping coordinates for mate

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

							if ((posmate + slen1 - 1) > span[1]) {
								mate_len -= span[1] - posmate - mate_offset + 1;
								mate_offset = span[2];
							} else {
								mate_len = 0;
								posmate = span[1] - span[2] - posmate - slen1 + 1;
							}

						} else {  // Positive strand

							if ((posmate + slen1 - 1) > span[1]) {
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
					sprintf(cig, "%dS", read_len);
					cig_string.insert(0, cig);
					posread = -read_offset - read_len;
				} else {  // Positive strand
					sprintf(cig, "%dS", read_len);
					cig_string += cig;
				}
			}

			if ((mate_len != 0) && (mate_offset < 0)) {  // Account for mate overhangs on negative strand
				posmate = -mate_offset - mate_len;
			}

			if ((f2 & 0x04) == 0) {
				sprintf(cig, "%s", cig_string.c_str());
			} else {
				sprintf(cig, "*");
			}

		}

        printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d", n2, f2 & 0xFFFF, trans_name, posread, cig, posmate, tlen, (f2 & 0x10) ? &buf1[0] : s2,  (f2 & 0x10) ? &buf2[0] : q2, nmap);
		if (strand == 0) {
			printf("\n");
		} else if (strand < 0) {
			printf("\tXS:A:-\n");
		} else {
			printf("\tXS:A:+\n");
		}

      }


    } else {
      // single end
      int nmap = (int) u.size();
      KmerEntry val1 = v1[0].first;
      int p1 = v1[0].second;
      for (auto &x : v1) {
        if (x.second < p1) {
          val1 = x.first;
          p1 = x.second;
        }
      }
      Kmer km1 = Kmer((s1+p1));

      bool revset = false;
      bool firstTr = true;
	  std::set<std::string> gene_list;
      for (auto tr : u) {
        int f1 = 0;
        auto x1 = index.findPosition(tr, km1, val1, p1);

        if (!x1.second) {
          f1 += 0x10;
          if (!revset) {
            revseq(&buf1[0], &buf2[0], s1, q1, slen1);
            revset = true;
          }
        }
        if (!firstTr) {
          f1 += 0x100; // secondary alignment
        }
        firstTr = false;

		int posread = (f1 & 0x10) ? (x1.first - slen1 + 1) : x1.first;
		int strand = 0;
		const char * trans_name = index.target_names_[tr].c_str();
		if (exmap.empty()) {  // Default calculations when exon coordinates map is empty

			int dummy = 1;
			getCIGARandSoftClip(cig, bool(f1 & 0x10), (f1 & 0x04) == 0, posread, dummy, slen1, index.target_lens_[tr]);

		} else {  // Convert to genome coordinates and build CIGAR string

			ExonMap::const_iterator map_entry = exmap.find(trans_name);
			if (map_entry == exmap.end()) {
				std::cerr << "Transcript name could not be found in exon coordinate file: " << trans_name << std::endl;
				exit(1);
			}

			const char * gene_name = std::get<0>(map_entry->second).c_str();
			if (gene_list.find(gene_name) == gene_list.end()) {
				gene_list.emplace(gene_name);
			} else {
				continue;
			}

			trans_name = std::get<1>(map_entry->second).c_str();
			strand = (std::get<2>(map_entry->second)[0][2] < 0) ? -1 : 1;

			std::string cig_string;
			int read_len = slen1;
			int read_offset = 0;
			for (auto &span : std::get<2>(map_entry->second)) {

				if (read_len != slen1) {  // In the process of mapping coordinates for read

					if (span[2] < 0) {  // Negative strand

						sprintf(cig, "%dN", span[2] - span[1] + span[0] - 1 - read_offset);
						cig_string.insert(0, cig);
						read_offset = span[1] - span[0] + 1;
						if (read_len > read_offset) {
							sprintf(cig, "%dM", read_offset);
							cig_string.insert(0, cig);
							read_len -= read_offset;
							read_offset = span[2];
						} else {
							sprintf(cig, "%dM", read_len);
							cig_string.insert(0, cig);
							posread = read_offset - span[2] - read_len;
							read_len = 0;
							break;
						}

					} else {  // Positive strand

						sprintf(cig, "%dN", span[2] - read_offset - 1);
						cig_string += cig;
						read_offset = span[1] - span[0] + 1;
						if (read_len > read_offset){
							sprintf(cig, "%dM", read_offset);
							cig_string += cig;
							read_len -= read_offset;
							read_offset += span[2] - 1;
						} else {
							sprintf(cig, "%dM", read_len);
							cig_string += cig;
							read_len = 0;
							break;
						}

					}

				} else if (posread <= span[1]) {  // Begin mapping coordinates for read

					if (posread < span[0]) {
						read_offset = span[0] - posread;
						sprintf(cig, "%dS", read_offset);
						cig_string = cig;
						read_len -= read_offset;
					}

					if (span[2] < 0) {  // Negative strand

						if ((posread + slen1 - 1) > span[1]) {
							read_offset = span[1] - posread - read_offset + 1;
							sprintf(cig, "%dM", read_offset);
							cig_string.insert(0, cig);
							read_len -= read_offset;
							read_offset = span[2];
						} else {
							sprintf(cig, "%dM", read_len);
							cig_string.insert(0, cig);
							posread = span[1] - span[2] - posread - slen1 + 1;
							read_len = 0;
							break;
						}

					} else {  // Positive strand

						if ((posread + slen1 - 1) > span[1]) {
							read_offset = span[1] - posread - read_offset + 1;
							sprintf(cig, "%dM", read_offset);
							cig_string += cig;
							read_len -= read_offset;
							read_offset = span[2] + span[1] - span[0];
							posread += span[2] - span[0];
						} else {
							sprintf(cig, "%dM", read_len);
							cig_string += cig;
							posread += span[2] - span[0];
							read_len = 0;
							break;
						}

					}

				}

			}

			if (read_len != 0) {  // Account for read overhangs
				if (read_offset < 0) {  // Negative strand
					sprintf(cig, "%dS", read_len);
					cig_string.insert(0, cig);
					posread = -read_offset - read_len;
				} else {  // Positive strand
					sprintf(cig, "%dS", read_len);
					cig_string += cig;
				}
			}

			if ((f1 & 0x04) == 0) {
				sprintf(cig, "%s", cig_string.c_str());
			} else {
				sprintf(cig, "*");
			}

		}

        printf("%s\t%d\t%s\t%d\t255\t%s\t*\t%d\t%d\t%s\t%s\tNH:i:%d", n1, f1 & 0xFFFF, trans_name, posread, cig, 0, 0, (f1 & 0x10) ? &buf1[0] : s1, (f1 & 0x10) ? &buf2[0] : q1, nmap);
		if (strand == 0) {
			printf("\n");
		} else if (strand < 0) {
			printf("\tXS:A:-\n");
		} else {
			printf("\tXS:A:+\n");
		}

      }
    }
  }
}


void revseq(char *b1, char *b2, const char *s, const char *q, int n) {
  b1[n] = 0;
  for (int i = 0; i < n; i++) {
    switch(s[i]) {
    case 'A': b1[n-1-i] = 'T'; break;
    case 'C': b1[n-1-i] = 'G'; break;
    case 'G': b1[n-1-i] = 'C'; break;
    case 'T': b1[n-1-i] = 'A'; break;
    default:  b1[n-1-i] = 'N';
    }
  }
  b2[n] = 0;
  for (int i = 0; i < n; i++) {
    b2[n-1-i] = q[i];
  }
}



void getCIGARandSoftClip(char* cig, bool strand, bool mapped, int &posread, int &posmate, int length, int targetlength) {
  int softclip = 1 - posread;
  int overhang = (posread + length) - targetlength - 1;

  if (posread <= 0) {
    posread = 1;
  }

  if (mapped) {
    if (softclip > 0) {
      if (overhang > 0) {
        sprintf(cig, "%dS%dM%dS",softclip, (length-overhang - softclip), overhang);
      } else {
        sprintf(cig, "%dS%dM",softclip,length-softclip);
      }
    } else if (overhang > 0) {
      sprintf(cig, "%dM%dS", length-overhang, overhang);
    } else {
      sprintf(cig, "%dM",length);
    }
  } else {
    sprintf(cig, "*");
  }


  if (posmate <= 0) {
    posmate = 1;
  }
}
