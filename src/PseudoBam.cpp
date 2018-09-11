#include "PseudoBam.h"

/** --- pseudobam functions -- **/

void outputPseudoBam(const KmerIndex &index, const std::vector<int> &u,
    const char *s1, const char *n1, const char *q1, int slen1, int nlen1, const std::vector<std::pair<KmerEntry,int>>& v1,
    const char *s2, const char *n2, const char *q2, int slen2, int nlen2, const std::vector<std::pair<KmerEntry,int>>& v2,
	bool paired, EnhancedOutput &output_handler, int id) {

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

	if (u.empty()) {  // no mapping
		if (output_handler.enhancedoutput) {
			if (!output_handler.outputunmapped) {
				return;
			}
			// HAVE TO ADD OPTION TO OUTPUT UNMAPPED READS FOR ENHANCED_OUTPUT!!!
		} else if (paired) {
			printf("%s\t77\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n1,s1,q1);
			printf("%s\t141\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n2,s2,q2);
	    } else {
			printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n1,s1,q1);
	    }
		return;
	}


	// Trim u to remove repeated genes
	std::vector<int> v;
	if (output_handler.enhancedoutput) {
		std::set<std::string> gene_list;
		for (auto tr : u) {
			if (gene_list.emplace(std::get<0>(output_handler.gene_map.at(index.target_names_[tr]))).second) {
				v.push_back(tr);
			}
		}
	} else {
		v = u;
	}

    if (paired) {  // paired end

		int flag1 = 0x01 + 0x40;
		int flag2 = 0x01 + 0x80;
		int p1 = -1, p2 = -1;
		KmerEntry val1, val2;
		int nmap = v.size();//index.ecmap[ec].size();
		Kmer km1, km2;

		if (v1.empty()) {
			flag1 += 0x04; // read unmapped
			flag2 += 0x08; // mate unmapped
		} else {
			val1 = v1[0].first;
			p1 = v1[0].second;
			for (auto &x : v1) {
				if (x.second < p1) {
					val1 = x.first;
					p1 = x.second;
				}
			}
			km1 = Kmer((s1 + p1));
		}

		if (v2.empty()) {
			flag1 += 0x08; // mate unmapped
			flag2 += 0x04; // read unmapped
		} else {
			val2 = v2[0].first;
			p2 = v2[0].second;
			for (auto &x : v2) {
				if (x.second < p2) {
					val2 = x.first;
					p2 = x.second;
				}
			}
			km2 = Kmer((s2 + p2));
		}

		if (!v1.empty() && !v2.empty()) {
			flag1 += 0x02; // proper pair
			flag2 += 0x02; // proper pair
		}

		bool revset = false;
		// output pseudoalignments for read 1
		for (auto tr : v) {

			std::string ref_name = index.target_names_[tr];
			bool read_sense = true;
			bool mate_sense = true;
			bool gene_sense = true;
			if (output_handler.enhancedoutput) {
				gene_sense = (std::get<2>(output_handler.gene_map.at(ref_name)) >= 0);
			}

			int f1 = flag1;
			std::pair<int, bool> x1 {-1, true};
			std::pair<int, bool> x2 {-1, true};
			if (p1 != -1) {
				x1 = index.findPosition(tr, km1, val1, p1);
				if (p2 == -1) {
					x2 = {x1.first,!x1.second};
				}
				read_sense = x1.second;
				if (read_sense ^ gene_sense) {
					f1 += 0x10; // read reverse
					if (!revset) {
						revseq(&buf1[0], &buf2[0], s1, q1, slen1);
						revset = true;
					}
				}
			}
			if (p2 != -1) {
				x2 = index.findPosition(tr, km2 , val2, p2);
				if (p1 == -1) {
					x1 = {x2.first, !x2.second};
				}
				mate_sense = x2.second;
				if (mate_sense ^ gene_sense) {
					f1 += 0x20; // mate reverse
				}
			}

			int tlen = x2.first - x1.first;
			tlen += sgn<int>(tlen);
			bool mapped = (f1 & 0x04) == 0;
			int flag = f1 & 0xFFFF;
			int posread = (read_sense) ? x1.first : (x1.first - slen1 + 1);
			int posmate = (mate_sense) ? x2.first : (x2.first - slen2 + 1);
			if (v1.empty()) {
				posread = posmate;
			}
			if (v2.empty()) {
				posmate = posread;
			}
			int mapq = (!v1.empty()) ? 255 : 0;
			const char *seq = (f1 & 0x10) ? &buf1[0] : s1;
			const char *qual = (f1 & 0x10) ? &buf2[0] : q1;

			// Temp check for sense!
			//bool fw = (km1 == km1.rep());
			//bool csense = (fw == val1.isFw());
			//bool trsense = (csense == read_sense);
			//std::string gene_name = std::get<0>(output_handler.gene_map.at(ref_name));
			//if (!gene_sense) {
			//	if (trsense) {
			//		if (csense) {
			//			if (output_handler.sense_types[0] == 0) {
			//				std::cerr << "0: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[0]++;
			//			}
			//		} else {
			//			if (output_handler.sense_types[1] == 0) {
			//				std::cerr << "1: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[1]++;
			//			}
			//		}
			//	} else {
			//		if (csense) {
			//			if (output_handler.sense_types[2] < 30) {
			//				std::cerr << "2: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[2]++;
			//			}
			//		} else {
			//			if (output_handler.sense_types[3] == 0) {
			//				std::cerr << "3: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[3]++;
			//			}
			//		}
			//	}
			//} else {
			//	if (trsense) {
			//		if (csense) {
			//			if (output_handler.sense_types[4] == 0) {
			//				std::cerr << "4: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[4]++;
			//			}
			//		} else {
			//			if (output_handler.sense_types[5] == 0) {
			//				std::cerr << "5: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[5]++;
			//			}
			//		}
			//	} else {
			//		if (csense) {
			//			if (output_handler.sense_types[6] == 0) {
			//				std::cerr << "6: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[6]++;
			//			}
			//		} else {
			//			if (output_handler.sense_types[7] == 0) {
			//				std::cerr << "7: " << n1 << " , " << gene_name << std::endl;
			//				output_handler.sense_types[7]++;
			//			}
			//		}
			//	}
			//}

			if (output_handler.enhancedoutput) {  // Convert to genome coordinates

				if (mapped) {  // Add option to output unmapped reads?!
					output_handler.processAlignment(ref_name, flag, posread, slen1, posmate, slen2, tlen, n1, mapq, seq, qual, nmap, id);
				}

			} else {  // Default calculations when not using enhanced output

				getCIGARandSoftClip(cig, bool(f1 & 0x10), mapped, posread, posmate, slen1, index.target_lens_[tr]);

				printf("%s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", n1, flag, ref_name.c_str(), posread, mapq, cig, posmate, tlen, seq, qual, nmap);
				if (v1.empty()) {
					break;  // Only report primary alignment
				}

			}

			flag1 += 0x100;  // for subsequent secondary alignments

		}

		revset = false;
		// output pseudoalignments for read 2
		for (auto tr : v) {

			std::string ref_name = index.target_names_[tr];
			bool read_sense = true;
			bool mate_sense = true;
			bool gene_sense = true;
			if (output_handler.enhancedoutput) {
				gene_sense = (std::get<2>(output_handler.gene_map.at(ref_name)) >= 0);
			}

			int f2 = flag2;
			std::pair<int, bool> x1 {-1,true};
			std::pair<int, bool> x2 {-1,true};
			if (p1 != -1) {
				x1 = index.findPosition(tr, km1, val1, p1);
				if (p2 == -1) {
					x2 = {x1.first,!x1.second};
				}
				read_sense = x1.second;
				if (read_sense ^ gene_sense) {
					f2 += 0x20; // mate reverse
				}
			}
			if (p2 != -1) {
				x2 = index.findPosition(tr, km2, val2, p2);
				if (p1 == -1) {
					x1 = {x2.first, !x2.second};
				}
				mate_sense = x2.second;
				if (mate_sense ^ gene_sense) {
					f2 += 0x10; // read reverse
					if (!revset) {
						revseq(&buf1[0], &buf2[0], s2, q2, slen2);
						revset = true;
					}
				}
			}

			int tlen = x1.first - x2.first;
			tlen += sgn<int>(tlen);
			bool mapped = (f2 & 0x04) == 0;
			int flag = f2 & 0xFFFF;
			int posread = (mate_sense) ? x2.first : (x2.first - slen2 + 1);
			int posmate = (read_sense) ? x1.first : (x1.first - slen1 + 1);
			if (v1.empty()) {
				posmate = posread;
			}
			if (v2.empty()) {
				posread = posmate;
			}
			int mapq = (!v2.empty()) ? 255 : 0;
			const char *seq = (f2 & 0x10) ? &buf1[0] : s2;
			const char *qual = (f2 & 0x10) ? &buf2[0] : q2;

			if (output_handler.enhancedoutput) {  // Convert to genome coordinates

				if (mapped) {  // Add option to output unmapped reads?!
					output_handler.processAlignment(ref_name, flag, posread, slen2, posmate, slen1, tlen, n2, mapq, seq, qual, nmap, id);
				}

			} else {  // Default calculations when not using enhanced output

				getCIGARandSoftClip(cig, bool(f2 & 0x10), mapped, posread, posmate, slen2, index.target_lens_[tr]);

				printf("%s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", n2, flag, ref_name.c_str(), posread, mapq, cig, posmate, tlen, seq, qual, nmap);
				if (v2.empty()) {
					break;  // Only print primary alignment
				}

			}

			flag2 += 0x100;  // for subsequent secondary alignments

		}

    } else {  // single end

		int flag1 = 0;
		int nmap = (int) v.size();
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
		for (auto tr : v) {

			std::string ref_name = index.target_names_[tr];
			bool gene_sense = true;
			if (output_handler.enhancedoutput) {
				gene_sense = (std::get<2>(output_handler.gene_map.at(ref_name)) >= 0);
			}

			int f1 = flag1;
			auto x1 = index.findPosition(tr, km1, val1, p1);
	        if (x1.second ^ gene_sense) {
				f1 += 0x10;
				if (!revset) {
					revseq(&buf1[0], &buf2[0], s1, q1, slen1);
					revset = true;
				}
			}

			bool mapped = (f1 & 0x04) == 0;
			int flag = f1 & 0xFFFF;
			int posread = (x1.second) ? x1.first : (x1.first - slen1 + 1);
			int dummy = 1;
			
			const char *seq = (f1 & 0x10) ? &buf1[0] : s1;
			const char *qual = (f1 & 0x10) ? &buf2[0] : q1;
			if (output_handler.enhancedoutput) {  // Convert to genome coordinates

				if (mapped) {  // Add option to output unmapped reads?!
					output_handler.processAlignment(ref_name, flag, posread, slen1, 0, 0, 0, n1, 255, seq, qual, nmap, id);
				}

			} else {  // Default calculations when not using enhanced output

				getCIGARandSoftClip(cig, bool(f1 & 0x10), mapped, posread, dummy, slen1, index.target_lens_[tr]);

				printf("%s\t%d\t%s\t%d\t255\t%s\t*\t%d\t%d\t%s\t%s\tNH:i:%d\n", n1, flag, ref_name.c_str(), posread, cig, 0, 0, seq, qual, nmap);

			}

			flag1 = 0x100;  // for subsequent secondary alignments

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
