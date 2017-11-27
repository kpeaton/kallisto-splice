#include "PseudoBam.h"

/** --- pseudobam functions -- **/

void outputPseudoBam(const KmerIndex &index, const std::vector<int> &u,
    const char *s1, const char *n1, const char *q1, int slen1, int nlen1, const std::vector<std::pair<KmerEntry,int>>& v1,
    const char *s2, const char *n2, const char *q2, int slen2, int nlen2, const std::vector<std::pair<KmerEntry,int>>& v2,
	bool paired, EnhancedOutput &output_handler) {

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
  } else {

	if (output_handler.enhancedoutput) {
		// Trim u to remove repeated genes!
	}

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
	  output_handler.gene_list.clear();
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

		bool mapped = (f1 & 0x04) == 0;
		int flag = f1 & 0xFFFF;
		int posread = (f1 & 0x10) ? (x1.first - slen1 + 1) : x1.first;
		int posmate = (f1 & 0x20) ? (x2.first - slen2 + 1) : x2.first;
		std::string ref_name = index.target_names_[tr];
		const char *seq = (f1 & 0x10) ? &buf1[0] : s1;
		const char *qual = (f1 & 0x10) ? &buf2[0] : q1;
		HighResTimer timer;
		if (!output_handler.enhancedoutput) {  // Default calculations when not using enhanced output

			getCIGARandSoftClip(cig, bool(f1 & 0x10), mapped, posread, posmate, slen1, index.target_lens_[tr]);

			printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", n1, flag, ref_name.c_str(), posread, cig, posmate, tlen, seq, qual, nmap);

		} else {  // Convert to genome coordinates and build CIGAR string

			int strand = 0;

			if (output_handler.getSamData(ref_name, cig, strand, mapped, posread, posmate, slen1, slen2)) {
				output_handler.output_time += timer.timeSinceReset();
				continue;
			}

			if (output_handler.sortedbam) {

				output_handler.outputBamAlignment(ref_name, posread, flag, slen1, posmate, tlen, n1, seq, qual, nmap, strand);

			} else {

				printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d", n1, flag, ref_name.c_str(), posread, cig, posmate, tlen, seq, qual, nmap);
				if (strand == 0) {
					printf("\n");
				} else if (strand < 0) {
					printf("\tXS:A:-\n");
				} else {
					printf("\tXS:A:+\n");
				}

			}

		}
		output_handler.output_time += timer.timeSinceReset();

      }

      revset = false;
      // output pseudoalignments for read 2
      firstTr = true;
	  output_handler.gene_list.clear();
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

		bool mapped = (f2 & 0x04) == 0;
		int flag = f2 & 0xFFFF;
		int posread = (f2 & 0x10) ? (x2.first - slen2 + 1) : x2.first;
		int posmate = (f2 & 0x20) ? (x1.first - slen1 + 1) : x1.first;
		std::string ref_name = index.target_names_[tr];
		const char *seq = (f2 & 0x10) ? &buf1[0] : s2;
		const char *qual = (f2 & 0x10) ? &buf2[0] : q2;
		HighResTimer timer;
		if (!output_handler.enhancedoutput) {  // Default calculations when not using enhanced output

			getCIGARandSoftClip(cig, bool(f2 & 0x10), mapped, posread, posmate, slen2, index.target_lens_[tr]);

			printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", n2, flag, ref_name.c_str(), posread, cig, posmate, tlen, seq, qual, nmap);

		} else {  // Convert to genome coordinates and build CIGAR string

			int strand = 0;

			if (output_handler.getSamData(ref_name, cig, strand, mapped, posread, posmate, slen2, slen1)) {
				continue;
			}

			if (output_handler.sortedbam) {

				output_handler.outputBamAlignment(ref_name, posread, flag, slen2, posmate, tlen, n2, seq, qual, nmap, strand);

			} else {

				printf("%s\t%d\t%s\t%d\t255\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d", n2, flag, ref_name.c_str(), posread, cig, posmate, tlen, seq, qual, nmap);
				if (strand == 0) {
					printf("\n");
				} else if (strand < 0) {
					printf("\tXS:A:-\n");
				} else {
					printf("\tXS:A:+\n");
				}

			}

		}
		output_handler.output_time += timer.timeSinceReset();

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
	  output_handler.gene_list.clear();
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

		bool mapped = (f1 & 0x04) == 0;
		int flag = f1 & 0xFFFF;
		int posread = (f1 & 0x10) ? (x1.first - slen1 + 1) : x1.first;
		int dummy = 1;
		std::string ref_name = index.target_names_[tr];
		const char *seq = (f1 & 0x10) ? &buf1[0] : s1;
		const char *qual = (f1 & 0x10) ? &buf2[0] : q1;
		if (!output_handler.enhancedoutput) {  // Default calculations when not using enhanced output

			getCIGARandSoftClip(cig, bool(f1 & 0x10), mapped, posread, dummy, slen1, index.target_lens_[tr]);

			printf("%s\t%d\t%s\t%d\t255\t%s\t*\t%d\t%d\t%s\t%s\tNH:i:%d\n", n1, flag, ref_name.c_str(), posread, cig, 0, 0, seq, qual, nmap);

		} else {  // Convert to genome coordinates and build CIGAR string

			int strand = 0;

			if (output_handler.getSamData(ref_name, cig, strand, mapped, posread, dummy, slen1, 0)) {
				continue;
			}

			if (output_handler.sortedbam) {

				output_handler.outputBamAlignment(ref_name, posread, flag, slen1, dummy, 0, n1, seq, qual, nmap, strand);

			} else {

				printf("%s\t%d\t%s\t%d\t255\t%s\t*\t%d\t%d\t%s\t%s\tNH:i:%d", n1, flag, ref_name.c_str(), posread, cig, 0, 0, seq, qual, nmap);
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
