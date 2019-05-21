#include <string>
#include <algorithm>
#include <vector>
#include <limits>
#include "OSALG_lib.h"

#define N_BORDER 30
#define L 2
#define MATCH_SCORE -2
#define MISMATCH_SCORE 4


namespace OSALG {

	std::vector<int> u{4, 2};
	std::vector<int> v{1, 13};

	struct save2 {
		int p;
		int q;

		save2(int p, int q) : p(p), q(q)
		{}

		save2()
		{}
	} typedef SAVE2_type;//secondary list entries

	struct save1 {
		int m;
		int n;
		SAVE2_type pointer;
	
		save1(int m, int n, SAVE2_type pointer) : m(m), n(n), pointer(pointer)
		{}

		save1(int m, int n, int p, int q) : m(m), n(n), pointer(p, q)
		{}

		save1()
		{}
	} typedef SAVE1_type;//primary list entries

	struct al_inf {
		int d_arr_val = 0;
		int f_arr[L + 2];
		bool e_arr[L + 2];
		SAVE2_type p_arr[L + 2];
	} typedef alignment_info;

	void init_first_row(std::vector<alignment_info> &row) {
		//setting 0,0
		row[0].f_arr[0] = 0;
		row[0].d_arr_val = 0;
		for (int i = 1; i <= L + 1; ++i) {
			row[0].f_arr[i] = std::numeric_limits<int>::max();
		}

		row[0].p_arr[0].p = 0;
		row[0].p_arr[0].q = -1;
		row[0].e_arr[0] = true;

		for (int i = 1; i <= L + 1; ++i) {
			row[0].p_arr[i].p = -1;
			row[0].p_arr[i].q = -1;
		}
		//--------------

		// setting 0, i
		for (int i = 1; i < row.size(); ++i) {
			row[i].f_arr[0] = std::numeric_limits<int>::max();

			row[i].f_arr[L + 1] = std::min(row[i - 1].d_arr_val + v[0], row[i - 1].f_arr[L + 1]) + u[0];

			for (int j = 1; j <= L; ++j) {
				row[i].f_arr[j] = std::numeric_limits<int>::max();
			}

			row[i].d_arr_val = row[i].f_arr[L + 1];

			for (int k = 0; k <= L + 1; ++k) {
				row[i].e_arr[k] = (row[i].d_arr_val == row[i].f_arr[k]);
				row[i].p_arr[k].p = (row[i].e_arr[k]) ? 0 : -1;
				row[i].p_arr[k].q = -1;
			}
		}
	}

	void init_first_column_element(std::vector<alignment_info> &current_row, std::vector<alignment_info> &previous_row) {
		current_row[0].f_arr[0] = std::numeric_limits<int>::max();

		int d_val_temp = std::numeric_limits<int>::max();
		
		current_row[0].f_arr[L + 1] = std::numeric_limits<int>::max();
		for (int j = 1; j <= L; ++j) {
			current_row[0].f_arr[j] = std::min(previous_row[0].d_arr_val + v[j - 1], previous_row[0].f_arr[j]) + u[j - 1];

			if (current_row[0].f_arr[j] < d_val_temp) {
				d_val_temp = current_row[0].f_arr[j];
			}
		}

		current_row[0].d_arr_val = d_val_temp;

		for (int k = 0; k <= L + 1; ++k) {
			current_row[0].e_arr[k] = (current_row[0].d_arr_val == current_row[0].f_arr[k]);
			current_row[0].p_arr[k].p = (current_row[0].e_arr[k]) ? 0 : -1;
			current_row[0].p_arr[k].q = -1;
		}
	}

	void adr_function(SAVE2_type &pp, int m, int n, std::vector<alignment_info> &row, std::vector<SAVE1_type> &primary_list) {
		pp.p = primary_list.size();
		pp.q = -1;

		for (int i = 0; i <= L + 1; ++i) {
			if (row[n].e_arr[i]) {
				primary_list.emplace_back(m, n, row[n].p_arr[i].p, row[n].p_arr[i].q);
			}
		}
	}

	void link_function(SAVE2_type &pp, int m, int n, std::vector<alignment_info> &row, std::vector<SAVE2_type> &secondary_list, int i) {
		
		pp.p = row[n].p_arr[0].p;
		pp.q = secondary_list.size();

		secondary_list.emplace_back(row[n].p_arr[i].p, row[n].p_arr[i].q);
	}

	int diff(char first, char second) {
		return (first == second) ? MATCH_SCORE : MISMATCH_SCORE;
	}

	bool check_for_truth(std::vector<alignment_info> const &row, int n) {
		for (int i = 1; i <= 2 * L; ++i) {
			if (row[n].e_arr[i]) return true;
		}

		return false;
	}

	//(i, j) -> (k, l) -> (m, n)
	void append_CIGAR_segment(int i, int j, int k, int l, int m, int n, std::string &cigar, std::string const &seq1, std::string const &seq2, bool extended_cigar) {
		int temp;
		
		if(k == m) {
			temp = n - l;

			if (temp != 0) {
				cigar = std::to_string(temp) + "I" + cigar;
			}
		} else {
			temp = m - k;
			
			if (temp != 0) {
				cigar = std::to_string(temp) + ((temp < N_BORDER) ? "D" : "N") + cigar;
			}
		}


		temp = k - i;
		if(temp != 0) {
			if(extended_cigar) {
				int i = 0;
				while(i < temp) {
					char c = (seq1[k - i - 1] == seq2[l - i - 1]) ? '=' : 'X';
					int counter = 0;

					while(i < temp) {
						if((c == '=' && seq1[k - i - 1] == seq2[l - i - 1]) || (c == 'X' && seq1[k - i - 1] != seq2[l - i - 1])) {
							++counter;
						} else {
							break;
						}
						++i;
					}

					cigar = std::to_string(counter) + c + cigar;
				}

			} else {
				cigar = std::to_string(temp) + "M" + cigar;
			}
		}
	}

	void construct_CIGAR(std::vector<SAVE1_type> const &primary_list, std::string &cigar, std::string const &seq1, std::string const &seq2, bool extended_cigar) {

		int current_index = primary_list.size() - 1;
		
		while(current_index > 0) {
			int next_index = primary_list[current_index].pointer.p;

			int k, l;
			if (primary_list[next_index].m - primary_list[next_index].n > primary_list[current_index].m - primary_list[current_index].n) {
				k = primary_list[current_index].m;
				l = primary_list[current_index].m - primary_list[next_index].m + primary_list[next_index].n;
			}
			else {
				k = primary_list[current_index].n + primary_list[next_index].m - primary_list[next_index].n;
				l = primary_list[current_index].n;
			}

			append_CIGAR_segment(primary_list[next_index].m, primary_list[next_index].n, k, l, primary_list[current_index].m, primary_list[current_index].n, cigar, seq1, seq2, extended_cigar);
			current_index = next_index;
		}
	}

	int long_gaps_alignment(std::string const &seq1, std::string const &seq2, std::string &cigar, bool extended_cigar) {

		std::vector<SAVE1_type> primary_list;
		primary_list.emplace_back(0, 0, -1, -1);
		std::vector<SAVE2_type> secondary_list;

		short current_row_index = 1;
		short previous_row_index = 0;
		std::vector<std::vector<alignment_info>> row{std::vector<alignment_info>(seq2.length() + 1), std::vector<alignment_info>(seq2.length() + 1)};

		init_first_row(row[previous_row_index]);

		for (unsigned int m = 1; m <= seq1.length(); ++m) {
			init_first_column_element(row[current_row_index], row[previous_row_index]);

			for (unsigned int n = 1; n <= seq2.length(); ++n) {
				//Deletions
				for (int i = 1; i <= L; ++i) {
					row[current_row_index][n].f_arr[i] = std::min(row[previous_row_index][n].d_arr_val + v[i - 1], row[previous_row_index][n].f_arr[i]) + u[i - 1];

					if (row[previous_row_index][n].d_arr_val + v[i - 1] < row[previous_row_index][n].f_arr[i]) {
						row[current_row_index][n].p_arr[i] = row[previous_row_index][n].p_arr[0];
					}
					else if (row[previous_row_index][n].d_arr_val + v[i - 1] == row[previous_row_index][n].f_arr[i]) {
						link_function(row[current_row_index][n].p_arr[i], m - 1, n, row[previous_row_index], secondary_list, i);
					}
					else {
						row[current_row_index][n].p_arr[i] = row[previous_row_index][n].p_arr[i];
					}
				}

				//Insertions
				row[current_row_index][n].f_arr[L + 1] = std::min(row[current_row_index][n - 1].d_arr_val + v[0], row[current_row_index][n - 1].f_arr[L + 1]) + u[0];

				if (row[current_row_index][n - 1].d_arr_val + v[0] < row[current_row_index][n - 1].f_arr[L + 1]) {
					row[current_row_index][n].p_arr[L + 1] = row[current_row_index][n - 1].p_arr[0];
				} 
				else if (row[current_row_index][n - 1].d_arr_val + v[0] == row[current_row_index][n - 1].f_arr[L + 1]) {
					link_function(row[current_row_index][n].p_arr[L + 1], m, n - 1, row[current_row_index], secondary_list, L + 1);
				}
				else {
					row[current_row_index][n].p_arr[L + 1] = row[current_row_index][n - 1].p_arr[L + 1];
				}

				//Matches/mismatches
				row[current_row_index][n].f_arr[0] = row[previous_row_index][n - 1].d_arr_val + diff(seq1[m - 1], seq2[n - 1]);

				row[current_row_index][n].d_arr_val = *std::min_element(row[current_row_index][n].f_arr, row[current_row_index][n].f_arr + L + 2);

				for (int i = 0; i <= L + 1; ++i) {
					row[current_row_index][n].e_arr[i] = (row[current_row_index][n].d_arr_val == row[current_row_index][n].f_arr[i]);
				}

				if (row[current_row_index][n].e_arr[0] && check_for_truth(row[previous_row_index], n - 1) && !row[previous_row_index][n - 1].e_arr[0]) {
					adr_function(row[current_row_index][n].p_arr[0], m - 1, n - 1, row[previous_row_index], primary_list);
				}
				else {
					row[current_row_index][n].p_arr[0] = row[previous_row_index][n - 1].p_arr[0];
				}

			}
			std::swap(current_row_index, previous_row_index);
			//previous_row = current_row;
		}

		SAVE2_type p_last;
		adr_function(p_last, seq1.length(), seq2.length(), row[current_row_index], primary_list);

		construct_CIGAR(primary_list, cigar, seq1, seq2, extended_cigar);

		return 0;
	}
}
