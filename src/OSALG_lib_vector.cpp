#include <string>
#include <algorithm>
#include <vector>
#include <limits>
#include <unordered_map>
#include <immintrin.h> //AVX
#include <stdlib.h>

#include "OSALG_lib.h"

#define N_BORDER 30
#define L 2
#define MATCH_SCORE -2
#define MISMATCH_SCORE 4
#define TRIANGLE_SIZE 8
#define VECTOR_SIZE 16

#define OVERFLOW_CONTROL_VALUE 200

#define STOP -1
#define MATCH 0
#define INSERT 1
#define DELETE 2
#define DELETE1 3
#define DELETE2 4

namespace OSALG_vector {

	std::vector<short int> u{ 4, 2 };
	std::vector<short int> v{ 1, 13 };

	std::unordered_map<int, char> CIGAR_map = {
			{MATCH, 'M'},
			{INSERT, 'I'},
			{DELETE, 'D'}
	};

	short int min_val(short int a, short int b) {
		return ((a < b) ? a : b);
	}

	short int diff(char first, char second) {
		return (first == second) ? MATCH_SCORE : MISMATCH_SCORE;
	}

	int get_last_index(int i, std::string const &seq1, std::string const &seq2) {
		int min = std::min(seq1.length() + 1, seq2.length() + 1);

		if(i < min) {
			return i + 1;
		} else if(i < min + labs(seq2.length() - seq1.length())){
			return min;
		}

		return seq1.length() + seq2.length() + 1 - i;
	}

	void init_first_triangle(short int **d_mat, short int **f1_mat, short int **f2_mat, std::string const &seq1, std::string const &seq2) {
		d_mat[0][0] = d_mat[0][2] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;
		d_mat[0][1] = 0;
		
		f1_mat[0][0] = f1_mat[0][1] = f1_mat[0][2] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;
		f2_mat[0][0] = f2_mat[0][1] = f2_mat[0][2] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;

		for(int i = 1; i < TRIANGLE_SIZE; ++i) {
			int last_ind = get_last_index(i, seq1, seq2);
			
			f1_mat[i][0] = f1_mat[i][last_ind + 1] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;
			f2_mat[i][0] = f2_mat[i][last_ind + 1] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;
			d_mat[i][0] = d_mat[i][last_ind + 1] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;

			for(int j = 1; j <= last_ind; ++j) {
				if(j == 1) {
					//Deletion
					f1_mat[i][j] = d_mat[i][j] = min_val(d_mat[i - 1][j] + v[0], f1_mat[i - 1][j]) + u[0];

					f2_mat[i][j] = min_val(d_mat[i - 1][j] + v[1], f2_mat[i - 1][j]) + u[1];

					d_mat[i][j] = min_val(d_mat[i][j], f2_mat[i][j]);
					
				} else if(j == last_ind) {
					f1_mat[i][j] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;
					f2_mat[i][j] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;

					d_mat[i][j] = d_mat[i - 1][j - 1] + u[0];//INSERTION
				} else {
					//Deletion
					f1_mat[i][j] = d_mat[i][j] = min_val(d_mat[i - 1][j] + v[0], f1_mat[i - 1][j]) + u[0];
					f2_mat[i][j] = min_val(d_mat[i - 1][j] + v[1], f2_mat[i - 1][j]) + u[1];

					d_mat[i][j] = min_val(d_mat[i][j], f2_mat[i][j]);

					//Insertion
					d_mat[i][j] = min_val(d_mat[i][j], d_mat[i - 1][j - 1] + u[0]);
					
					//Match/mismatch
					d_mat[i][j] = min_val(d_mat[i][j], d_mat[i - 2][j - 1] + diff(seq1[last_ind - j - 1], seq2[j - 2]));
				}
			}
		}
	}

	bool test_coord_eligibility(int m, int n, std::string const &seq1, std::string const &seq2) {
		return m > 0 && m <= seq1.length() && n > 0 && n <= seq2.length();
	}

	void updatePosition(unsigned int& i, unsigned int& j, int parent, bool first_half_diags, std::string const &seq1) {
		if(first_half_diags) {
			switch (parent) {
				case MATCH:
					i -= 2;
					j -= 1;
					break;
				case INSERT:
					i -= 1;
					j -= 1;
					break;
				case DELETE:
					i -= 1;
					break;
			}
		} else {
			switch (parent) {
				case MATCH:
					if(i != seq1.length() + 1) {
						j += 1;
					}

					i -= 2;
					break;
				case INSERT:
					i -= 1;
					break;
				case DELETE:
					i -= 1;
					j += 1;
					break;
			}
		}
	}

	void construct_CIGAR(short int **d_mat, short int **f1_mat, short int **f2_mat, int matrix_row_num, std::string &cigar, bool extended_cigar, std::string const &seq1, std::string const &seq2) {
		unsigned int i = matrix_row_num - 1;
		unsigned int j = 1;

		bool firstIdentified = false;
		int counter;
		char lastChar, c;

		bool del_mode = false;

		bool first_half_diags = false;

		while (true) {
			if(i == 0 && j == 1) break;

			if(i <= seq1.length()) first_half_diags = true;

			short parent;
			int last_ind = get_last_index(i, seq1, seq2);		

			if(first_half_diags && j == 1) {
				parent = DELETE;
			} else if(first_half_diags && j == get_last_index(i, seq1, seq2)) {
				parent = INSERT;
			} else {
				if(d_mat[i][j] == f2_mat[i][j]) {
					parent = DELETE;
					del_mode = true;
				} else if(d_mat[i][j] == f1_mat[i][j]) {
					parent = DELETE;
					del_mode = false;
				} else if(d_mat[i][j] == (((first_half_diags) ? d_mat[i - 1][j - 1] : d_mat[i - 1][j]) + u[0])) {
					parent = INSERT;
				} else {
					parent = MATCH;
				}

				if(del_mode) {
					parent = DELETE;
				}	
			}

			if(extended_cigar) {
				c = CIGAR_map.at(parent);

				if(c == 'M') {
					bool match;

					if(first_half_diags) {
						match = (seq1[last_ind - j + ((i > seq2.length()) ? (i - seq2.length()) : 0) - 1] == seq2[j - 2]);
					} else {
						match = (seq1[seq1.length() - j] == seq2[seq2.length() - 1 - last_ind + j - ((i < seq2.length()) ? (seq2.length() - i) : 0)]);
					}

					if(match) {
						c = '=';
					} else {
						c = 'X';
					}
				}
			} else {
				c = CIGAR_map.at(parent);
			}

			if (firstIdentified && c == lastChar) {
				counter++;
			} else {
				if (firstIdentified) {
					if(counter >= 30 && lastChar == 'D') lastChar = 'N';
					cigar = std::to_string(counter) + lastChar + cigar;
				}
				else {
					firstIdentified = true;
				}

				lastChar = c;
				counter = 1;
			}

			updatePosition(i, j, parent, first_half_diags, seq1);
		}

		cigar = std::to_string(counter) + lastChar + cigar;
	}

	int long_gaps_alignment(std::string const &seq1, std::string const &seq2, std::string &cigar, bool extended_cigar) {
		int matrix_row_num = seq1.length() + seq2.length() + 1;

		short int **d_mat = new short int*[matrix_row_num];
		short int **f1_mat = new short int*[matrix_row_num];
		short int **f2_mat = new short int*[matrix_row_num];

		int diagonal_size = ((std::min(seq1.length() + 1, seq2.length() + 1) - 1) / VECTOR_SIZE) * VECTOR_SIZE + VECTOR_SIZE + 2;

		//TO DO ALLOCATION
		for(int i = 0; i < matrix_row_num; ++i) {
			d_mat[i] = new short int[diagonal_size];
			f1_mat[i] = new short int[diagonal_size];
			f2_mat[i] = new short int[diagonal_size];
		}

		init_first_triangle(d_mat, f1_mat, f2_mat, seq1, seq2);
		__m256i insert_penalty_vec = _mm256_set1_epi16(u[0]);

		__m256i u_1_vec = _mm256_set1_epi16(u[0]);
		__m256i u_2_vec = _mm256_set1_epi16(u[1]);
		__m256i v_1_vec = _mm256_set1_epi16(v[0]);
		__m256i v_2_vec = _mm256_set1_epi16(v[1]);

		//first half of diagonals
		for(int i = TRIANGLE_SIZE; i <= seq1.length(); ++i) {
			int last_ind = get_last_index(i, seq1, seq2);

			for(int j = 1; j <= last_ind; j += VECTOR_SIZE) {
				__m256i second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j]);
				//First part of convex function
				__m256i second_diagonal_f1 = _mm256_loadu_si256((__m256i *)&f1_mat[i - 1][j]);

				__m256i third_diagonal_f1 = _mm256_add_epi16(second_diagonal_d, v_1_vec);
				third_diagonal_f1 = _mm256_min_epi16(third_diagonal_f1, second_diagonal_f1);
				third_diagonal_f1 = _mm256_add_epi16(third_diagonal_f1, u_1_vec);
				__m256i third_diagonal_d = third_diagonal_f1;
				//Second part of convex function
				__m256i second_diagonal_f2 = _mm256_loadu_si256((__m256i *)&f2_mat[i - 1][j]);

				__m256i third_diagonal_f2 = _mm256_add_epi16(second_diagonal_d, v_2_vec);
				third_diagonal_f2 = _mm256_min_epi16(third_diagonal_f2, second_diagonal_f2);
				third_diagonal_f2 = _mm256_add_epi16(third_diagonal_f2, u_2_vec);
				third_diagonal_d = _mm256_min_epi16(third_diagonal_f2, third_diagonal_d);

				//Match/mismatch
				int m = last_ind - j + ((i > seq2.length()) ? (i - seq2.length()) : 0);
				int n = j - 1;

				__m256i diff_vec = _mm256_setr_epi16(test_coord_eligibility(m, n, seq1, seq2) ? diff(seq1[m - 1], seq2[n - 1]) : 0, test_coord_eligibility(m - 1, n + 1, seq1, seq2) ? diff(seq1[m - 2], seq2[n]) : 0,
						test_coord_eligibility(m - 2, n + 2, seq1, seq2) ? diff(seq1[m - 3], seq2[n + 1]) : 0, test_coord_eligibility(m - 3, n + 3, seq1, seq2) ? diff(seq1[m - 4], seq2[n + 2]) : 0,
						test_coord_eligibility(m - 4, n + 4, seq1, seq2) ? diff(seq1[m - 5], seq2[n + 3]) : 0, test_coord_eligibility(m - 5, n + 5, seq1, seq2) ? diff(seq1[m - 6], seq2[n + 4]) : 0,
						test_coord_eligibility(m - 6, n + 6, seq1, seq2) ? diff(seq1[m - 7], seq2[n + 5]) : 0, test_coord_eligibility(m - 7, n + 7, seq1, seq2) ? diff(seq1[m - 8], seq2[n + 6]) : 0,
						test_coord_eligibility(m - 8, n + 8, seq1, seq2) ? diff(seq1[m - 9], seq2[n + 7]) : 0, test_coord_eligibility(m - 9, n + 9, seq1, seq2) ? diff(seq1[m - 10], seq2[n + 8]) : 0,
						test_coord_eligibility(m - 10, n + 10, seq1, seq2) ? diff(seq1[m - 11], seq2[n + 9]) : 0, test_coord_eligibility(m - 11, n + 11, seq1, seq2) ? diff(seq1[m - 12], seq2[n + 10]) : 0,
						test_coord_eligibility(m - 12, n + 12, seq1, seq2) ? diff(seq1[m - 13], seq2[n + 11]) : 0, test_coord_eligibility(m - 13, n + 13, seq1, seq2) ? diff(seq1[m - 14], seq2[n + 12]) : 0,
						test_coord_eligibility(m - 14, n + 14, seq1, seq2) ? diff(seq1[m - 15], seq2[n + 13]) : 0, test_coord_eligibility(m - 15, n + 15, seq1, seq2) ? diff(seq1[m - 16], seq2[n + 14]) : 0);

				__m256i first_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 2][j - 1]);

				__m256i third_diagonal_f0 = _mm256_add_epi16(first_diagonal_d, diff_vec);

				third_diagonal_d = _mm256_min_epi16(third_diagonal_f0, third_diagonal_d);

				//Insertion
				second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j - 1]);

				__m256i third_diagonal_ins = _mm256_add_epi16(second_diagonal_d, insert_penalty_vec);
				third_diagonal_d = _mm256_min_epi16(third_diagonal_ins, third_diagonal_d);

				//storing results in matrices
				_mm256_storeu_si256((__m256i *)&d_mat[i][j], third_diagonal_d);
				_mm256_storeu_si256((__m256i *)&f1_mat[i][j], third_diagonal_f1);
				_mm256_storeu_si256((__m256i *)&f2_mat[i][j], third_diagonal_f2);
				
			}
			
			d_mat[i][0] = f1_mat[i][0] = f2_mat[i][0] = d_mat[i][last_ind + 1] = f1_mat[i][last_ind + 1] = f2_mat[i][last_ind + 1] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;

		}

		//second half of diagonals
		for(int i = seq1.length() + 1; i < matrix_row_num; ++i) {
			int last_ind = get_last_index(i, seq1, seq2);

			for(int j = 1; j <= last_ind; j += VECTOR_SIZE) {
				__m256i second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j + 1]);
				//First part of convex function
				__m256i second_diagonal_f1 = _mm256_loadu_si256((__m256i *)&f1_mat[i - 1][j + 1]);

				__m256i third_diagonal_f1 = _mm256_add_epi16(second_diagonal_d, v_1_vec);
				third_diagonal_f1 = _mm256_min_epi16(third_diagonal_f1, second_diagonal_f1);
				third_diagonal_f1 = _mm256_add_epi16(third_diagonal_f1, u_1_vec);
				__m256i third_diagonal_d = third_diagonal_f1;
				//Second part of convex function
				__m256i second_diagonal_f2 = _mm256_loadu_si256((__m256i *)&f2_mat[i - 1][j + 1]);

				__m256i third_diagonal_f2 = _mm256_add_epi16(second_diagonal_d, v_2_vec);
				third_diagonal_f2 = _mm256_min_epi16(third_diagonal_f2, second_diagonal_f2);
				third_diagonal_f2 = _mm256_add_epi16(third_diagonal_f2, u_2_vec);
				third_diagonal_d = _mm256_min_epi16(third_diagonal_f2, third_diagonal_d);

				//Match/mismatch
				int m = seq1.length() + 1 - j;
				int n = seq2.length() - last_ind + j - ((i < seq2.length()) ? (seq2.length() - i) : 0);

				__m256i diff_vec = _mm256_setr_epi16(test_coord_eligibility(m, n, seq1, seq2) ? diff(seq1[m - 1], seq2[n - 1]) : 0, test_coord_eligibility(m - 1, n + 1, seq1, seq2) ? diff(seq1[m - 2], seq2[n]) : 0,
						test_coord_eligibility(m - 2, n + 2, seq1, seq2) ? diff(seq1[m - 3], seq2[n + 1]) : 0, test_coord_eligibility(m - 3, n + 3, seq1, seq2) ? diff(seq1[m - 4], seq2[n + 2]) : 0,
						test_coord_eligibility(m - 4, n + 4, seq1, seq2) ? diff(seq1[m - 5], seq2[n + 3]) : 0, test_coord_eligibility(m - 5, n + 5, seq1, seq2) ? diff(seq1[m - 6], seq2[n + 4]) : 0,
						test_coord_eligibility(m - 6, n + 6, seq1, seq2) ? diff(seq1[m - 7], seq2[n + 5]) : 0, test_coord_eligibility(m - 7, n + 7, seq1, seq2) ? diff(seq1[m - 8], seq2[n + 6]) : 0,
						test_coord_eligibility(m - 8, n + 8, seq1, seq2) ? diff(seq1[m - 9], seq2[n + 7]) : 0, test_coord_eligibility(m - 9, n + 9, seq1, seq2) ? diff(seq1[m - 10], seq2[n + 8]) : 0,
						test_coord_eligibility(m - 10, n + 10, seq1, seq2) ? diff(seq1[m - 11], seq2[n + 9]) : 0, test_coord_eligibility(m - 11, n + 11, seq1, seq2) ? diff(seq1[m - 12], seq2[n + 10]) : 0,
						test_coord_eligibility(m - 12, n + 12, seq1, seq2) ? diff(seq1[m - 13], seq2[n + 11]) : 0, test_coord_eligibility(m - 13, n + 13, seq1, seq2) ? diff(seq1[m - 14], seq2[n + 12]) : 0,
						test_coord_eligibility(m - 14, n + 14, seq1, seq2) ? diff(seq1[m - 15], seq2[n + 13]) : 0, test_coord_eligibility(m - 15, n + 15, seq1, seq2) ? diff(seq1[m - 16], seq2[n + 14]) : 0);

				__m256i first_diagonal_d = _mm256_loadu_si256((__m256i *)((i == seq1.length() + 1) ? &d_mat[i - 2][j] : &d_mat[i - 2][j + 1]));

				__m256i third_diagonal_f0 = _mm256_add_epi16(first_diagonal_d, diff_vec);

				third_diagonal_d = _mm256_min_epi16(third_diagonal_f0, third_diagonal_d);

				//Insertion
				second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j]);

				__m256i third_diagonal_ins = _mm256_add_epi16(second_diagonal_d, insert_penalty_vec);
				third_diagonal_d = _mm256_min_epi16(third_diagonal_ins, third_diagonal_d);

				//storing results in matrices
				_mm256_storeu_si256((__m256i *)&d_mat[i][j], third_diagonal_d);
				_mm256_storeu_si256((__m256i *)&f1_mat[i][j], third_diagonal_f1);
				_mm256_storeu_si256((__m256i *)&f2_mat[i][j], third_diagonal_f2);
			}

			d_mat[i][0] = f1_mat[i][0] = f2_mat[i][0] = d_mat[i][last_ind + 1] = f1_mat[i][last_ind + 1] = f2_mat[i][last_ind + 1] = std::numeric_limits<short int>::max() - OVERFLOW_CONTROL_VALUE;

		}

		construct_CIGAR(d_mat, f1_mat, f2_mat, matrix_row_num, cigar, extended_cigar, seq1, seq2);

		//Deleting allocated memory
		for(int i = 0; i < matrix_row_num; ++i) {
			delete[] d_mat[i];
			delete[] f1_mat[i];
			delete[] f2_mat[i];
		}
		delete[] d_mat;
		delete[] f1_mat;
		delete[] f2_mat;

		return 1;
	}
}
