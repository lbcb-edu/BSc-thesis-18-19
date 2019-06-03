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
#define VECTOR_SIZE 8

#define OVERFLOW_CONTROL_VALUE 100

#define STOP -1
#define MATCH 0
#define INSERT 1
#define DELETE 2
#define DELETE1 3
#define DELETE2 4

namespace OSALG_vector {

	std::vector<int> u{ 4, 2 };
	std::vector<int> v{ 1, 13 };

	int diff(char first, char second) {
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

	void init_first_triangle(int **d_mat, int **f1_mat, int **f2_mat, std::string const &seq1, std::string const &seq2) {
		d_mat[0][0] = d_mat[0][2] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
		d_mat[0][1] = 0;
		
		f1_mat[0][0] = f1_mat[0][1] = f1_mat[0][2] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
		f2_mat[0][0] = f2_mat[0][1] = f2_mat[0][2] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;

		for(int i = 1; i < TRIANGLE_SIZE; ++i) {
			int last_ind = get_last_index(i, seq1, seq2);
			
			f1_mat[i][0] = f1_mat[i][last_ind + 1] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
			f2_mat[i][0] = f2_mat[i][last_ind + 1] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
			d_mat[i][0] = d_mat[i][last_ind + 1] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;

			for(int j = 1; j <= last_ind; ++j) {
				if(j == 1) {
					//Deletion
					f1_mat[i][j] = d_mat[i][j] = std::min(d_mat[i - 1][j] + v[0], f1_mat[i - 1][j]) + u[0];

					f2_mat[i][j] = std::min(d_mat[i - 1][j] + v[1], f2_mat[i - 1][j]) + u[1];

					d_mat[i][j] = std::min(d_mat[i][j], f2_mat[i][j]);
					
				} else if(j == last_ind) {
					f1_mat[i][j] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
					f2_mat[i][j] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;

					d_mat[i][j] = d_mat[i - 1][j - 1] + u[0];//INSERTION
				} else {
					//Deletion
					f1_mat[i][j] = d_mat[i][j] = std::min(d_mat[i - 1][j] + v[0], f1_mat[i - 1][j]) + u[0];
					f2_mat[i][j] = std::min(d_mat[i - 1][j] + v[1], f2_mat[i - 1][j]) + u[1];

					d_mat[i][j] = std::min(d_mat[i][j], f2_mat[i][j]);

					//Insertion
					d_mat[i][j] = std::min(d_mat[i][j], d_mat[i - 1][j - 1] + u[0]);
					
					//Match/mismatch
					d_mat[i][j] = std::min(d_mat[i][j], d_mat[i - 2][j - 1] + diff(seq1[last_ind - j - 1], seq2[j - 2]));
				}
			}
		}
	}

	bool test_coord_eligibility(int m, int n, std::string const &seq1, std::string const &seq2) {
		return m > 0 && m <= seq1.length() && n > 0 && n <= seq2.length();
	}

	int long_gaps_alignment(std::string const &seq1, std::string const &seq2, std::string &cigar, bool extended_cigar) {
		int matrix_row_num = seq1.length() + seq2.length() + 1;

		int **d_mat = new int*[matrix_row_num];
		int **f1_mat = new int*[matrix_row_num];
		int **f2_mat = new int*[matrix_row_num];

		int diagonal_size = ((std::min(seq1.length() + 1, seq2.length() + 1) - 1) / 8) * 8 + 8 + 2;

		//TO DO ALLOCATION
		for(int i = 0; i < matrix_row_num; ++i) {
			d_mat[i] = new int[diagonal_size];
			f1_mat[i] = new int[diagonal_size];
			f2_mat[i] = new int[diagonal_size];
		}

		init_first_triangle(d_mat, f1_mat, f2_mat, seq1, seq2);
		__m256i insert_penalty_vec = _mm256_set1_epi32(u[0]);

		__m256i u_1_vec = _mm256_set1_epi32(u[0]);
		__m256i u_2_vec = _mm256_set1_epi32(u[1]);
		__m256i v_1_vec = _mm256_set1_epi32(v[0]);
		__m256i v_2_vec = _mm256_set1_epi32(v[1]);

		//first half of diagonals
		for(int i = TRIANGLE_SIZE; i <= seq1.length(); ++i) {
			int last_ind = get_last_index(i, seq1, seq2);

			for(int j = 1; j <= last_ind; j += VECTOR_SIZE) {
				__m256i second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j]);
				//First part of convex function
				__m256i second_diagonal_f1 = _mm256_loadu_si256((__m256i *)&f1_mat[i - 1][j]);

				__m256i third_diagonal_f1 = _mm256_add_epi32(second_diagonal_d, v_1_vec);
				third_diagonal_f1 = _mm256_min_epi32(third_diagonal_f1, second_diagonal_f1);
				third_diagonal_f1 = _mm256_add_epi32(third_diagonal_f1, u_1_vec);
				__m256i third_diagonal_d = third_diagonal_f1;
				//Second part of convex function
				__m256i second_diagonal_f2 = _mm256_loadu_si256((__m256i *)&f2_mat[i - 1][j]);

				__m256i third_diagonal_f2 = _mm256_add_epi32(second_diagonal_d, v_2_vec);
				third_diagonal_f2 = _mm256_min_epi32(third_diagonal_f2, second_diagonal_f2);
				third_diagonal_f2 = _mm256_add_epi32(third_diagonal_f2, u_2_vec);
				third_diagonal_d = _mm256_min_epi32(third_diagonal_f2, third_diagonal_d);

				//Match/mismatch
				int m = last_ind - j + ((i > seq2.length()) ? (i - seq2.length()) : 0);
				int n = j - 1;

				__m256i diff_vec = _mm256_setr_epi32(test_coord_eligibility(m, n, seq1, seq2) ? diff(seq1[m - 1], seq2[n - 1]) : 0, test_coord_eligibility(m - 1, n + 1, seq1, seq2) ? diff(seq1[m - 2], seq2[n]) : 0,
						test_coord_eligibility(m - 2, n + 2, seq1, seq2) ? diff(seq1[m - 3], seq2[n + 1]) : 0, test_coord_eligibility(m - 3, n + 3, seq1, seq2) ? diff(seq1[m - 4], seq2[n + 2]) : 0,
						test_coord_eligibility(m - 4, n + 4, seq1, seq2) ? diff(seq1[m - 5], seq2[n + 3]) : 0, test_coord_eligibility(m - 5, n + 5, seq1, seq2) ? diff(seq1[m - 6], seq2[n + 4]) : 0,
						test_coord_eligibility(m - 6, n + 6, seq1, seq2) ? diff(seq1[m - 7], seq2[n + 5]) : 0, test_coord_eligibility(m - 7, n + 7, seq1, seq2) ? diff(seq1[m - 8], seq2[n + 6]) : 0);

				__m256i first_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 2][j - 1]);

				__m256i third_diagonal_f0 = _mm256_add_epi32(first_diagonal_d, diff_vec);

				third_diagonal_d = _mm256_min_epi32(third_diagonal_f0, third_diagonal_d);

				//Insertion
				second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j - 1]);

				__m256i third_diagonal_ins = _mm256_add_epi32(second_diagonal_d, insert_penalty_vec);
				third_diagonal_d = _mm256_min_epi32(third_diagonal_ins, third_diagonal_d);

				//storing results in matrices
				_mm256_storeu_si256((__m256i *)&d_mat[i][j], third_diagonal_d);
				_mm256_storeu_si256((__m256i *)&f1_mat[i][j], third_diagonal_f1);
				_mm256_storeu_si256((__m256i *)&f2_mat[i][j], third_diagonal_f2);
				
			}
			
			d_mat[i][0] = f1_mat[i][0] = f2_mat[i][0] = d_mat[i][last_ind + 1] = f1_mat[i][last_ind + 1] = f2_mat[i][last_ind + 1] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;

		}

		//second half of diagonals
		for(int i = seq1.length() + 1; i < matrix_row_num; ++i) {
			int last_ind = get_last_index(i, seq1, seq2);

			for(int j = 1; j <= last_ind; j += VECTOR_SIZE) {
				__m256i second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j + 1]);

				//First part of convex function
				__m256i second_diagonal_f1 = _mm256_loadu_si256((__m256i *)&f1_mat[i - 1][j + 1]);

				__m256i third_diagonal_f1 = _mm256_add_epi32(second_diagonal_d, v_1_vec);
				third_diagonal_f1 = _mm256_min_epi32(third_diagonal_f1, second_diagonal_f1);
				third_diagonal_f1 = _mm256_add_epi32(third_diagonal_f1, u_1_vec);
				__m256i third_diagonal_d = third_diagonal_f1;

				//Second part of convex function
				__m256i second_diagonal_f2 = _mm256_loadu_si256((__m256i *)&f2_mat[i - 1][j + 1]);

				__m256i third_diagonal_f2 = _mm256_add_epi32(second_diagonal_d, v_2_vec);
				third_diagonal_f2 = _mm256_min_epi32(third_diagonal_f2, second_diagonal_f2);
				third_diagonal_f2 = _mm256_add_epi32(third_diagonal_f2, u_2_vec);
				third_diagonal_d = _mm256_min_epi32(third_diagonal_f2, third_diagonal_d);

				//Match/mismatch
				int m = seq1.length() + 1 - j;
				int n = seq2.length() - last_ind + j;

				__m256i diff_vec = _mm256_setr_epi32(test_coord_eligibility(m, n, seq1, seq2) ? diff(seq1[m - 1], seq2[n - 1]) : 0, test_coord_eligibility(m - 1, n + 1, seq1, seq2) ? diff(seq1[m - 2], seq2[n]) : 0,
						test_coord_eligibility(m - 2, n + 2, seq1, seq2) ? diff(seq1[m - 3], seq2[n + 1]) : 0, test_coord_eligibility(m - 3, n + 3, seq1, seq2) ? diff(seq1[m - 4], seq2[n + 2]) : 0,
						test_coord_eligibility(m - 4, n + 4, seq1, seq2) ? diff(seq1[m - 5], seq2[n + 3]) : 0, test_coord_eligibility(m - 5, n + 5, seq1, seq2) ? diff(seq1[m - 6], seq2[n + 4]) : 0,
						test_coord_eligibility(m - 6, n + 6, seq1, seq2) ? diff(seq1[m - 7], seq2[n + 5]) : 0, test_coord_eligibility(m - 7, n + 7, seq1, seq2) ? diff(seq1[m - 8], seq2[n + 6]) : 0);


				__m256i first_diagonal_d = _mm256_loadu_si256((__m256i *)((i == seq1.length() + 1) ? &d_mat[i - 2][j] : &d_mat[i - 2][j + 1]));

				__m256i third_diagonal_f0 = _mm256_add_epi32(first_diagonal_d, diff_vec);
				third_diagonal_d = _mm256_min_epi32(third_diagonal_f0, third_diagonal_d);

				//Insertion
				second_diagonal_d = _mm256_loadu_si256((__m256i *)&d_mat[i - 1][j]);

				__m256i third_diagonal_ins = _mm256_add_epi32(second_diagonal_d, insert_penalty_vec);
				third_diagonal_d = _mm256_min_epi32(third_diagonal_ins, third_diagonal_d);

				//storing results in matrices
				_mm256_storeu_si256((__m256i *)&d_mat[i][j], third_diagonal_d);
				_mm256_storeu_si256((__m256i *)&f1_mat[i][j], third_diagonal_f1);
				_mm256_storeu_si256((__m256i *)&f2_mat[i][j], third_diagonal_f2);
			}

		}

		//printf("res je: %d\n", d_mat[matrix_row_num - 1][1]);

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
