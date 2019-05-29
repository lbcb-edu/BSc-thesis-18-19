#include <string>
#include <algorithm>
#include <vector>
#include <limits>
#include <unordered_map>
#include <immintrin.h> //AVX

#include "OSALG_lib.h"

#define N_BORDER 30
#define L 2
#define MATCH_SCORE -2
#define MISMATCH_SCORE 4
#define DEFAULT_BANDWITH 16
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

	struct mat_elem {
		int d_val;
		int f_arr[L + 2];
		short parent;
	} typedef matrix_element;

	struct diag_vec {
		int x = 0;
		int y = 0;
	} typedef diagonal_vector;

	std::unordered_map<int, char> CIGAR_map = {
			{MATCH, 'M'},
			{INSERT, 'I'},
			{DELETE, 'D'}
	};

	int diff(char first, char second) {
		return (first == second) ? MATCH_SCORE : MISMATCH_SCORE;
	}

	void calculate_matrix_element(std::vector<std::vector<matrix_element>> &mat, int i, int j, std::string const &seq1, std::string const &seq2) {
		if (i == 0) {
			for (int k = 0; k < L + 1; ++k) {
				mat[i][j].f_arr[k] = std::numeric_limits<int>::max();
			}
			mat[i][j].f_arr[L + 1] = mat[i][j].d_val = std::min(mat[i][j - 1].d_val + v[0], mat[i][j - 1].f_arr[L + 1]) + u[0];

			mat[i][j].parent = INSERT;
		}
		else if (j == 0) {
			mat[i][j].f_arr[0] = mat[i][j].f_arr[L + 1] = mat[i][j].d_val = std::numeric_limits<int>::max();

			for (int k = 1; k < L + 1; ++k) {
				mat[i][j].f_arr[k] = mat[i][j].d_val = std::min(mat[i - 1][j].d_val + v[k - 1], mat[i - 1][j].f_arr[k]) + u[k - 1];
				if (mat[i][j].d_val > mat[i][j].f_arr[k]) mat[i][j].d_val = mat[i][j].f_arr[k];
			}

			mat[i][j].parent = DELETE;
		}
		else {
			mat[i][j].d_val = std::numeric_limits<int>::max();

			for (int k = 1; k < L + 1; ++k) {
				mat[i][j].f_arr[k] = std::min(mat[i - 1][j].d_val + v[k - 1], mat[i - 1][j].f_arr[k]) + u[k - 1];

				if(mat[i][j].d_val >= mat[i][j].f_arr[k]) {
					mat[i][j].d_val = mat[i][j].f_arr[k];
					mat[i][j].parent = (k == 1) ? DELETE1 : DELETE2;
				}
			}

			mat[i][j].f_arr[L + 1] = mat[i][j - 1].d_val + u[0];

			if(mat[i][j].d_val > mat[i][j].f_arr[L + 1]) {
				mat[i][j].d_val = mat[i][j].f_arr[L + 1];
				mat[i][j].parent = INSERT;
			}

			mat[i][j].f_arr[0] = mat[i - 1][j - 1].d_val + diff(seq1[i - 1], seq2[j - 1]);

			if(mat[i][j].d_val > mat[i][j].f_arr[0]) {
				mat[i][j].d_val = mat[i][j].f_arr[0];
				mat[i][j].parent = MATCH;
			}

		}
	}

	void init_first_triangle(std::vector<std::vector<matrix_element>> &mat, std::string const &seq1, std::string const &seq2) {
		mat[0][0].f_arr[0] = mat[0][0].d_val = 0;
		for (int i = 1; i < L + 2; ++i) {
			mat[0][0].f_arr[i] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
		}
		mat[0][0].parent = STOP;
		for (int i = 0; i <= TRIANGLE_SIZE; ++i) {
			for (int j = 0; j <= TRIANGLE_SIZE - i; ++j) {
				if (i == 0 && j == 0) continue;
				calculate_matrix_element(mat, i, j, seq1, seq2);
			}
		}
	}

	void init_last_triangle(std::vector<std::vector<matrix_element>> &mat, std::string const &seq1, std::string const &seq2) {
		int c = 0;
		for (int i = seq1.length() - TRIANGLE_SIZE + 1; i <= seq1.length(); ++i) {
			for (int j = seq2.length() - c; j <= seq2.length(); ++j) {
				calculate_matrix_element(mat, i, j, seq1, seq2);
			}
			++c;
		}
	}

	void updatePosition(unsigned int& i, unsigned int& j, int parent) {
		switch (parent) {
		case MATCH:
			i = i - 1;
			j = j - 1;
			break;
		case INSERT:
			j = j - 1;
			break;
		case DELETE:
			i = i - 1;
			break;
		}
	}

	void construct_CIGAR(std::vector<std::vector<matrix_element>> &matrix, std::string const &seq1, std::string const &seq2, std::string &cigar) {
		unsigned int i = seq1.length();
		unsigned int j = seq2.length();

		bool firstIdentified = false;
		unsigned int counter;
		char lastChar, c;

		bool del_mode = false;

		while (matrix[i][j].parent != STOP) {

			if(matrix[i][j].parent == DELETE1) {
				del_mode = false;
				matrix[i][j].parent = DELETE;
			} else if(matrix[i][j].parent == DELETE2) {
				del_mode = true;
				matrix[i][j].parent = DELETE;
			}

			if(del_mode && i > 0) {
				matrix[i][j].parent = DELETE;
			}

			if (matrix[i][j].parent == MATCH) {
				if (seq2[j - 1] == seq1[i - 1]) {
					c = '=';
				}
				else {
					c = 'X';
				}
			}
			else
				c = CIGAR_map.at(matrix[i][j].parent);

			if (firstIdentified && c == lastChar) {
				counter++;
			}
			else {
				if (firstIdentified) {
					cigar = std::to_string(counter) + lastChar + cigar;
				}
				else {
					firstIdentified = true;
				}

				lastChar = c;
				counter = 1;
			}
			updatePosition(i, j, matrix[i][j].parent);
		}

		cigar = std::to_string(counter) + lastChar + cigar;
	}

	void fill_array_from_diagonal_dvals(int i, int j, std::vector<std::vector<matrix_element>> const &mat, int *arr, std::string const &seq1, std::string const &seq2) {
		
		for(int k = 0; k < VECTOR_SIZE; ++k) {
			int i_pos = i + k;
			int j_pos = j - k;

			if(i_pos < 0 || i_pos  > seq1.length() || j_pos < 0 || j_pos > seq2.length()) {
				arr[k] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
			} else {
				arr[k] = mat[i + k][j - k].d_val;
			}
		}
	}

	void fill_array_from_diagonal_farr(int i, int j, std::vector<std::vector<matrix_element>> const &mat, int *arr, int level, std::string const &seq1, std::string const &seq2) {

		for(int k = 0; k < VECTOR_SIZE; ++k) {
			int i_pos = i + k;
			int j_pos = j - k;

			if(i_pos < 0 || i_pos  > seq1.length() || j_pos < 0 || j_pos > seq2.length()) {
				arr[k] = std::numeric_limits<int>::max() - OVERFLOW_CONTROL_VALUE;
			} else {
				arr[k] = mat[i + k][j - k].f_arr[level];
			}
		}
	}

	void fill_array_with_diff(int i, int j, int *arr, std::string const &seq1, std::string const &seq2) {

		for(int k = 0; k < VECTOR_SIZE; ++k) {
			int i_pos = i + k;
			int j_pos = j - k;

			if(i_pos < 0 || i_pos  > seq1.length() || j_pos < 0 || j_pos > seq2.length()) {
				arr[k] = 0;
			} else {
				arr[k] = (seq1[i_pos] == seq2[j_pos]) ? MATCH_SCORE : MISMATCH_SCORE;
			}
		}
	}

	int long_gaps_alignment(std::string const &seq1, std::string const &seq2, std::string &cigar, bool extended_cigar) {

		std::vector<std::vector<matrix_element>> mat = std::vector<std::vector<matrix_element>>(seq1.length() + 1);
		for (int i = 0; i < seq1.length() + 1; ++i) {
			mat[i] = std::vector<matrix_element>(seq2.length() + 1);
		}

		init_first_triangle(mat, seq1, seq2);

		__m256i third_diagonal_d;

		__m256i third_diagonal_f[L];

		//i and j store third diagonal "location"
		int i = 0;
		int j = TRIANGLE_SIZE + 1;

		__m256i insert_penalty_vec = _mm256_set1_epi32(u[0]);

		__m256i u_0_vec = _mm256_set1_epi32(u[0]);
		__m256i u_1_vec = _mm256_set1_epi32(u[1]);
		__m256i v_0_vec = _mm256_set1_epi32(v[0]);
		__m256i v_1_vec = _mm256_set1_epi32(v[1]);

		__m256i v_vecs[] = {v_0_vec, v_1_vec};
		__m256i u_vecs[] = {u_0_vec, u_1_vec};

		int temp_arr[VECTOR_SIZE];

		//Loops through diagonals
		while(i != seq1.length() - TRIANGLE_SIZE) {
			int k = 0;
			//Loop through specific diagonal
			while(i + k <= seq1.length() && j - k >= 0) {

				fill_array_from_diagonal_dvals(i + k - 1, j - k, mat, temp_arr, seq1, seq2);
				__m256i second_diagonal_d = _mm256_load_si256((__m256i *)&temp_arr[0]);

				//Going through F[0], ..., F[L]
				//DELETION
				for(int h = 0; h < L; ++h) {
					fill_array_from_diagonal_farr(i + k - 1, j - k, mat, temp_arr, h + 1, seq1, seq2);
					__m256i second_diagonal_f = _mm256_load_si256((__m256i *)&temp_arr[0]);


					third_diagonal_f[h] = _mm256_add_epi32(second_diagonal_d, v_vecs[h]);
					third_diagonal_f[h] = _mm256_min_epi32(third_diagonal_f[h], second_diagonal_f);
					third_diagonal_f[h] = _mm256_add_epi32(third_diagonal_f[h], u_vecs[h]);

					if(h == 0) {
						third_diagonal_d = third_diagonal_f[h];
					} else {
						third_diagonal_d = _mm256_min_epi32(third_diagonal_f[h], third_diagonal_d);
					}
				}

				//MATCH/MISMATCH
				fill_array_with_diff(i + k - 1, j - k - 1, temp_arr, seq1, seq2);
				__m256i diff_vec = _mm256_load_si256((__m256i *)&temp_arr[0]);

				fill_array_from_diagonal_dvals(i + k - 1, j - k - 1, mat, temp_arr, seq1, seq2);
				__m256i first_diagonal_d = _mm256_load_si256((__m256i *)&temp_arr[0]);

				__m256i third_diagonal_f0 = _mm256_add_epi32(first_diagonal_d, diff_vec);
				third_diagonal_d = _mm256_min_epi32(third_diagonal_f0, third_diagonal_d);
				
				//for(int l = 0; l < 8; ++l) {
				//	printf("temp1: %d\n", ((int*)&third_diagonal_d)[l]);
				//}

				//INSERTION
				fill_array_from_diagonal_dvals(i + k, j - k - 1, mat, temp_arr, seq1, seq2);
				second_diagonal_d = _mm256_load_si256((__m256i *)&temp_arr[0]);

				__m256i third_diagonal_ins = _mm256_add_epi32(second_diagonal_d, insert_penalty_vec);
				third_diagonal_d = _mm256_min_epi32(third_diagonal_ins, third_diagonal_d);
				
				int* third_diag_d_arr = (int*)&third_diagonal_d;
				int* third_diag_f_arr_0 = (int*)&third_diagonal_f[0];
				int* third_diag_f_arr_1 = (int*)&third_diagonal_f[1];
			
				for(int h = 0; h < VECTOR_SIZE; ++h) {
					if(i + k + h <= seq1.length() && j - k - h >= 0) {
						mat[i + k + h][j - k - h].d_val = third_diag_d_arr[h];
						mat[i + k + h][j - k - h].f_arr[1] = third_diag_f_arr_0[h];
						mat[i + k + h][j - k - h].f_arr[2] = third_diag_f_arr_1[h];
					} else {
						break;
					}
				}

				k += VECTOR_SIZE;
			}
			if(j < seq2.length()) {
				++j;
			} else {
				++i;
			}
		}
		
		init_last_triangle(mat, seq1, seq2);

		//construct_CIGAR(mat, seq1, seq2, cigar);
		return 0;
	}
}
