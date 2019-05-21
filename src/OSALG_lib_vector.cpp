#include <string>
#include <algorithm>
#include <vector>
#include <limits>
#include <unordered_map>
#include "OSALG_lib.h"

#define N_BORDER 30
#define L 2
#define MATCH_SCORE -2
#define MISMATCH_SCORE 4
#define DEFAULT_BANDWITH 16

#define STOP -1
#define MATCH 0
#define INSERT 1
#define DELETE 2


namespace OSALG {

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
			mat[i][j].d_val = mat[i][j].f_arr[0] = mat[i - 1][j - 1].d_val + diff(seq1[i - 1], seq2[j - 1]);
			mat[i][j].parent = MATCH;

			mat[i][j].f_arr[L + 1] = std::min(mat[i][j - 1].d_val + v[0], mat[i][j - 1].f_arr[L + 1]) + u[0];
			if (mat[i][j].d_val > mat[i][j].f_arr[L + 1]) {
				mat[i][j].d_val = mat[i][j].f_arr[L + 1];
				mat[i][j].parent = INSERT;
			}

			for (int k = 1; k < L + 1; ++k) {
				mat[i][j].f_arr[k] = std::min(mat[i - 1][j].d_val + v[k - 1], mat[i - 1][j].f_arr[k]) + u[k - 1];

				if(k > 1) {
					if (mat[i][j].d_val >= mat[i][j].f_arr[k]) {
						mat[i][j].d_val = mat[i][j].f_arr[k];
						mat[i][j].parent = DELETE;
					}
				} else {
					if (mat[i][j].d_val >= mat[i][j].f_arr[k]) {
						mat[i][j].d_val = mat[i][j].f_arr[k];
						mat[i][j].parent = DELETE;
					}	
				}
			}

		}
	}

	void init_first_triangle(std::vector<std::vector<matrix_element>> &mat, int bandwidth, std::string const &seq1, std::string const &seq2, diagonal_vector &first_diag, diagonal_vector &second_diag) {
		mat[0][0].f_arr[0] = mat[0][0].d_val = 0;
		for (int i = 1; i < L + 2; ++i) {
			mat[0][0].f_arr[i] = std::numeric_limits<int>::max();
		}
		mat[0][0].parent = STOP;
		for (int i = 0; i <= seq1.length(); ++i) {
			for (int j = 0; j <= seq2.length(); ++j) {
				if (i == 0 && j == 0) continue;
				calculate_matrix_element(mat, i, j, seq1, seq2);
				if(i == 54 && j == 64)
					printf("mat %d %d %d %d %d mat\n", mat[i][j].f_arr[0], mat[i][j].f_arr[1], mat[i][j].f_arr[2], mat[i][j].f_arr[3], mat[i][j].d_val);
			}
		}

		first_diag.x = 0;
		first_diag.y = bandwidth - 1;

		if ((mat[first_diag.x][first_diag.y].d_val <= mat[first_diag.x + bandwidth - 1][first_diag.y - bandwidth + 1].d_val && first_diag.x + bandwidth <= seq1.length()) || first_diag.y == seq2.length()) {
			second_diag.x = first_diag.x + 1;
			second_diag.y = first_diag.y;
		}
		else {
			second_diag.x = first_diag.x;
			second_diag.y = first_diag.y + 1;
		}
	}

	void init_last_triangle(std::vector<std::vector<matrix_element>> &mat, int bandwidth, std::string const &seq1, std::string const &seq2) {
		int c = 0;
		for (int i = seq1.length() - bandwidth + 2; i <= seq1.length(); ++i) {
			for (int j = seq2.length() - c; j <= seq2.length(); ++j) {
				calculate_matrix_element(mat, i, j, seq1, seq2);
			}
			c++;
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

	void construct_CIGAR(std::vector<std::vector<matrix_element>> const &matrix, std::string const &seq1, std::string const &seq2, std::string &cigar) {
		unsigned int i = seq1.length();
		unsigned int j = seq2.length();

		bool firstIdentified = false;
		unsigned int counter;
		char lastChar, c;

		while (matrix[i][j].parent != STOP) {
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

	int long_gaps_alignment(std::string const &seq1, std::string const &seq2, std::string &cigar, bool extended_cigar) {

		std::vector<std::vector<matrix_element>> mat = std::vector<std::vector<matrix_element>>(seq1.length() + 1);
		for (int i = 0; i < seq1.length() + 1; ++i) {
			mat[i] = std::vector<matrix_element>(seq2.length() + 1);
		}

		int bandwidth = std::min(seq1.length(), seq2.length()) + 1;
		diagonal_vector first_diag;
		diagonal_vector second_diag;
		diagonal_vector third_diag;

		init_first_triangle(mat, bandwidth, seq1, seq2, first_diag, second_diag);

		/*
		while (!(third_diag.x == seq1.length() - bandwidth + 1 && third_diag.y == seq2.length())) {
			bool downward = false;

			if ((mat[second_diag.x][second_diag.y].d_val <= mat[second_diag.x + bandwidth - 1][second_diag.y - bandwidth + 1].d_val &&
				(second_diag.x + bandwidth) <= seq1.length()) || (second_diag.y) == seq2.length()) {
				downward = true;
				third_diag.x = second_diag.x + 1;
				third_diag.y = second_diag.y;
			}
			else {
				third_diag.x = second_diag.x;
				third_diag.y = second_diag.y + 1;
			}

			for (int h = 0; h < bandwidth; ++h) {
				int i = third_diag.x + h;
				int j = third_diag.y - h;

				if (h == 0) {
					if (downward) {
						calculate_matrix_element(mat, i, j, seq1, seq2);
					}
					else {
						mat[i][j].f_arr[0] = mat[i][j].d_val = (first_diag.x == third_diag.x - 1 && first_diag.y == third_diag.y - 1) ? mat[i - 1][j - 1].d_val + diff(seq1[i - 1], seq2[j - 1]) : std::numeric_limits<int>::max();
						mat[i][j].parent = MATCH;

						for (int k = 1; k <= L; ++k) {
							mat[i][j].f_arr[k] = std::numeric_limits<int>::max();
						}

						mat[i][j].f_arr[L + 1] = std::min(mat[i][j - 1].d_val + v[0], mat[i][j - 1].f_arr[L + 1]) + u[0];

						if (mat[i][j].d_val > mat[i][j].f_arr[L + 1]) {
							mat[i][j].d_val = mat[i][j].f_arr[L + 1];
							mat[i][j].parent = INSERT;
						}
					}
				}
				else if (h == bandwidth - 1) {
					if (downward) {
						mat[i][j].f_arr[0] = mat[i][j].d_val = (first_diag.x + bandwidth - 1 == third_diag.x - 1 + bandwidth - 1 && first_diag.y - bandwidth + 1 == third_diag.y - 1 - bandwidth + 1) ? mat[i - 1][j - 1].d_val + diff(seq1[i - 1], seq2[j - 1]) : std::numeric_limits<int>::max();

						mat[i][j].parent = MATCH;

						for (int k = 1; k <= L; ++k) {
							mat[i][j].f_arr[k] = std::min(mat[i - 1][j].d_val + v[k - 1], mat[i - 1][j].f_arr[k]) + u[k - 1];

							if (mat[i][j].d_val >= mat[i][j].f_arr[k]) {
								mat[i][j].d_val = mat[i][j].f_arr[k];
								mat[i][j].parent = DELETE;
							}
						}

						mat[i][j].f_arr[L + 1] = std::numeric_limits<int>::max();
					}
					else {
						calculate_matrix_element(mat, i, j, seq1, seq2);
					}
				}
				else {
					calculate_matrix_element(mat, i, j, seq1, seq2);
				}
			}
			first_diag = second_diag;
			second_diag = third_diag;
		}
		*/
		//init_last_triangle(mat, bandwidth, seq1, seq2);

		construct_CIGAR(mat, seq1, seq2, cigar);

		return 0;
	}
}
