#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <algorithm>
#include <map>
#include <fstream>
#include <memory>
#include <cmath>

#include "white_minimizers.h"

namespace white {

//pomocna fja samo za dobit stari prikaz kmera da mogu provjerit jel dobro
//ne koristi se osim za provjeru
/*unsigned int get_kmer(unsigned int number, unsigned int k) {
	unsigned int mask = 0x00000003;
	unsigned int result = 0;
	unsigned int digit = 0;

	for (unsigned int i = 0; i < k; i++) {
		digit = (number & mask) + 1;

		result = result + digit * std::pow (10, i);

		number = number >> 2;

	}

	return result;
}*/

std::vector<std::tuple<unsigned long long int, unsigned long long int, bool>> minimizers(
	const char* sequence,
	unsigned long long int sequence_length,
	unsigned long long int k,
	unsigned long long int window_length) {

	if (k > sequence_length) {
		std::cout << "K is bigger than the actual sequence." << std::endl;
		exit(1);
	}
	else if (k > 32) {
		std::cout << "Program doesn't work for k > 32," << std::endl;
		exit(1);
	}
	else if (window_length > sequence_length - k + 1) {
		std::cout << "Window size is bigger than the number of k-mers." << std::endl;
		exit(1);
	}

	//if (sequence_length >= 4294967295)
	//	exit(1);

	std::vector<std::tuple<unsigned long long int, unsigned long long int, bool>> minimizers;
	minimizers.reserve(sequence_length);
	unsigned long long int current_kmer_original = 0;
	unsigned long long int current_kmer_reverse_complement = 0;
	unsigned long long int next_window_starting_kmer_original = 0;
	unsigned long long int next_window_starting_kmer_reverse_complement = 0;
	unsigned long long int mask = 0;

	//stvorimo masku
	for (unsigned long long int i = 0; i < k; i++)
		mask = (mask << 2) + 3;

	//stvorimo pocetne k-mere

	for (unsigned long long int i = 0; i < k; i++) {
		switch (sequence[i]) {
			case 'A':
				current_kmer_original = (current_kmer_original << 2) & mask;
				break;
			case 'C':
				current_kmer_original = ((current_kmer_original << 2) | 1) & mask;
				break;
			case 'G':
				current_kmer_original = ((current_kmer_original << 2) | 2) & mask;
				break;
			case 'T':
				current_kmer_original = ((current_kmer_original << 2) | 3) & mask;
				break;

		}

		switch (sequence[sequence_length - 1 - i]) {
			case 'A':
				current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 3)& mask;
				break;
			case 'C':
				current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 2) & mask;
				break;
			case 'G':
				current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 1) & mask;
				break;
			case 'T':
				current_kmer_reverse_complement = (current_kmer_reverse_complement << 2) & mask;
				break;

		}
	}

	unsigned long long int last_minimizer = 0;
	long long int position_of_last_minimizer = -1;
	bool isOriginal = true;

	for (long long int i = 0; (unsigned) i < sequence_length - (window_length - 1 + k) + 1; i++) { // i mi treba kao int zbog pocetne pozicije koja je -1

		//pocetni kmeri za slijedeci prozor

		switch (sequence[i + k]) {
			case 'A':
				next_window_starting_kmer_original = (current_kmer_original << 2) & mask;
				break;
			case 'C':
				next_window_starting_kmer_original = ((current_kmer_original << 2) | 1) & mask;
				break;
			case 'G':
				next_window_starting_kmer_original = ((current_kmer_original << 2) | 2) & mask;
				break;
			case 'T':
				next_window_starting_kmer_original = ((current_kmer_original << 2) | 3) & mask;
				break;
		}

		switch (sequence[sequence_length - 1 - (i + k)]) {
			case 'A':
				next_window_starting_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 3) & mask;
				break;
			case 'C':
				next_window_starting_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 2) & mask;
				break;
			case 'G':
				next_window_starting_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 1) & mask;
				break;
			case 'T':
				next_window_starting_kmer_reverse_complement = (current_kmer_reverse_complement << 2) & mask;
				break;
		}

		//pomakom na slijedeci prozor minimizator ispada -> treba proci cijeli prozor
		if ((isOriginal && position_of_last_minimizer < i) || (!isOriginal && position_of_last_minimizer > (sequence_length-1) - i - (k-1))) {
			//postavimo prvi kmer kao minimizator
			if (current_kmer_original <= current_kmer_reverse_complement) {
				last_minimizer = current_kmer_original;
				position_of_last_minimizer = i;
				isOriginal = true;
			}
			else {
				last_minimizer = current_kmer_reverse_complement;
				position_of_last_minimizer = (sequence_length - 1) - i - (k - 1);
				isOriginal = false;
			}

			unsigned long long int end = i + k + window_length - 1;

			for (unsigned long long int j = i + k; j < end; j++) {
				//radimo slijedeci kmer u nizu (za original i reverse_complement)
				switch (sequence[j]) {
					case 'A':
						current_kmer_original = (current_kmer_original << 2) & mask;
						break;
					case 'C':
						current_kmer_original = ((current_kmer_original << 2) | 1) & mask;
						break;
					case 'G':
						current_kmer_original = ((current_kmer_original << 2) | 2) & mask;
						break;
					case 'T':
						current_kmer_original = ((current_kmer_original << 2) | 3) & mask;
						break;
				}

				switch (sequence[sequence_length - 1 - j]) {
					case 'A':
						current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 3) & mask;
						break;
					case 'C':
						current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 2) & mask;
						break;
					case 'G':
						current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 1) & mask;
						break;
					case 'T':
						current_kmer_reverse_complement = (current_kmer_reverse_complement << 2) & mask;
						break;
				}

				//usporedba
				if (last_minimizer >= current_kmer_original) { // radimo zamjenu ako su jednaki zato sto na taj nacin uzimamo minimizator
					last_minimizer = current_kmer_original;	// koji je na vecem indeksu, pa ce nam u buduce on ostati dulje unutar prozora,
					position_of_last_minimizer = j - (k - 1);	//tj. trebat ce proci vise prozora da bi on ispao
					isOriginal = true;
				}
				if (last_minimizer >= current_kmer_reverse_complement) {
					last_minimizer = current_kmer_reverse_complement;
					position_of_last_minimizer = (sequence_length - 1) - j /*+ (k - 1) - (k - 1)*/; //jednako kao gore jer bi trebali gledat obrnuto
					isOriginal = false;
				}
			}

		}
		//u suprotnom treba samo usporedit sa kranje desnim kmerom prozora
		else {
			current_kmer_original = 0;
			current_kmer_reverse_complement = 0;

			unsigned long long int start = i + window_length - 1; //trebat ce posli
			unsigned long long int end = i + k + window_length - 1; //zbog brzine

			//napravimo zadnji kmer (original i reverse complement)
			for (unsigned long long int j = start; j < end; j++) {
				switch (sequence[j]) {
					case 'A':
						current_kmer_original = (current_kmer_original << 2) & mask;
						break;
					case 'C':
						current_kmer_original = ((current_kmer_original << 2) | 1) & mask;
						break;
					case 'G':
						current_kmer_original = ((current_kmer_original << 2) | 2) & mask;
						break;
					case 'T':
						current_kmer_original = ((current_kmer_original << 2) | 3) & mask;
						break;

				}

				switch (sequence[sequence_length - 1 - j]) {
					case 'A':
						current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 3) & mask;
						break;
					case 'C':
						current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 2) & mask;
						break;
					case 'G':
						current_kmer_reverse_complement = ((current_kmer_reverse_complement << 2) | 1) & mask;
						break;
					case 'T':
						current_kmer_reverse_complement = (current_kmer_reverse_complement << 2) & mask;
						break;
				}
			}

			if (last_minimizer >= current_kmer_original) {
				last_minimizer = current_kmer_original;
				position_of_last_minimizer = start;
				isOriginal = true;
			}
			if (last_minimizer >= current_kmer_reverse_complement) {
				last_minimizer = current_kmer_reverse_complement;
				position_of_last_minimizer = (sequence_length - 1) - start - (k - 1);
				isOriginal = false;
			}
		}

		//pripremimo se za slijedeci prozor
		current_kmer_original = next_window_starting_kmer_original;
		current_kmer_reverse_complement = next_window_starting_kmer_reverse_complement;

		//iz pronadenog minimizatora stvorimo tuple i dodamo ga u vektor
		//std::tuple<unsigned int, unsigned int, bool> current_minimizer = std::make_tuple(get_kmer(last_minimizer, k), (unsigned int) position_of_last_minimizer, isOriginal); -> stari ispis
		std::tuple<unsigned long long int, unsigned long long int, bool> current_minimizer = std::make_tuple(last_minimizer, (unsigned long long int) position_of_last_minimizer, isOriginal);

		if (minimizers.empty())
			minimizers.emplace_back(current_minimizer);

		else if (minimizers[minimizers.size() - 1] != current_minimizer)
			minimizers.emplace_back(current_minimizer);

	}

	return minimizers;
}

}

