#include <iostream>
#include <memory>
#include <vector>
#include <getopt.h>
#include <stdlib.h>
#include <set>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <future>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>

#include "../src/white_minimizers.h"
#include "sequence_to_graph_mapping.h"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

#define SEQUENCE_ENDS_SIZE 2000
#define REGION_SIZE 4000
#define K_MER_LENGTH 10
#define WINDOW_SIZE 15
#define PERCENTAGE 0.001
#define ABUNDANCE_MIN 2
#define NUMBER_OF_THREADS 4
#define ACCURACY_OF_MAPPING 0.3

/*double time_spent_on_match = 0;
double time_spent_on_sorting = 0;
double time_spent_on_creating_regions = 0;
double time_spent_on_matching = 0;
double time_spent_on_dfs = 0;
double time_spent_on_creating_minimizers = 0;
*/
//---------------------------------FUNKCIJA ZA MJERENJE PROTEKLOG VREMENA------------------------------------

    /*typedef unsigned long long timestamp_t;

    static timestamp_t
    get_timestamp ()
    {
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }*/


//----------------------------------------------------------------------------------------------------------


unsigned long long int reference_genome_size = 0;

class SequenceFormat
{
	 public:
                std::string name;
                std::string sequence;
		std::string quality;

		friend void statistics (std::vector<std::unique_ptr<SequenceFormat>>&);
		friend bioparser::FastaParser<SequenceFormat>;
		friend bioparser::FastqParser<SequenceFormat>;

	private:
		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length
			) : name (name, name_length), sequence (sequence, sequence_length), quality (quality, quality_length)
                        {
			}

		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length
                        ) : name (name, name_length), sequence (sequence, sequence_length), quality ("")
                        {
                        }

};

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};


const struct option long_options[] = {
{"help",	no_argument,	0,	'h' },
{"version",	no_argument,	0,	'v' },
{0,		0,		0,	 0  },
};

void statistics (std::vector<std::unique_ptr<SequenceFormat>> &v, std::string file)
{
        uint32_t min_length = v.at(0) -> sequence.size();
        uint32_t max_length = min_length;
        uint32_t total_length = 0;

        for (auto &ptr : v)
        {
                if (ptr -> sequence.size() < min_length)
                        min_length = ptr -> sequence.size();

                if (ptr -> sequence.size() > max_length)
                        max_length = ptr -> sequence.size();

                total_length += ptr -> sequence.size();
        }

	printf( "\nStatistics for file: %s\n"
		"	Number of sequences: %ld\n"
		"	Minimum length: %d\n"
		"	Maximum length: %d\n"
		"	Average length: %.2f\n\n" ,
		file.c_str(), v.size(), min_length, max_length, (float)total_length/v.size() );

}

//TODO fix help and version
// Imported from a previous projects
void help() {
/*	std::cout << 		"--Help:"
				"\n	This program is a mapper that requires two files. First file should contain a set of	\n"
				"	fragments in FASTA or FASTQ format and the second one a corresponding reference genome	\n"
                                "	in FASTA format. The mapper parses the files and stores them in	memory while also	\n"
				"	providing statistics of given files:							\n"
                                "	-> number of sequences									\n"
                                "	-> average length									\n"
                                "	-> minimal and maximal length.								\n\n"

				"	Afterwards, two random sequences from the first file are selected and aligned using	\n"
				"	global, semi-global and local alignment. The resulting cigar string is printed.		\n"
				"	Default for match, mismatch ang gap values is 2, -1, -2 respectively. The default	\n"
				"	values can be changed using appropriate options.					\n\n"

				"	Seed and extend:									\n"
				"	In order to alleviate the alignment proces we use seed and extend approach. Among all	\n"
				"	k-mers we choose a set of minimizers which will be used for fast detection of similar	\n"
				"	regions between two sequences prior the exact alignment. The mapper will find minimizers\n"
				"	for every sequence in the first file and print the number of distinct minimizers and the\n"
				"	number of occurences of the most frequent minimizer when the top f frequent minimizers 	\n"
				"	are not taken in account.								\n\n"

				"	The proces of finding minimizers uses three variables:					\n"
				"	-> k-mer length k									\n"
				"	-> window size w									\n"
				"	-> percentage of top minimizers to disregard f						\n\n"

				"	Their default values are 15, 5 and 0.001 respecively.					\n\n"

                                "	File extensions accepted:								\n"
                                "	-> .fasta             -> .fastq								\n"
                                "	-> .fa                -> .fq								\n"
                                "	-> .fasta.gz          -> .fastq.gz							\n"
                                "	-> .fa.gz             -> .fq.gz								\n\n"

				"	Usage: white_mapper [OPTIONS] sequences_file reference_genome_file match mismatch gap	\n\n"

                                "	There are five options:									\n"
                                "	-> \"-h\" or \"--help\" for help							\n"
                                "	-> \"-v\" or \"--version\" for displaying the current version.				\n"
				"	-> \"-m ARG\" sets match value to ARG							\n"
				"       -> \"-s ARG\" sets mismatch value to ARG                                                \n"
				"       -> \"-g ARG\" sets gap value to ARG                                                     \n"
				"       -> \"-k ARG\" sets length k of k-mers to ARG						\n"
				"       -> \"-w ARG\" sets window size to ARG							\n"
				"	-> \"-f ARG\" sets f to ARG								\n\n";
*/}

void version() {
/*	std::cout << "\n--Version:";
	std::cout << " " << sequence_to_graph_mapping_VERSION_MAJOR << "." << sequence_to_graph_mapping_VERSION_MINOR << "\n\n";
*/}

bool check_format (std::string file, std::set<std::string> list_of_formats) {
	for (auto &format : list_of_formats)
		if(file.size() > format.size())
			if (file.compare (file.size() - format.size(), format.size(), format) == 0)
				return true;
	return false;
}

template < template<class> class T>
std::vector<std::unique_ptr<SequenceFormat>> parse_file (std::string file_path)
{
	std::vector<std::unique_ptr<SequenceFormat>> objects;
	auto parser = bioparser::createParser<T, SequenceFormat> (file_path);
	parser -> parse (objects, -1);

	statistics (objects, file_path);

	return objects;
}

//---------------------------------------------------------------OSTAVI_OVO_IZNAD------------------------------------------------------------------------------------------------------
/*bool comp(std::vector< std::tuple <unsigned long long int, unsigned long long int, bool>> a, std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> b) {
	return a.size() < b.size();
}

void get_reference_genome_minimizers (
	std::unique_ptr<SequenceFormat> &ref,
	unsigned int k,
	unsigned int window_size,
	float f,
	std::unordered_map <unsigned long long int, unsigned long long int > &minimizer_position_in_vector,
	std::vector <std::vector <std::tuple<unsigned long long int, unsigned long long int, bool>>> &occurrences
	) {

	std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> minimizers;

	minimizers = white::minimizers(ref -> sequence.c_str(), ref -> sequence.size(), k, window_size);

	unsigned long long int i = 0;

	std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> temp;

	for (auto minimizer_tuple : minimizers) {
		if (minimizer_position_in_vector[std::get<0>(minimizer_tuple)] == 0){
			i++;
			temp.emplace_back(minimizer_tuple);
			minimizer_position_in_vector[std::get<0>(minimizer_tuple)] = i;
			occurrences.emplace_back(temp);
			temp.clear();
		}
		else
			occurrences[minimizer_position_in_vector[std::get<0>(minimizer_tuple)] - 1].emplace_back(minimizer_tuple);
	}

	unsigned long long int number_of_minimizers_to_delete = (unsigned long long int) std::round(occurrences.size() * f);

	std::sort(occurrences.begin(), occurrences.end(), comp);

	for (unsigned int i = occurrences.size() - 1; i >= occurrences.size() - number_of_minimizers_to_delete; i++) {
		minimizer_position_in_vector.erase(std::get<0>(occurrences[i].at(0)));
		occurrences.erase(occurrences.begin() + i);

	}

}*/

bool comp(std::tuple <unsigned long long int, unsigned long long int> a, std::tuple <unsigned long long int, unsigned long long int> b) {
        return std::get<1>(a) < std::get<1>(b);
}

//funkcija koja stvara hash za matchanje minimizatora i regija

std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> create_hash_for_matching(
	std::vector <std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>>> &unitigs_and_transitions,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> &unitig_u_kojim_regijama,
	unsigned int k,
	unsigned int window_size
) {

	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> hash_minimizer_region_number_to_add;
	std::unordered_map <unsigned long long int, unsigned long long int> minimizer_count;


	unsigned long long int j = 0;

	for (auto tuple : unitigs_and_transitions) {

		std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> minimizers;

		minimizers = white::minimizers(std::get<0>(tuple).c_str(), std::get<0>(tuple).size(), k, window_size);

		for (auto minimizer_tuple : minimizers) {
			minimizer_count[std::get<0>(minimizer_tuple)]++;
			for (auto region : unitig_u_kojim_regijama[j])
				hash_minimizer_region_number_to_add[std::get<0>(minimizer_tuple)][region.first]++;
		}

		j++;

	}


	unsigned int number_of_minimizers_to_keep = (unsigned int) std::round(minimizer_count.size() * (1 - PERCENTAGE));

	//printf("\n%d\n", minimizer_count.size() - number_of_minimizers_to_keep);

	std::vector <std::tuple <unsigned long long int, unsigned long long int>> map_vector;

	for(auto elem : minimizer_count)
		map_vector.emplace_back(std::make_tuple (elem.first, elem.second));

        std::sort(map_vector.begin(), map_vector.end(), comp);

	for (unsigned long long int i = number_of_minimizers_to_keep; i < map_vector.size(); i++)
		hash_minimizer_region_number_to_add.erase(std::get<0>(map_vector[i]));

	return hash_minimizer_region_number_to_add;

}


std::vector <std::tuple < std::unordered_map <unsigned long long int, unsigned long long int >, std::vector < std::vector < std::tuple <unsigned long long int, unsigned long long int, bool>>>>> create_unitig_minimizers (
	std::vector <std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>>> &unitigs_and_transitions,
	unsigned int k,
	unsigned int window_size
) {

	std::vector< std::tuple < std::unordered_map <unsigned long long int, unsigned long long int >, std::vector < std::vector < std::tuple <unsigned long long int, unsigned long long int, bool>>>>> unitig_minimizer_map;

	for (auto tuple : unitigs_and_transitions) {

		std::unordered_map <unsigned long long int, unsigned long long int > unitig_minimizers_position_in_vector;
		std::vector <std::vector <std::tuple<unsigned long long int, unsigned long long int, bool>>> unitig_minimizers;

		std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> minimizers;

		minimizers = white::minimizers(std::get<0>(tuple).c_str(), std::get<0>(tuple).size(), k, window_size);

		unsigned long long int j = 0;

		std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> temp;

		for (auto minimizer_tuple : minimizers) {
			if (unitig_minimizers_position_in_vector[std::get<0>(minimizer_tuple)] == 0) {
				j++;
				temp.emplace_back(minimizer_tuple);
				unitig_minimizers_position_in_vector[std::get<0>(minimizer_tuple)] = j;
				unitig_minimizers.emplace_back(temp);
				temp.clear();
			}
			else
				unitig_minimizers[unitig_minimizers_position_in_vector[std::get<0>(minimizer_tuple)] - 1].emplace_back(minimizer_tuple);
		}

		unitig_minimizer_map.emplace_back(std::make_tuple(unitig_minimizers_position_in_vector, unitig_minimizers));

	}

	return unitig_minimizer_map;
}



bool compare_minimizers(std::tuple <unsigned long long int, unsigned long long int, bool> a, std::tuple <unsigned long long int, unsigned long long int, bool> b) {
	return std::get<0>(a) < std::get<0>(b);
}

int my_binary_search(std::vector <std::tuple <unsigned int, unsigned int, bool>> &array, unsigned int number_to_find) {
	int min = 0;
	int max = (array.size() - 1);
	int guess;
	int result = -1;

	while (min <= max)
	{
		guess = (int)(((max + min) / 2) + 0.5);

		if (number_to_find == std::get<0>(array[guess]))
		{

			//posto moze biti vise istih minimizatora (samo na razlicitim pozicijama) moramo se vratiti natrag do prvoga
			result = guess;

			while (result > 0 && std::get<0>(array[result - 1]) == number_to_find)
				result--;

			return result;
		}
		else if (std::get<0>(array[guess]) < number_to_find) {
			min = guess + 1;
		}
		else {
			max = guess - 1;
		}
	}

	return result;
}

bool minimizers_equal(unsigned long long int a, unsigned long long int b) {

	unsigned long long int temp_a;
	unsigned long long int temp_b;
	unsigned long long int mask = 0x0000000000000003;
	unsigned long long int number_of_hits = 0;

	for (int i = 0; i < K_MER_LENGTH; i++) {
		temp_a = a & mask;
		temp_b = b & mask;

		a = a >> 2;
		b = b >> 2;

		if (temp_a == temp_b)
			number_of_hits++;
	}

	if ((double) number_of_hits >= (double) K_MER_LENGTH * 0.9)
		return true;

	return false;

}

void minimizer_matches(
	std::unordered_map <unsigned long long int, unsigned long long int> &reference_genome_minimizers_positions_in_vector,
        std::vector < std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>>> &reference_genome_minimizers,
	std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> &sequence_minimizers,
	std::vector <std::tuple <unsigned long long int, unsigned long long int>> &matches_same_strand,
	std::vector <std::tuple <unsigned long long int, unsigned long long int>> &matches_different_strand
) {

	for (auto tuple : sequence_minimizers) {
		if (reference_genome_minimizers_positions_in_vector.find(std::get<0>(tuple)) != reference_genome_minimizers_positions_in_vector.end())
			for (auto tupleR : reference_genome_minimizers[reference_genome_minimizers_positions_in_vector[std::get<0>(tuple)] - 1]) {
				if (minimizers_equal(std::get<2>(tupleR), std::get<2>(tuple))) //if (std::get<2>(tupleR) == std::get<2>(tuple))
					matches_same_strand.emplace_back(std::make_tuple(std::get<1>(tupleR), std::get<1>(tuple)));
				else
					matches_different_strand.emplace_back(std::make_tuple(std::get<1>(tupleR), std::get<1>(tuple)));
			}

	}

}


//funkcija koja vraca indeks elementa koji je najmanji element veci ili jednak key (naravno funkcija funkcionira jedino ako je vektor
//sortiran i provjerili smo da je key nije veci od najveceg elementa vektora (tj. ne vrijedi: key > najveci element))
unsigned int CeilIndex(std::vector <std::tuple <long long int, long long int, long long int>> &v, int l, int r, long long int key) {
	while (r - l > 1) {
		int m = l + (r - l) / 2;

		if (std::get<2>(v[m]) >= key)
			r = m;
		else
			l = m;
	}

	return r;
}

//funkcija koja trazi LIS za zadane matcheve minimizatora (funkcija pazi da u LIS-u ne budu matchevi s istom vrijednoscu query -> pomocno polje)
std::pair<int, int> find_LIS(std::vector <std::tuple <long long int, long long int>> &minimizer_matches) {
	std::vector <std::tuple <long long int, long long int, long long int>> tail(minimizer_matches.size()); //polje za biljezit LIS-a
	//polja su redom:
	//	1. pozicija prvog matcha LIS-a u minimizer_matches
	//	2. pozicija zadnjeg matcha LIS-a u minimizer_matches
	//	3. zadnji element LIS-a

	std::vector <std::tuple <long long int, long long int, long long int>> auxiliary_tail(minimizer_matches.size()); //pomocno polje, ovo polje sluzi kako bi eliminirali
														//pojavljivanje matcheva koji imaju istu vrijednost query

	std::vector <long long int> changed_indexes; //vektor koji biljezi gdje su se dogodile promjene u tail tako da ih mozemo prepisati u auxiliary_tail
						//umjesto uvijek da kad predemo sve matcheve koji imaju isti query_position prepisemo cijeli tail u
						//auxiliary_tail mozemo samo prepisati promijenjene elemente
	if (minimizer_matches.size() == 0)
		return std::make_pair(-1, -1); //ako nema matcheva ne treba radit alignment

	long long int length = 1; //broji duljinu LIS-a
	long long int auxiliary_tail_length = 1; //velicina pomocnog polja

	//minimizer_matches su sortirani po query_position(prva vrijednost) pa po target_position(druga vrijednost)
	//kao prvi element LIS-a postavljamo match koji ima najmanju vrijednost target_position, ali
	//gledamo samo one elemente s prvom (najmanjom) vrijednost query_position

	long long int query_value_min = std::get<0>(minimizer_matches[0]);
	long long int first_target_value_position = 0;

	for (long long int i = 1; query_value_min == std::get<0>(minimizer_matches[i]); i++)
		if (std::get<1>(minimizer_matches[first_target_value_position]) < std::get<1>(minimizer_matches[i]))
			first_target_value_position = i;

	tail[0] = std::make_tuple(first_target_value_position, first_target_value_position, std::get<1>(minimizer_matches[first_target_value_position]));
	auxiliary_tail[0] = tail[0];

	//	for (int i = 0; i < length; i++)
	//		std::cout << std::get<0>(tail[i]) << std::get<1>(tail[i]) << std::get<2>(tail[i]) << std::endl;

	//std::cout << std::endl;

	//sad prolazimo kroz ostale matcheve i racunamo LIS
	//da bi eliminari pojavu elemenata s istom vrijednosti query koristimo pomocno polje

	for (unsigned int i = 1; i < minimizer_matches.size(); i++) {

		long long int target_position = std::get<1>(minimizer_matches[i]);

		//manji od najmanjeg -> samo zamjenimo najmanji element
		if (target_position <= std::get<2>(tail[0])) {
			tail[0] = std::make_tuple(i, i, target_position);
			changed_indexes.emplace_back(0);
		}
		//veci od najveceg: dvije opcije
		//1 -> ako najveci element nema isti query_position kao i novi element, onda dodamo novi element u vektor
		//2 -> ako najveci element ima isti query_position kao i novi element, onda provjerimo postoji li element na istoj poziciji
		//	-> kao i najveci element u pomocnom vektoru:
		//	-> ako da, onda stvorimo novi element u tail ako je on veci od elemnta u pomocnom vektoru
		//  -> ako ne, onda ne dodajemo nista
		else if (target_position > std::get<2>(tail[length - 1])) {
			if (std::get<0>(minimizer_matches[std::get<1>(tail[length - 1])]) != std::get<0>(minimizer_matches[i])) {
				tail[length] = std::make_tuple(std::get<0>(tail[length - 1]), i, target_position);
				changed_indexes.emplace_back(length);
				length++;
			}
			else if (length == auxiliary_tail_length && std::get<2>(auxiliary_tail[length - 1]) < target_position) {
				tail[length] = std::make_tuple(std::get<0>(auxiliary_tail[length - 1]), i, target_position);
				changed_indexes.emplace_back(length);
				length++;
			}
		}

		//ako nije nista od prethodnog, onda trazimo indeks najmanjeg veceg ili jednakog broja od novog elementa -> dobijemo indeks
		//promatramo element na indeks - 1
		//dvije opcije:
		//1 -> ako taj element nema isti query_position kao i novi element, onda zamijenimo element na indeksu: indeks sa novim elementom
		//2 -> ako taj element ima isti query_position kao i novi element, onda pogledamo u pomocnom vektoru element na poziciji indeks - 1
		//dvije opcije:
		//1 -> ako je taj element manji od novog elementa, ako zamijenimo u vektoru tail element na indeksu: indeks s novim elementom
		//2 -> ako ne, ne radimo nista

		else {
			long long int index = CeilIndex(tail, -1, length - 1, target_position);

			if (std::get<0>(minimizer_matches[std::get<1>(tail[index - 1])]) != std::get<0>(minimizer_matches[i])) {
				tail[index] = std::make_tuple(std::get<0>(tail[index - 1]), i, target_position);
				changed_indexes.emplace_back(index);
			}

			else if (std::get<2>(auxiliary_tail[index - 1]) < target_position) {
				tail[index] = std::make_tuple(std::get<0>(auxiliary_tail[index - 1]), i, target_position);
				changed_indexes.emplace_back(index);
			}

		}

		//nakon sto prodemo sve elemente s istim query_position moramo prekopirati tail u auxiliary_tail
		if (i + 1 < minimizer_matches.size() && std::get<0>(minimizer_matches[i]) != std::get<0>(minimizer_matches[i + 1])) {
			for (auto index : changed_indexes)
				auxiliary_tail[index] = tail[index];

			changed_indexes.clear();

			auxiliary_tail_length = length;
		}


		/*		for (int i = 0; i < length; i++)
					std::cout << std::get<0>(tail[i]) << std::get<1>(tail[i]) << std::get<2>(tail[i]) << std::endl;

				std::cout << std::endl;

				for (int i = 0; i < auxiliary_tail_length; i++)
					std::cout << std::get<0>(auxiliary_tail[i]) << std::get<1>(auxiliary_tail[i]) << std::get<2>(auxiliary_tail[i]) << std::endl;

				std::cout << std::endl;
		*/
	}

	if (length >= 4)
		return std::make_pair(std::get<0>(tail[length - 1]), std::get<1>(tail[length - 1]));
	else
		return std::make_pair(-1, -1);
}


std::string reverse_complement (std::string normal_strand) {
	std::string rc_string;

	for (int i = normal_strand.size() - 1; i >= 0; i--) {
		switch (normal_strand[i]) {
			case 'A':
				rc_string.push_back('T');
				break;

			case 'T':
                                rc_string.push_back('A');
                                break;

			case 'C':
                                rc_string.push_back('G');
                                break;

			case 'G':
                                rc_string.push_back('C');
                                break;
		}
	}

	return rc_string;
}


std::vector<unsigned long long int> create_regions_puna_regija(
	std::vector < std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>>> &unitigs_and_transitions,
	std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>> &unitig_tuple,
	unsigned long long int unitig_index,
	unsigned long long int starting_position_in_unitig,
	std::string strand,
	unsigned long long int regions_total_size,
	unsigned long long int current_size,
	std::unordered_map <unsigned long long int, std::tuple<bool, bool>> &processed_unitigs,
	std::vector <std::tuple <unsigned long long int, unsigned long long int, std::string, unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int>> &queue,
	bool processed,
	unsigned long long int &index_counter,
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> &regions,
	//std::unordered_map <std::vector<unsigned long long int>, unsigned long long int> &regions_auxiliary_map,
	unsigned long long int level,
	std::unordered_map <unsigned long long int, bool> &regija_processed,
	std::unordered_map <unsigned long long int, std::string> &regija_zadnjiStrand,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> &unitig_u_kojim_regijama
) {

	std::vector<unsigned long long int> region_unitigs;
	std::vector<unsigned long long int> neighbouring_regions;

	std::vector<unsigned long long int> created_regions;
	std::vector<unsigned long long int> all_created_regions;

	unsigned long long int pom = current_size;

	if (level == 0) {
		current_size += std::get<0>(unitig_tuple).size();
	}
	else {
		current_size += std::get<0>(unitig_tuple).size() - K_MER_LENGTH;
	}



	//std::string region_string;
	std::tuple <bool, bool> temp = std::make_tuple(false, false);

	/*for (auto it = processed_unitigs.begin(); it != processed_unitigs.end(); it++)
		printf("%llu\n", it->first);*/

		//printf("%llu\n", unitig_index);
		//printf("%llu\n", level);



	if (processed_unitigs.find(unitig_index) != processed_unitigs.end()) {
		processed = true;
		temp = processed_unitigs[unitig_index];
	}

	if (strand.compare("+") == 0) {
		processed_unitigs[unitig_index] = std::make_tuple(true, std::get<1>(temp));
	}
	else {
		processed_unitigs[unitig_index] = std::make_tuple(std::get<0>(temp), true);
	}

	//---------------------------------------------------------------------------------------------------

	if (processed) {
		index_counter++;

		regija_zadnjiStrand[index_counter] = strand;
		regija_processed[index_counter] = true;

		//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
		for (unsigned long long int i = 0; i <= level; i++)
			region_unitigs.emplace_back(555);

		//-------------------
		unsigned long long int rSize = 0;

		if (level == 0)
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;
		else if (starting_position_in_unitig < K_MER_LENGTH)
			rSize = pom + std::get<0>(unitig_tuple).size() - K_MER_LENGTH;
		else
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;

		if (rSize > REGION_SIZE)
			rSize = REGION_SIZE;

		//-------------------

		unsigned long long int remove = 0;

		if (current_size > REGION_SIZE)
			remove = current_size - REGION_SIZE;

		//-------------------SUSJEDNOST
		/*for (auto it = regions.begin(); it != regions.end(); it++)
			for (auto unitig : std::get<0>(it->second))
				if (unitig == unitig_index) {
					neighbouring_regions.emplace_back(it->first);
					std::get<3>(it->second).emplace_back(index_counter);
				}*/
		//-------------------\SUSJEDNOST

		region_unitigs[level] = unitig_index;
		regions[index_counter] = std::make_tuple(region_unitigs, 0, std::get<0>(unitig_tuple).size() - 1 - remove, neighbouring_regions, rSize);

		if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
			unitig_u_kojim_regijama[unitig_index][index_counter] = 1;

		all_created_regions.emplace_back(index_counter);

		//if (strand.compare("+") == 0) {
		//	all_created_regions.emplace_back(std::make_tuple(std::get<0>(unitig_tuple), index_counter));
		//}
		//else {
		//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple)), index_counter));
		//}

		return all_created_regions;
	}

	//---------------------------------------------------------------------------------------------------

	if (level != 0) {
		if (current_size > regions_total_size) {
			unsigned long long int remaining_size = regions_total_size - pom;
			index_counter++;

			regija_zadnjiStrand[index_counter] = strand;

			if (processed == true)
				regija_processed[index_counter] = true;

			//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
			for (unsigned long long int i = 0; i <= level; i++)
				region_unitigs.emplace_back(555);

			//printf("%llu %llu\n", regions_total_size, current_size);
			//printf("%llu\n", remaining_size);

			//if (strand.compare("+") == 0) {
			region_unitigs[level] = unitig_index;
			regions[index_counter] = std::make_tuple(region_unitigs, 0, remaining_size + K_MER_LENGTH, neighbouring_regions, REGION_SIZE);

			if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
				unitig_u_kojim_regijama[unitig_index][index_counter] = 1;

			all_created_regions.emplace_back(index_counter);
			//}
			//else {
			//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple).substr(0, remaining_size)), index_counter));
			//}

			if (pom < regions_total_size / 2 /*&& current_size >= regions_total_size / 2*/ && processed == false) {
				//if (processed_unitigs.find(unitig_index) != processed_unitigs.end())
				processed_unitigs.erase(unitig_index);

				queue.emplace_back(std::make_tuple(
					//std::get<0>(regions[region])[level],
					unitig_index,
					regions_total_size / 2 - pom + K_MER_LENGTH,
					strand,
					index_counter,
					0,
					0,
					regions_total_size - pom + K_MER_LENGTH));
			}

			return all_created_regions;
		}
	}

	else {
		if (current_size > regions_total_size) {
			unsigned long long int remaining_size = regions_total_size;
			index_counter++;

			regija_zadnjiStrand[index_counter] = strand;

			if (processed == true)
				regija_processed[index_counter] = true;

			//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
			for (unsigned long long int i = 0; i <= level; i++)
				region_unitigs.emplace_back(555);

			//printf("%llu %llu\n", regions_total_size, current_size);
			//printf("%llu\n", remaining_size);

			//if (strand.compare("+") == 0) {
			region_unitigs[level] = unitig_index;
			regions[index_counter] = std::make_tuple(region_unitigs, 0, remaining_size, neighbouring_regions, REGION_SIZE);

			if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
				unitig_u_kojim_regijama[unitig_index][index_counter] = 1;

			all_created_regions.emplace_back(index_counter);
			//}
			//else {
			//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple).substr(0, remaining_size)), index_counter));
			//}

			if (pom < regions_total_size / 2/*&& current_size >= regions_total_size / 2 && processed == false*/) {
				//if (processed_unitigs.find(unitig_index) != processed_unitigs.end())
				processed_unitigs.erase(unitig_index);

				queue.emplace_back(std::make_tuple(
					//std::get<0>(regions[region])[level],
					unitig_index,
					regions_total_size / 2,
					strand,
					index_counter,
					0,
					0,
					regions_total_size));
			}

			return all_created_regions;
		}
	}


	bool no_connections = true;

	for (auto &connected_unitig : std::get<1>(unitig_tuple)) {
		if (strand.compare(std::get<0>(connected_unitig)) == 0) {
			no_connections = false;
			break;
		}
	}
	if (no_connections) {
		index_counter++;

		regija_zadnjiStrand[index_counter] = strand;
		regija_processed[index_counter] = true;

		//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
		for (unsigned long long int i = 0; i <= level; i++)
			region_unitigs.emplace_back(555);

		//-------------------
		unsigned long long int rSize = 0;

		if (level == 0)
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;
		else if (starting_position_in_unitig < K_MER_LENGTH)
			rSize = pom + std::get<0>(unitig_tuple).size() - K_MER_LENGTH;
		else
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;

		//-------------------

		region_unitigs[level] = unitig_index;
		regions[index_counter] = std::make_tuple(region_unitigs, 0, std::get<0>(unitig_tuple).size() - 1, neighbouring_regions, rSize);

		if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
			unitig_u_kojim_regijama[unitig_index][index_counter] = 1;

		all_created_regions.emplace_back(index_counter);

		//if (strand.compare("+") == 0) {
		//	all_created_regions.emplace_back(std::make_tuple(std::get<0>(unitig_tuple), index_counter));
		//}
		//else {
		//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple)), index_counter));
		//}

		return all_created_regions;
	}



	for (auto &connected_unitig : std::get<1>(unitig_tuple)) {

		//mislim da ne treba ovaj dio
		/*temp = std::make_tuple (false, false);

		if (processed_unitigs.find(std::get<1>(connected_unitig)) != processed_unitigs.end()) {
					temp = processed_unitigs[connected_unitig];
			}

			if (std::get<2>(connected_unitig).compare("+") == 0) {
					processed_unitigs[connected_unitig] = std::make_tuple (true, std::get<1>(temp));
			}else {
					processed_unitigs[connected_unitig] = std::make_tuple (std::get<0>(temp), true);
			}*/

		if (strand.compare(std::get<0>(connected_unitig)) == 0) {

			created_regions = create_regions_puna_regija(
				unitigs_and_transitions,
				unitigs_and_transitions[std::get<1>(connected_unitig)],
				std::get<1>(connected_unitig),
				0,
				std::get<2>(connected_unitig),
				regions_total_size,
				current_size,
				processed_unitigs,
				queue,
				processed,
				index_counter,
				regions,
				//regions_auxiliary_map,
				level + 1,
				regija_processed,
				regija_zadnjiStrand,
				unitig_u_kojim_regijama
			);

			for (auto &region : created_regions) {

				std::get<0>(regions[region])[level] = unitig_index;

				if (unitig_u_kojim_regijama[unitig_index][region] == 0)
					unitig_u_kojim_regijama[unitig_index][region] = 1;

				all_created_regions.emplace_back(region);

				//printf("%llu %llu %llu %d\n", pom, current_size, regions_total_size / 2, std::get<0>(unitig_tuple).size());
				if (pom < regions_total_size / 2 && current_size >= regions_total_size / 2 && regija_processed.find(region) == regija_processed.end()) {

					unsigned long long int unitig_index_temp = std::get<0>(regions[region])[std::get<0>(regions[region]).size() - 1];

					//for (int i = level; i <= std::get<0>(regions[region]).size() - 1; i++)
					//	processed_unitigs.erase(std::get<0>(regions[region])[i]);


					if (processed_unitigs.find(unitig_index_temp) != processed_unitigs.end())
						processed_unitigs.erase(unitig_index_temp);
					unsigned long long int starting_pos_temp;

					if (level == 0)
						starting_pos_temp = regions_total_size / 2 - pom;
					else
						starting_pos_temp = regions_total_size / 2 - pom + K_MER_LENGTH;

					//unsigned long long int pom = std::get<0>(regions[region]).size() - level - 1;

					queue.emplace_back(std::make_tuple(
						//std::get<0>(regions[region])[level],
						unitig_index_temp,
						starting_pos_temp,
						regija_zadnjiStrand[region],
						region,
						std::get<0>(regions[region]).size() - level - 1,
						level,
						std::get<2>(regions[region])));
				}
			}
		}

	}

	return all_created_regions;



}

std::vector<unsigned long long int> create_regions_pola_regije(
	std::vector < std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>>> &unitigs_and_transitions,
	std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>> &unitig_tuple,
	unsigned long long int unitig_index,
	unsigned long long int starting_position_in_unitig,
	std::string strand,
	unsigned long long int regions_total_size,
	unsigned long long int current_size,
	std::unordered_map <unsigned long long int, std::tuple<bool, bool>> &processed_unitigs,
	std::vector <std::tuple <unsigned long long int, unsigned long long int, std::string, unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int>> &queue,
	bool processed,
	unsigned long long int &index_counter,
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> &regions,
	//std::unordered_map <std::vector<unsigned long long int>, unsigned long long int> &regions_auxiliary_map,
	unsigned long long int level,
	unsigned long long int pocetni_level,
	std::unordered_map <unsigned long long int, bool> &regija_processed,
	std::unordered_map <unsigned long long int, std::string> &regija_zadnjiStrand,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> &unitig_u_kojim_regijama
) {

	std::vector<unsigned long long int> region_unitigs;
	std::vector<unsigned long long int> neighbouring_regions;

	std::vector<unsigned long long int> created_regions;
	std::vector<unsigned long long int> all_created_regions;

	unsigned long long int pom = current_size;

	if (starting_position_in_unitig < K_MER_LENGTH)
		current_size += std::get<0>(unitig_tuple).size() - K_MER_LENGTH;
	else
		current_size += std::get<0>(unitig_tuple).size() - starting_position_in_unitig;


	//std::string region_string;
	std::tuple <bool, bool> temp = std::make_tuple(false, false);

	/*for (auto it = processed_unitigs.begin(); it != processed_unitigs.end(); it++)
		printf("%llu\n", it->first);*/

		//printf("%llu\n", unitig_index);
		//printf("%llu\n", level);



	if (processed_unitigs.find(unitig_index) != processed_unitigs.end()) {
		processed = true;
		temp = processed_unitigs[unitig_index];
	}

	if (strand.compare("+") == 0) {
		processed_unitigs[unitig_index] = std::make_tuple(true, std::get<1>(temp));
	}
	else {
		processed_unitigs[unitig_index] = std::make_tuple(std::get<0>(temp), true);
	}


	//---------------------------------------------------------------------------------------------------

	if (processed) {
		index_counter++;

		regija_zadnjiStrand[index_counter] = strand;
		regija_processed[index_counter] = true;

		//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
		for (unsigned long long int i = 0; i <= level; i++)
			region_unitigs.emplace_back(555);

		//-------------------
		unsigned long long int rSize = 0;

		if (level == 0)
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;
		else if (starting_position_in_unitig < K_MER_LENGTH)
			rSize = pom + std::get<0>(unitig_tuple).size() - K_MER_LENGTH;
		else
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;

		if (rSize > REGION_SIZE)
			rSize = REGION_SIZE;
		//-------------------

		unsigned long long int remove = 0;

		if (current_size > REGION_SIZE)
			remove = current_size - REGION_SIZE;

		//-------------------SUSJEDNOST
		/*for (auto it = regions.begin(); it != regions.end(); it++)
			for (auto unitig : std::get<0>(it->second))
				if (unitig == unitig_index) {
					neighbouring_regions.emplace_back(it->first);
					std::get<3>(it->second).emplace_back(index_counter);
				}*/
		//-------------------\SUSJEDNOST

		region_unitigs[level] = unitig_index;
		regions[index_counter] = std::make_tuple(region_unitigs, 0, std::get<0>(unitig_tuple).size() - 1 - remove, neighbouring_regions, rSize);

		if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
			unitig_u_kojim_regijama[unitig_index][index_counter] = 1;

		all_created_regions.emplace_back(index_counter);

		//if (strand.compare("+") == 0) {
		//	all_created_regions.emplace_back(std::make_tuple(std::get<0>(unitig_tuple), index_counter));
		//}
		//else {
		//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple)), index_counter));
		//}

		return all_created_regions;
	}

	//---------------------------------------------------------------------------------------------------


	if (level != 0) {
		if (std::get<0>(unitig_tuple).size() - K_MER_LENGTH + pom > regions_total_size) {
			unsigned long long int remaining_size = regions_total_size - pom;
			index_counter++;

			regija_zadnjiStrand[index_counter] = strand;

			if (processed == true)
				regija_processed[index_counter] = true;

			//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
			for (unsigned long long int i = 0; i <= level; i++)
				region_unitigs.emplace_back(555);

			//printf("%llu %llu\n", regions_total_size, current_size);
			//printf("%llu\n", remaining_size);

			//if (strand.compare("+") == 0) {
			unsigned long long int temp;

			if (starting_position_in_unitig < K_MER_LENGTH)
				temp = K_MER_LENGTH + remaining_size;
			else
				temp = starting_position_in_unitig + remaining_size;

			region_unitigs[level] = unitig_index;
			regions[index_counter] = std::make_tuple(region_unitigs, 0, temp, neighbouring_regions, REGION_SIZE);

			if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
				unitig_u_kojim_regijama[unitig_index][index_counter] = 1;

			all_created_regions.emplace_back(index_counter);
			//}
			//else {
			//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple).substr(0, remaining_size)), index_counter));
			//}

			if (level == pocetni_level) {

				if (/*pom < regions_total_size / 2 && current_size >= regions_total_size / 2 &&*/ processed == false) {

					processed_unitigs.erase(unitig_index);

					//if (processed_unitigs.find(unitig_index_temp) != processed_unitigs.end())
					//	processed_unitigs.erase(unitig_index_temp);

					queue.emplace_back(std::make_tuple(
						//std::get<0>(regions[region])[level],
						unitig_index,
						starting_position_in_unitig,
						strand,
						index_counter,
						0,
						0,
						remaining_size + starting_position_in_unitig)
					);
				}
			}

			return all_created_regions;
		}
	}

	else {
		if (std::get<0>(unitig_tuple).size() /*+ pom -> pom = 0*/ - starting_position_in_unitig + pom > regions_total_size) {
			unsigned long long int remaining_size = regions_total_size - pom;
			index_counter++;

			regija_zadnjiStrand[index_counter] = strand;

			if (processed == true)
				regija_processed[index_counter] = true;

			//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
			for (unsigned long long int i = 0; i <= level; i++)
				region_unitigs.emplace_back(555);

			//printf("%llu %llu\n", regions_total_size, current_size);
			//printf("%llu\n", remaining_size);

			//if (strand.compare("+") == 0) {
			region_unitigs[level] = unitig_index;
			regions[index_counter] = std::make_tuple(region_unitigs, starting_position_in_unitig, starting_position_in_unitig + remaining_size, neighbouring_regions, REGION_SIZE);

			if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
				unitig_u_kojim_regijama[unitig_index][index_counter] = 1;


			all_created_regions.emplace_back(index_counter);
			//}
			//else {
			//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple).substr(0, remaining_size)), index_counter));
			//}

			if (processed == false) {
				//if (processed_unitigs.find(unitig_index) != processed_unitigs.end())
				processed_unitigs.erase(unitig_index);

				queue.emplace_back(std::make_tuple(
					//std::get<0>(regions[region])[level],
					unitig_index,
					starting_position_in_unitig,
					strand,
					index_counter,
					0,
					0,
					starting_position_in_unitig + remaining_size));
			}

			return all_created_regions;
		}
	}



	bool no_connections = true;

	for (auto &connected_unitig : std::get<1>(unitig_tuple)) {
		if (strand.compare(std::get<0>(connected_unitig)) == 0)
			no_connections = false;
	}
	if (no_connections) {
		index_counter++;

		regija_zadnjiStrand[index_counter] = strand;
		regija_processed[index_counter] = true;

		//samo ubaci nesto tako da popunis vektor sa level elemenata kako bi mogao napraviti vector[level] = nesto
		for (unsigned long long int i = 0; i <= level; i++)
			region_unitigs.emplace_back(555);

		//-------------------
		unsigned long long int rSize = 0;

		if (level == 0)
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;
		else if (starting_position_in_unitig < K_MER_LENGTH)
			rSize = pom + std::get<0>(unitig_tuple).size() - K_MER_LENGTH;
		else
			rSize = pom + std::get<0>(unitig_tuple).size() - starting_position_in_unitig;

		//-------------------

		region_unitigs[level] = unitig_index;
		regions[index_counter] = std::make_tuple(region_unitigs, 0, std::get<0>(unitig_tuple).size(), neighbouring_regions, rSize);

		if (unitig_u_kojim_regijama[unitig_index][index_counter] == 0)
			unitig_u_kojim_regijama[unitig_index][index_counter] = 1;

		all_created_regions.emplace_back(index_counter);

		//if (strand.compare("+") == 0) {
		//	all_created_regions.emplace_back(std::make_tuple(std::get<0>(unitig_tuple), index_counter));
		//}
		//else {
		//	all_created_regions.emplace_back(std::make_tuple(reverse_complement(std::get<0>(unitig_tuple)), index_counter));
		//}

		return all_created_regions;
	}



	for (auto &connected_unitig : std::get<1>(unitig_tuple)) {

		//mislim da ne treba ovaj dio
		/*temp = std::make_tuple (false, false);

		if (processed_unitigs.find(std::get<1>(connected_unitig)) != processed_unitigs.end()) {
					temp = processed_unitigs[connected_unitig];
			}

			if (std::get<2>(connected_unitig).compare("+") == 0) {
					processed_unitigs[connected_unitig] = std::make_tuple (true, std::get<1>(temp));
			}else {
					processed_unitigs[connected_unitig] = std::make_tuple (std::get<0>(temp), true);
			}*/

		if (strand.compare(std::get<0>(connected_unitig)) == 0) {

			created_regions = create_regions_pola_regije(
				unitigs_and_transitions,
				unitigs_and_transitions[std::get<1>(connected_unitig)],
				std::get<1>(connected_unitig),
				0,
				std::get<2>(connected_unitig),
				regions_total_size,
				current_size,
				processed_unitigs,
				queue,
				processed,
				index_counter,
				regions,
				//regions_auxiliary_map,
				level + 1,
				pocetni_level,
				regija_processed,
				regija_zadnjiStrand,
				unitig_u_kojim_regijama
			);

			for (auto &region : created_regions) {

				std::get<0>(regions[region])[level] = unitig_index;

				if (unitig_u_kojim_regijama[unitig_index][region] == 0)
					unitig_u_kojim_regijama[unitig_index][region] = 1;

				if (level == pocetni_level) {
					std::get<1>(regions[region]) = starting_position_in_unitig;

					if (/*pom < regions_total_size / 2 && current_size >= regions_total_size / 2 &&*/regija_processed.find(region) == regija_processed.end()) {
						unsigned long long int unitig_index_temp = std::get<0>(regions[region])[std::get<0>(regions[region]).size() - 1];

						//for (unsigned long long int i = level; i <= std::get<0>(regions[region]).size() - 1; i++)
						//	processed_unitigs.erase(std::get<0>(regions[region])[i]);

						if (processed_unitigs.find(unitig_index_temp) != processed_unitigs.end())
							processed_unitigs.erase(unitig_index_temp);

						queue.emplace_back(std::make_tuple(
							//std::get<0>(regions[region])[level],
							unitig_index_temp,
							starting_position_in_unitig,
							regija_zadnjiStrand[region],
							region,
							std::get<0>(regions[region]).size() - level - 1,
							level,
							std::get<2>(regions[region])));
					}
				}

				all_created_regions.emplace_back(region);

				//printf("%llu %llu %llu %d\n", pom, current_size, regions_total_size / 2, std::get<0>(unitig_tuple).size());


			}
		}

	}

	return all_created_regions;

}


std::vector<unsigned long long int> create_regions(
	std::vector < std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>>> &unitigs_and_transitions,
	std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>> &unitig_tuple,
	unsigned long long int unitig_index,
	unsigned long long int starting_position_in_unitig,
	std::string strand,
	unsigned long long int regions_total_size,
	unsigned long long int current_size,
	std::unordered_map <unsigned long long int, std::tuple<bool, bool>> &processed_unitigs,
	std::vector <std::tuple <unsigned long long int, unsigned long long int, std::string, unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int>> &queue,
	bool processed,
	unsigned long long int &index_counter,
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> &regions,
	//std::unordered_map <std::vector<unsigned long long int>, unsigned long long int> &regions_auxiliary_map,
	unsigned long long int level,
	unsigned long long int verzija,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> &unitig_u_kojim_regijama
) {

	std::unordered_map <unsigned long long int, bool> regija_processed;
	std::unordered_map <unsigned long long int, std::string> regija_zadnjiStrand;

	if (verzija == 0)
		return create_regions_puna_regija(
			unitigs_and_transitions,
			unitig_tuple,
			unitig_index,
			starting_position_in_unitig,
			strand,
			regions_total_size,
			current_size,
			processed_unitigs,
			queue,
			processed,
			index_counter,
			regions,
			level,
			regija_processed, //ako slucajno zadnji unitig nema prijelaza ovdje to biljezimo kako ne bi zadali pravljenje nove regije od polovice
			regija_zadnjiStrand,	//mapa koja biljezi strand zadnjeg unitiga kako bi na polovici to mogli zapisati u queue
			unitig_u_kojim_regijama
		);
	else
		return create_regions_pola_regije(
			unitigs_and_transitions,
			unitig_tuple,
			unitig_index,
			starting_position_in_unitig,
			strand,
			regions_total_size,
			current_size,
			processed_unitigs,
			queue,
			processed,
			index_counter,
			regions,
			level,	//pocetni level (sluzi za pracenje na kojem smo trenutno levelu)
			level,	//pocetni level (ovaj sluzi kako bi zapamtili level od kojeg smo krenuli)
			regija_processed,
			regija_zadnjiStrand,
			unitig_u_kojim_regijama
		);
}

void create_regions_outer_function(
	std::vector <std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>>> &unitigs_and_transitions,
	unsigned long long int region_size,
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> &regions,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> &unitig_u_kojim_regijama
	//std::unordered_map <std::vector<unsigned long long int>, unsigned long long int> &regions_auxiliary_map
) {
	std::unordered_map <unsigned long long int, std::tuple <bool, bool>> processed_unitigs;
	std::vector <std::tuple <unsigned long long int, unsigned long long int, std::string, unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int>> queue; // queue for remaining unitgs that need to be processed
	unsigned long long int i = 0;
	unsigned long long int index_counter = -1;

	std::vector <unsigned long long> created_regions;

	for (auto &unitig_tuple : unitigs_and_transitions) {
		if (processed_unitigs.find(i) != processed_unitigs.end())
			if (std::get<0>(processed_unitigs[i]) == true) {
				i++;
				continue;
			}

		created_regions = create_regions(
			unitigs_and_transitions,
			unitig_tuple,
			i,
			0,
			"+",
			region_size,
			0,
			processed_unitigs,
			queue,
			false,
			index_counter,
			regions,
			//regions_auxiliary_map,
			0,
			0, //verzija puna regija
			unitig_u_kojim_regijama
		);


		while (queue.size() > 0) {
			//0 - pocetni unitig | 1 - pocetna pozicija u prvom unitigu | 2 - strand | 3 - susjedna regija
			//4 - level od kojeg pocinje slijedeca regija | 5 - level na kojem smo bili kad smo dodali regiju na queue
			// 6 - pocetna pozicija u unitigu od kojeg krecemo izradu
			std::tuple <unsigned long long int, unsigned long long int, std::string, unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int> unitig_and_starting_position_and_strand_and_index = queue[0];
			queue.erase(queue.begin());

			created_regions = create_regions(
				unitigs_and_transitions,
				unitigs_and_transitions[std::get<0>(unitig_and_starting_position_and_strand_and_index)],
				std::get<0>(unitig_and_starting_position_and_strand_and_index),
				std::get<6>(unitig_and_starting_position_and_strand_and_index),
				std::get<2>(unitig_and_starting_position_and_strand_and_index),
				region_size,	//std::get<6>(unitig_and_starting_position_and_strand_and_index),//region_size / 2
				region_size / 2,
				processed_unitigs,
				queue,
				false,
				index_counter,
				regions,
				//regions_auxiliary_map,
				std::get<4>(unitig_and_starting_position_and_strand_and_index),
				1, //verzija pola regije
				unitig_u_kojim_regijama
			);


			for (auto &region : created_regions) {

				std::get<3>(regions[std::get<3>(unitig_and_starting_position_and_strand_and_index)]).emplace_back(region);
				std::get<3>(regions[region]).emplace_back(std::get<3>(unitig_and_starting_position_and_strand_and_index));

				unsigned long long int k = 0;
				for (unsigned long long int i = std::get<5>(unitig_and_starting_position_and_strand_and_index); i < std::get<4>(unitig_and_starting_position_and_strand_and_index) + std::get<5>(unitig_and_starting_position_and_strand_and_index); i++) {
					std::get<0>(regions[region])[k] = std::get<0>(regions[std::get<3>(unitig_and_starting_position_and_strand_and_index)])[i];
					k++;

					if (unitig_u_kojim_regijama[std::get<0>(regions[std::get<3>(unitig_and_starting_position_and_strand_and_index)])[i]][region] == 0)
						unitig_u_kojim_regijama[std::get<0>(regions[std::get<3>(unitig_and_starting_position_and_strand_and_index)])[i]][region] = 1;

					//unitig_u_kojim_regijama[std::get<0>(regions[std::get<3>(unitig_and_starting_position_and_strand_and_index)])[i]].emplace_back(region);
				}

				std::get<1>(regions[region]) = std::get<1>(unitig_and_starting_position_and_strand_and_index);

				//ispis za provjeru
				/*printf("%llu |", region);

				for (auto unitig : std::get<0>(regions[region])) {
					printf("%llu ", unitig);
				}
				printf(" |");

				printf("%llu %llu |", std::get<1>(regions[region]), std::get<2>(regions[region]));

				for (auto neighbour : std::get<3>(regions[region])) {
					printf("%llu ", neighbour);
				}
				printf("\n");*/
			}

		}

		i++;
	}

}

struct node {
	node* parent;
	unsigned long long int currentRegion;
	unsigned long long int regionSize;
	unsigned long long int depth;
};


//funkcija koja kreira put od pocetne do zadnje regije
std::vector <unsigned long long int> create_path(
	node* n
) {

	std::vector <unsigned long long int> path;

	if (n->parent == NULL) {
		path.emplace_back(n->currentRegion);

		return path;
	}
	else {
		path = create_path(n->parent);
		path.emplace_back(n->currentRegion);

		return path;
	}
}


//bfs koji trazi put izmedu krajnjih regija koji ce po prilici biti jednak duzini reed-a
std::vector <unsigned long long int> find_optimal_path(
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> &regions,
	unsigned long long int starting_region,
	unsigned long long int ending_region,
	unsigned long long int sequence_length,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> &unitig_u_kojim_regijama,
	unsigned long long int &distance
) {

	std::vector <unsigned long long int> path;
	std::vector <unsigned long long int> optimal_path;
	unsigned long long int best_distance = 1000000;//sequence_length;

	std::unordered_map <unsigned long long int, unsigned long long int> visited_regions;

	if ( sequence_length <= REGION_SIZE || starting_region == ending_region) {
		optimal_path.emplace_back(starting_region);
		return optimal_path;
	}


	node* n;
	node* last_node;
	node* auxil;
	std::vector <node*> open;
	unsigned long long int max_depth = sequence_length / REGION_SIZE + 2;
	bool going_forward = true;

	//inicijalizacija pocetnog cvora
	node startNode;
	startNode.parent = NULL;
	startNode.currentRegion = starting_region;
	startNode.regionSize = std::get<4>(regions[starting_region]);
	startNode.depth = 0;

	open.emplace_back(&startNode);

	while (open.size() > 0) {
		n = open.back();
		open.pop_back();

		if (!going_forward) {
			auxil = last_node;
			while ( auxil -> parent != n -> parent) {
				visited_regions.erase(auxil->parent->currentRegion);
				auxil = auxil -> parent;
			}
		}

		last_node = n;

		/*if (n->regionSize > 10000 * sequence_length) {
			printf("Uf\n");
			break;
		}*/

		/*if (n->currentRegion == ending_region) {
			if ((unsigned long long int)abs((long long int) sequence_length - (long long int) n->regionSize) <= best_distance) {
				path = create_path(n);
				optimal_path = path;
				best_distance = (unsigned long long int) abs((long long int) sequence_length - (long long int)n->regionSize);
			}
			else
				break;
		}*/

		unsigned long long int last_unitig_in_region = std::get<0>(regions[n->currentRegion])[std::get<0>(regions[n->currentRegion]).size() - 1];
		unsigned long long int first_unitig_in_region = std::get<0>(regions[n->currentRegion])[0];

		going_forward = false;

		if (n->depth < max_depth) {

		for (auto next_region : unitig_u_kojim_regijama[last_unitig_in_region]) {
			if (visited_regions[next_region.first] == 1)
				continue;

			visited_regions[next_region.first] = 1;

			if (next_region.first > n->currentRegion)
                        if (n->parent == nullptr or (next_region.first != n->parent->currentRegion and n->currentRegion != next_region.first)) {

				if (n->currentRegion == ending_region) {
                        	if ((unsigned long long int)abs((long long int) sequence_length - (long long int) n->regionSize) <= best_distance) {
                                	path = create_path(n);
                                	optimal_path = path;
                                	best_distance = (unsigned long long int) abs((long long int) sequence_length - (long long int)n->regionSize);
                        	}
                        	else
                                	break;
                		}


                                node* nextNode = new node();
                                nextNode->parent = n;
                                nextNode->currentRegion = next_region.first;
                                nextNode->regionSize = n->regionSize + std::get<4>(regions[next_region.first]) / 2;
				nextNode->depth = n->depth + 1;

                                open.emplace_back(nextNode);
				going_forward = true;
                        }
                }

		 for (auto next_region : unitig_u_kojim_regijama[first_unitig_in_region]) {
                        if (visited_regions[next_region.first] == 1)
                                continue;

                        visited_regions[next_region.first] = 1;

			if (next_region.first < n->currentRegion)
                        if (n->parent == nullptr or (next_region.first != n->parent->currentRegion and n->currentRegion != next_region.first)) {

				if (n->currentRegion == ending_region) {
                        	if ((unsigned long long int)abs((long long int) sequence_length - (long long int) n->regionSize) <= best_distance) {
                                	path = create_path(n);
                                	optimal_path = path;
                                	best_distance = (unsigned long long int) abs((long long int) sequence_length - (long long int)n->regionSize);
                        	}
                        	else
                                	break;
                		}

				node* nextNode = new node();
                                nextNode->parent = n;
                                nextNode->currentRegion = next_region.first;
                                nextNode->regionSize = n->regionSize + std::get<4>(regions[next_region.first]) / 2;
				nextNode->depth = n->depth + 1;

                                open.emplace_back(nextNode);
				going_forward = true;
                        }
                }

		}

		/*for (auto next_region : std::get<3>(regions[n->currentRegion])) {
			if (n->parent == nullptr or next_region != n->parent->currentRegion) {
				node* nextNode = new node();
				nextNode->parent = n;
				nextNode->currentRegion = next_region;
				nextNode->regionSize = n->regionSize + std::get<4>(regions[next_region]);

				open.emplace_back(nextNode);
			}
		}*/
	}

	std::vector <unsigned long long int> empty;
	distance = best_distance;
	//if ((float)best_distance < ACCURACY_OF_MAPPING * sequence_length)
		return optimal_path;
	//else
	//	return empty;
}


//globalna varijabla
unsigned long long int number_of_discarded_sequences = 0;

//pomocna funkcija kako bi mogli napraviti efikasniji paralelizam

void map_sequence_to_reference(
	//std::unordered_map <unsigned long long int, unsigned long long int> &reference_genome_minimizers_positions_in_vector,
        //std::vector < std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>>> &reference_genome_minimizers,
	std::vector <std::unique_ptr <SequenceFormat>> &sequences,
	//std::vector <std::unique_ptr <SequenceFormat>> &reference_genome,
	//std::vector <std::tuple <std::unordered_map <unsigned long long int, unsigned long long int>, std::vector <std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>>>>> &unitig_minimizers,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> hash_minimizer_region_number_to_add,
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> regions,
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> unitig_u_kojim_regijama,
	int k,
	int window_size,
	unsigned int start,
	unsigned int end,
	std::unordered_map <unsigned long long int, std::vector <unsigned long long int>> &for_mapping
) {

	//unsigned int no_matches = 0;
	//unsigned int no_match_groups_original = 0;
	//unsigned int no_match_groups_reverse_complement = 0;


	for (unsigned int i = start; i < end; i++) {
		//printf("Bok! Evo me na iteraciji %u\n", i);

		//timestamp_t t0 = get_timestamp();

		std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> sequence_start_minimizer_index;
		std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> sequence_end_minimizer_index;

        	/*std::vector <std::tuple <unsigned long long int, unsigned long long int>> matches_same_strand_start;
		std::vector <std::tuple <unsigned long long int, unsigned long long int>> matches_different_strand_start;

		std::vector <std::tuple <unsigned long long int, unsigned long long int>> matches_same_strand_end;
                std::vector <std::tuple <unsigned long long int, unsigned long long int>> matches_different_strand_end;
		*/

		std::unordered_map <unsigned long long int, unsigned long long int> matches_per_region_start;
		std::unordered_map <unsigned long long int, unsigned long long int> matches_per_region_end;

		unsigned long long int best_region_start = regions.begin()->first; //indeks regije koja ima najvise matcheva sa pocetkom sekvence
		unsigned long long int best_region_end = regions.begin()->first; //indeks regije koja ima najvise matcheva sa krajem sekvence


		//-------4x version-------

		unsigned long long int second_best_region_start = regions.begin()->first; //indeks regije koja ima najvise matcheva sa pocetkom sekvence, ali zanemarimo li najbolju
                unsigned long long int second_best_region_end = regions.begin()->first; //indeks regije koja ima najvise matcheva sa krajem sekvence, ali zanemarimo li najbolju

		//-------4x version-------

		std::string sequence = sequences[i] -> sequence;
		unsigned long long int sequence_size = sequence.size();

		if(sequence_size < SEQUENCE_ENDS_SIZE * 2) {

			//minimizatori sekvence

			//timestamp_t t3 = get_timestamp();

			std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>> sequence_minimizer_index;
			sequence_minimizer_index = white::minimizers(sequence.c_str(), sequence_size, k, window_size);

			//timestamp_t t4 = get_timestamp();

                        //time_spent_on_creating_minimizers += (t4 - t3) / 1000000.0L;

			//najbolja regija
			unsigned long long int best_region = regions.begin()->first;

			//uzimamo samo najbolju regiju na koju se mapirala sekvenca

			//matchiranje

			//timestamp_t t5 = get_timestamp();
                	for (auto minimizer_tuple : sequence_minimizer_index)
                        	for (auto region_number_to_add : hash_minimizer_region_number_to_add[std::get<0>(minimizer_tuple)])
                        	        matches_per_region_start[region_number_to_add.first] += region_number_to_add.second;

			//timestamp_t t6 = get_timestamp();

                	//time_spent_on_match += (t6 - t5) / 1000000.0L;

			//najbolja regija

			//timestamp_t t7 = get_timestamp();

			for ( auto it = regions.begin(); it != regions.end(); it++ )
				if (matches_per_region_start[best_region] < matches_per_region_start[it->first])
                                        best_region = it->first;

			//timestamp_t t8 = get_timestamp();

                        //time_spent_on_sorting += (t8 - t7) / 1000000.0L;

			//spremimo najbolju regiju u for mapping

			std::vector <unsigned long long int> temp;
			temp.emplace_back(best_region);

			for_mapping[i] = temp;

			//timestamp_t t1 = get_timestamp();

    			//time_spent_on_matching += (t1 - t0) / 1000000.0L;

			continue;
		}

		std::string sequence_start = sequence.substr (0, SEQUENCE_ENDS_SIZE);
		std::string sequence_end = sequence.substr (sequence_size - SEQUENCE_ENDS_SIZE, SEQUENCE_ENDS_SIZE);

		//timestamp_t t3 = get_timestamp();

		//za sekvencu racunamo minimizatore
		sequence_start_minimizer_index = white::minimizers(sequence_start.c_str(), sequence_start.size(), k, window_size);
		sequence_end_minimizer_index = white::minimizers(sequence_end.c_str(), sequence_end.size(), k, window_size);

		//timestamp_t t4 = get_timestamp();

                //time_spent_on_creating_minimizers += (t4 - t3) / 1000000.0L;

		//printf("Pogotci i najbolja regija!\n");
		//racunamo najbolju regiju za pocetak i kraj sekvence

		//matchiranje

		//timestamp_t t5 = get_timestamp();

		for (auto minimizer_tuple : sequence_start_minimizer_index)
			for (auto region_number_to_add : hash_minimizer_region_number_to_add[std::get<0>(minimizer_tuple)])
				matches_per_region_start[region_number_to_add.first] += region_number_to_add.second;

		for (auto minimizer_tuple : sequence_end_minimizer_index)
                        for (auto region_number_to_add : hash_minimizer_region_number_to_add[std::get<0>(minimizer_tuple)])
                                matches_per_region_end[region_number_to_add.first] += region_number_to_add.second;

		//timestamp_t t6 = get_timestamp();

                //time_spent_on_match += (t6 - t5) / 1000000.0L;

		/*for ( auto it = unitig_u_kojim_regijama.begin(); it != unitig_u_kojim_regijama.end(); it++ ) {

			minimizer_matches(
                                        std::get<0>(unitig_minimizers[it->first]),
                                        std::get<1>(unitig_minimizers[it->first]),
                                        sequence_start_minimizer_index,
                                        matches_same_strand_start,
                                        matches_different_strand_start
                                );


                        minimizer_matches(
                                        std::get<0>(unitig_minimizers[it->first]),
                                        std::get<1>(unitig_minimizers[it->first]),
                                        sequence_end_minimizer_index,
                                        matches_same_strand_end,
                                        matches_different_strand_end
                                );


			for ( auto region_iterator = it->second.begin(); region_iterator != it->second.end(); region_iterator++ ) {
                                matches_per_region_start[region_iterator->first] += matches_same_strand_start.size() + matches_different_strand_start.size();
				matches_per_region_end[region_iterator->first] += matches_same_strand_end.size() + matches_different_strand_end.size();
			}

                        matches_same_strand_start.clear();
                        matches_different_strand_start.clear();
			matches_same_strand_end.clear();
                        matches_different_strand_end.clear();

		}*/

		for ( auto it = regions.begin(); it != regions.end(); it++ ) {

		//----------
		//radimo matcheve minimizatora (i za minimizatore koji dolaze iz originalnih sekvenci i reverznih komplementa)

			/*for (auto unitig : std::get<0>(it->second)) {
				minimizer_matches(
					std::get<0>(unitig_minimizers[unitig]),
					std::get<1>(unitig_minimizers[unitig]),
					sequence_start_minimizer_index,
					matches_same_strand,
					matches_different_strand
				);


				matches_per_region_start[it->first] += matches_same_strand.size() + matches_different_strand.size();

				matches_same_strand.clear();
				matches_different_strand.clear();


				minimizer_matches(
					std::get<0>(unitig_minimizers[unitig]),
                                        std::get<1>(unitig_minimizers[unitig]),
                                	sequence_end_minimizer_index,
                                	matches_same_strand,
                                	matches_different_strand
                        	);

				matches_per_region_end[it->first] += matches_same_strand.size() + matches_different_strand.size();


				matches_same_strand.clear();
                                matches_different_strand.clear();
			}*/

			//printf("M: %llu, %llu\n", matches_per_region_start[it->first], matches_per_region_end[it->first]);



			//poastavimo najbolju regiju
			/*if (matches_per_region_start[best_region_start] < matches_per_region_start[it->first])
                                        best_region_start = it->first;


			if (matches_per_region_end[best_region_end] < matches_per_region_end[it->first])
                                        best_region_end = it->first;*/


			//timestamp_t t7 = get_timestamp();

			//poastavimo najbolju regiju
			if (matches_per_region_start[best_region_start] < matches_per_region_start[it->first])
                                        best_region_start = it->first;

			//-------4x version-------

			else if (matches_per_region_start[second_best_region_start] <= matches_per_region_start[it->first])
				second_best_region_start = it->first;

			else if (best_region_start == second_best_region_start)
				second_best_region_start = it->first;

			//-------4x version-------


			if (matches_per_region_end[best_region_end] < matches_per_region_end[it->first])
                                        best_region_end = it->first;

			//-------4x version-------

			else if (matches_per_region_end[second_best_region_end] <= matches_per_region_end[it->first])
                                second_best_region_end = it->first;

                        else if (best_region_end == second_best_region_end)
                                second_best_region_end = it->first;

			//-------4x version-------

			//timestamp_t t8 = get_timestamp();

                        //time_spent_on_match += (t8 - t7) / 1000000.0L;

		}

		//timestamp_t t1 = get_timestamp();

    		//time_spent_on_matching += (t1 - t0) / 1000000.0L;

		//---------trazimo regije koje cemo iskoristiti za alignment

		//algoritam za naci sve moguce regije za poravnanje
		//double minimizer_percentage = 2 / (double)(window_size + 1);
		//double expected_number_of_minimizers = (double) (sequence_start.size() - k + 1) * minimizer_percentage;

		//t0 = get_timestamp();

		//-------4x version-------

		std::vector <unsigned long long int> path1;
		std::vector <unsigned long long int> path2;
		std::vector <unsigned long long int> path3;
		std::vector <unsigned long long int> path4;

		//printf("%llu, %llu, %llu, %llu, %llu\n", matches_per_region_start[best_region_start], matches_per_region_end[best_region_end], best_region_start, best_region_end, sequence_size);

		//-------4x version-------

		//printf("%llu, %.2f, %.2f\n", matches_per_region_start[best_region_start], expected_number_of_minimizers, (double) matches_per_region_start[best_region_start] / expected_number_of_minimizers);

		//std::vector <std::vector <unsigned long long int>> paths_for_alignment;

		/*std::vector <unsigned long long int> temp;

		for_mapping[i] = temp;

		for (auto it = matches_per_region_start.begin(); it != matches_per_region_start.end(); it++ ) {
			if (it -> second >= expected_number_of_minimizers * 0.01)
				for (auto it2 = matches_per_region_end.begin(); it2 != matches_per_region_end.end(); it2++) {
					if (it2 -> second >= expected_number_of_minimizers * 0.01) {
						path = find_optimal_path (regions, it -> first, it2 -> first, sequence_size);

						if (path.size() != 0)
							for (auto region : path)
								for_mapping[i].emplace_back(region);
					}
				}
		}

		if(for_mapping[i].size() == 0) {
			for_mapping.erase(i);
			number_of_discarded_sequences++;
			continue;
		}*/

		//printf("Put!\n");
		/*//algoritam za naci samo najbolje poravnanje
		//if (matches_per_region_start[best_region_start] >= expected_number_of_minimizers * 0.9 && matches_per_region_end[best_region_end] >= expected_number_of_minimizers * 0.9) {
			printf("Unutra!\n");
			path = find_optimal_path(regions, best_region_start, best_region_end, sequence_size, unitig_u_kojim_regijama);

			if (path.size() != 0) {
				for_mapping[i] = path;
			}
			else {
				printf("Continue3!\n");
				number_of_discarded_sequences++;
                        	continue;
			}

		//}
		//else {
			//printf("Continue2!\n");
			//number_of_discarded_sequences++;
			//continue;
		//}
		printf("\n");
			for (auto region : path)
				printf("%llu ", region);

		printf("\n");

		printf("Van!\n");*/


				//printf("Put!\n");
		//algoritam za naci samo najbolje poravnanje
		//if (matches_per_region_start[best_region_start] >= expected_number_of_minimizers * 0.9 && matches_per_region_end[best_region_end] >= expected_number_of_minimizers * 0.9
			//-------4x version-------

			//&& matches_per_region_start[second_best_region_start] >= expected_number_of_minimizers * 0.8 && matches_per_region_end[second_best_region_end] >= expected_number_of_minimizers * 0.8) {
			//printf("Unutra!\n");
			//printf("%llu %llu %llu %llu\n", best_region_start, best_region_end, second_best_region_start, second_best_region_end);
			std::vector <unsigned long long int> temp;
			for_mapping[i] = temp;
			unsigned long long int best_distance;
			unsigned long long int distance;
			unsigned long long int which_path = 1;

			//start_prva && end_prva
			path1 = find_optimal_path(regions, best_region_start, best_region_end, sequence_size, unitig_u_kojim_regijama, best_distance);
			//printf("Bok\n");
			//start_prva && end_druga

			path2 = find_optimal_path(regions, best_region_start, second_best_region_end, sequence_size, unitig_u_kojim_regijama, distance);

			if (distance < best_distance) {
				best_distance = distance;
				which_path = 2;
			}
			//printf("Bok\n");
			//start_druga && end_prva

			path3 = find_optimal_path(regions, second_best_region_start, best_region_end, sequence_size, unitig_u_kojim_regijama, distance);

                        if (distance < best_distance) {
                                best_distance = distance;
                                which_path = 3;
                        }
			//printf("Bok\n");
			//start_druga && end_druga

			path4 = find_optimal_path(regions, second_best_region_start, second_best_region_end, sequence_size, unitig_u_kojim_regijama, distance);

                        if (distance < best_distance) {
                                best_distance = distance;
                                which_path = 4;
                        }

			//printf("Bok\n");
			if (which_path == 1) {
                        	if (path1.size() != 0) {
                                	for (auto region : path1)
                                        	for_mapping[i].emplace_back(region);
                        	}
			}
			else if (which_path == 2) {
                                if (path2.size() != 0) {
                                        for (auto region : path2)
                                                for_mapping[i].emplace_back(region);
                                }
			}
			else if (which_path == 3) {
                                if (path3.size() != 0) {
                                        for (auto region : path3)
                                                for_mapping[i].emplace_back(region);
                                }
			}
			else if (which_path == 4) {
                                if (path4.size() != 0) {
                                        for (auto region : path4)
                                                for_mapping[i].emplace_back(region);
                                }
			}


			//t1 = get_timestamp();

    			//time_spent_on_dfs += (t1 - t0) / 1000000.0L;

			if (for_mapping[i].size() == 0) {
                                //printf("Continue3!\n");
				for_mapping.erase(i);
                                number_of_discarded_sequences++;
                                continue;
                        }



			//-------4x version-------


	}
}

//-----------------------------------------------------DO OVDJE-----------------------------------------------------------------------------------------------------------------------------

/*		//za svaku match grupu racunamo LIS, zatim radimo alignment dobivenih regija i printamo PAF
		if (!match_groups_same_strand.empty())
			for (auto match_group : match_groups_same_strand)
				ispis_za_pojedinu_dretvu.emplace_back (LIS_alignment_PAF (match_group, sequences[i], reference_genome[0], '+', match, mismatch, gap, k, print_cigar));
		//else
		//	no_match_groups_original++;

		if(!match_groups_different_strand.empty())
			for (auto match_group : match_groups_different_strand)
+				ispis_za_pojedinu_dretvu.emplace_back (LIS_alignment_PAF (match_group, sequences[i], reference_genome[0], '-', match, mismatch, gap, k, print_cigar));
		//else
*/		//	no_match_groups_reverse_complement++;
/*	}
		//std::cout << no_matches << std::endl;
		//std::cout << no_match_groups_original << " " << no_match_groups_reverse_complement << std::endl;
}*/


//-------------------------------------------------------------------------------------------------------------------------

//funkcija za zapisivanje jednog rezultantnog grafova u .fasta file

void create_fasta_graph_file(
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> &regions,
	std::vector <std::string> &unitig_strings_for_creating_fasta_file,
	std::unordered_map <unsigned long long int, std::vector <unsigned long long int>> &single_thread_mapping,
	unsigned long long int sequence
) {
	std::unordered_map <unsigned long long int, unsigned long long int> already_added_unitigs;

	std::string graph_file_name = "output/tmp/graph_for_sequence_" + std::to_string(sequence) + ".fa";
	std::ofstream outfile (graph_file_name);

	for (auto region : single_thread_mapping[sequence])
		for (auto unitig : std::get<0>(regions[region]))
			if (already_added_unitigs[unitig] == 0) {
				outfile << unitig_strings_for_creating_fasta_file[unitig];
				already_added_unitigs[unitig] = 1;
			}

	outfile.close();
}

//funkcija za zapisivanje svih rezultantnih grafova u .fasta filove

void create_fasta_graph_files(
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> &regions,
	std::vector <std::string> &unitig_strings_for_creating_fasta_file,
	std::vector <std::unordered_map <unsigned long long int, std::vector <unsigned long long int>>> &for_mapping
) {

	/*const int dir_err = mkdir("output", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	if (dir_err == -1)
	{
	        printf("Error creating directory!\n");
        	exit (1);
	}*/

	struct stat buffer;
  	if (stat ("output", &buffer) != 0)
		system("mkdir output");

	if (stat ("output", &buffer) != 0)
		system("mkdir output/tmp");

	for (auto single_thread_mapping : for_mapping)
		for (auto it = single_thread_mapping.begin(); it != single_thread_mapping.end(); it++)
			create_fasta_graph_file(
				regions,
				unitig_strings_for_creating_fasta_file,
				single_thread_mapping,
				it->first
			);

}

//---------------------------------------------------FUNKCIJA_RAZDVAJANJE_SEKVENCI_ZA_MAPIRANJE_U_ZASEBEN_FILEOVE----------------------------------------------

std::vector <std::string> seperate_sequences_for_mapping_into_seperate_files (std::string sequences_file_path) {

	std::ifstream graph_file(sequences_file_path);
	std::string line;
	std::vector <std::string> file_paths;

	struct stat buffer;

        if (stat ("output/tmp3", &buffer) != 0)
		system("mkdir output/tmp3");

        if (graph_file.is_open()) {
		int i = 0;
                while (getline(graph_file, line)) {
			std::stringstream ss;
			ss << "output/tmp3/sequence_" << i << ".fq";
			i++;

			file_paths.emplace_back(ss.str());
			std::ofstream outfile (ss.str());

			for (int j = 0; j < 2; j++) {
				outfile << line;
				getline(graph_file, line);
			}
			outfile << line;
		}

	} else {
                std::cout << "Unable to open sequences file!\n";
                exit(1);
        }

	return file_paths;

}



//-----------------------------------------------------------FUNKCIJA_ZA_POZIV_GRAPH_ALIGNERA------------------------------------------------------------------

void run_GraphAligner (std::string path_for_GraphAligner, std::vector<std::string> &sequences_file_path, unsigned long long int number_of_sequences) {

	struct stat buffer;

	for (unsigned long long int i = 0; i < number_of_sequences; i++) {

		std::stringstream ss;
		ss << path_for_GraphAligner << " -g output/tmp2/graph_for_sequences_" << i << ".gfa" << " -f " << sequences_file_path[i] << " -a output/alignment_for_sequence_" << i << ".json";
		std::string new_file_name = ss.str();

		if (stat (ss.str().c_str(), &buffer) == 0) {

			//std::cout << ss.str();
        		system(ss.str().c_str());

		}


	}

}

//-----------------------------------------------------------FUNKCIJE_ZA_PREBACIT_FA_FORMAT_U_GFA--------------------------------------------------------------

void fa_to_gfa_single_file (std::string path_for_bcalm_convert_fa_to_gfa, std::string file_to_convert, unsigned long long int k_mer_size) {

        std::string new_file_name;
        std::string string_to_find = ".fa";

	struct stat buffer;
        if (stat ("output/tmp2", &buffer) != 0)
		system("mkdir output/tmp2");

        new_file_name = "output/tmp2" + file_to_convert.substr( 10, file_to_convert.find(string_to_find) ) + ".gfa";

        // command line arguments
        std::stringstream ss;
        ss << path_for_bcalm_convert_fa_to_gfa << " " << file_to_convert << " " << new_file_name << " " << k_mer_size;
        //std::string arguments = ss.str();

        //std::cout << ss.str();
        system(ss.str().c_str());

}

void fa_to_gfa (std::string path_for_bcalm_convert_fa_to_gfa, unsigned long long int number_of_sequences) {

	for (unsigned long long int i = 0; i < number_of_sequences; i++) {
		std::stringstream graph_file_name;
		graph_file_name << "output/tmp/graph_for_sequence_" << i << ".fa";
		struct stat buffer;

		if (stat (graph_file_name.str().c_str(), &buffer) == 0)
			fa_to_gfa_single_file (path_for_bcalm_convert_fa_to_gfa, graph_file_name.str(), K_MER_LENGTH);
	}

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------

int main (int argc, char* argv[])
{
	int option = 0;
	std::vector<std::unique_ptr<SequenceFormat>> sequences;
	std::vector<std::unique_ptr<SequenceFormat>> reference_genome;

	unsigned int abundance_min = ABUNDANCE_MIN;

	unsigned int k = K_MER_LENGTH;
	unsigned int window_size = WINDOW_SIZE;

	unsigned int number_of_threads = NUMBER_OF_THREADS;

        while ((option = getopt_long(argc, argv, "hvk:w:f:a:t:c", long_options, NULL)) != -1) {
                switch (option) {
                        case 'h':
				help();
                                exit(0);

                        case 'v':
				version();
                                exit(0);

			case 'k':
				k = atoi (optarg);
				break;

			case 'w':
				window_size = atoi (optarg);
				break;

			case 'a':
				abundance_min = std::atoi (optarg);
				break;

			case 't':
				number_of_threads = std::atoi (optarg);
				break;

                        default:
                                fprintf (stderr, "\n--Error:\n\tUnsupported option. Usage: %s [-h] [--help] [-v] [--version] [-m ARG] [-s ARG] [-g ARG] [-k ARG] [-w ARG] [-f ARG] "
						"sequence_file reference_genome_file match mismatch gap\n\n", argv[0]);
                                exit(0);

                }
	}

	if (argc - optind == 2)
	{
		if (check_format (argv[optind], fasta_formats))
			sequences = parse_file<bioparser::FastaParser> (std::string (argv[optind]));

		else if (check_format (argv[optind], fastq_formats))
                        sequences = parse_file<bioparser::FastqParser> (std::string (argv[optind]));

		else
			fprintf (stderr, "\nError:													\n"
					 "	Unsuported format. File containing reeds (1st file) needs to have one of the following extensions:	\n"
                                     	 "	-> .fasta             -> .fastq										\n"
                                     	 "	-> .fa                -> .fq										\n"
                                     	 "	-> .fasta.gz          -> .fastq.gz									\n"
                                     	 "	-> .fa.gz             -> .fq.gz										\n\n");

		/*if (check_format (argv[++optind], fasta_formats)) {
                        reference_genome = parse_file<bioparser::FastaParser> (std::string (argv[optind]));
			reference_genome_size = reference_genome[0] -> sequence.size();
		}
		else

			fprintf (stderr, "\n --Error:														\n"
					 "	Unsuported format. File containing reference genome (2nd file) needs to have one of the following extensions:	\n"
                                     	 "	-> .fasta													\n"
                                     	 "	-> .fa														\n"
                                     	 "	-> .fasta.gz													\n"
                                     	 "	-> .fa.gz													\n\n");
		*/
	}

	 else
                        fprintf (stderr, "\n--Error:\n"
					 "	Usage: white_mapper [OPTIONS] reeds_file reference_genome_file				\n\n"

					 "	Number of nonoption arguments must be 2:						\n"
                                         "		1st file contains reeds in FASTA or FASTQ format				\n"
                                         "		2nd file contains a reference genome in FASTA fromat				\n\n");



//-----------------------------------------------UCITAVANJE_PATHOVA-------------------------------------------------------------------------------------------------------------------

/*	std::string path_for_GraphAligner;
        std::string path_for_bcalm;
        std::string path_for_bcalm_convert_to_gfa;
        std::string path_for_pbsim;
        std::string path_for_pbsim_model;

	std::string line;
	std::string path_for_what;
	std::string path;

        std::ifstream path_file("./paths.txt");

        if (path_file.is_open()) {

                std::string delimiter = " ";

                while (getline(path_file, line)) {

                        if (line.find(delimiter) != std::string::npos) {

                                path_for_what = line.substr(0, line.find(delimiter));
                                path = line.substr(line.find(delimiter) + 1, line.size() - line.find(delimiter));

                                if (path_for_what.compare("aligner_bin_path:") == 0)
                                        path_for_GraphAligner = path;

                                else if (path_for_what.compare("bcalm_bin_path:") == 0)
                                        path_for_bcalm = path;

                                else if (path_for_what.compare("bcalm_convert_path:") == 0)
                                        path_for_bcalm_convert_to_gfa = path;

                                else if (path_for_what.compare("pbsim_bin_path:") == 0)
                                        path_for_pbsim = path;

                                else if (path_for_what.compare("pbsim_config_path") == 0)
                                        path_for_pbsim_model = path;
                        }
                }

                path_file.close();

        } else {
                std::cout << "Unable to open path file!\n";
                exit(1);
        }

	//printf("\nPaths:\n%s\n%s\n%s\n%s\n%s\n", path_for_GraphAligner.c_str(), path_for_bcalm.c_str(), path_for_bcalm_convert_to_gfa.c_str(), path_for_pbsim.c_str(), path_for_pbsim_model.c_str());
*/
//-----------------------------------------------KONSTRUKCIJA_GRAFA-------------------------------------------------------------------------------------------------------------------

	// program name
	/*std::string name = std::string("~/GitHub/BSc-thesis/vendor/bcalm/build/bcalm");

	// command line arguments
	std::stringstream ss;
	ss << name << " -in " << argv[optind] << " -kmer-size " << k << " -abundance-min " << abundance_min;
	std::string arguments = ss.str();

	std::cout << ss.str();
	system(ss.str().c_str());*/

//-----------------------------------------------OBRADA_DOBIVENOG_GRAPHA--------------------------------------------------------------------------------------------------------------

	//std::string line; --> gore deklarirana
	std::ifstream graph_file(argv[++optind]);

	unsigned long long int position;
	std::vector <std::tuple <std::string, std::vector <std::tuple <std::string, unsigned long long int, std::string>>>> unitigs_and_transitions;

	std::string transition;

	//--
	std::vector <std::string> unitig_strings_for_creating_fasta_file;
	std::string unitig_string;
	//--

	if (graph_file.is_open()) {
		while (getline(graph_file, line)) {

			//--
			unitig_string = line + "\n";
			//--

			std::vector <std::tuple <std::string, unsigned long long int, std::string>> transitions;
			std::string delimiter = " ";
			//printf ("\n%s\n", line.c_str());
			//printf ("%d\n", line.find(delimiter));
			for (unsigned int k = 0; k < 4; k++) {
				line.erase(0, line.find(delimiter) + delimiter.length());
			}


			//printf ("%s\n", line.c_str());
			while (line.find(delimiter) != std::string::npos) {
				while (line.find(delimiter) == 0 || line.compare("\r") == 0) {
					line.erase(0, 1);
				}

				if (line.size() == 0) break;
				//printf("%d\n", line.size());
				//printf ("%s\n", line.c_str());

				std::string delimiter2 = ":";

				if (line.find(delimiter) != std::string::npos) {
					transition = line.substr(0, line.find(delimiter));
					line.erase(0, line.find(delimiter));
				}
				else {
					transition = line;
					line == std::string("");
				}
				//printf ("%s\n", transition.c_str());
				std::vector <std::string> temp;
				std::tuple <std::string, unsigned long long int, std::string> trans;

				while ((position = transition.find(delimiter2)) != std::string::npos) {
					transition.erase(0, position + delimiter2.length());

					if (transition.find(delimiter2) != std::string::npos) {
						temp.emplace_back(transition.substr(0, transition.find(delimiter2)));
					}
					else {
						temp.emplace_back(transition);
					}
				}

				trans = std::make_tuple(temp[0], std::atoi(temp[1].c_str()), temp[2]);
				transitions.emplace_back(trans);
			}

			getline(graph_file, line);
			unitigs_and_transitions.emplace_back(std::make_tuple(line, transitions));

			//--
			unitig_string += line + "\n";
			//printf("%s", unitig_string.c_str());
			unitig_strings_for_creating_fasta_file.emplace_back(unitig_string);
			//--
		}

		graph_file.close();
	}
	else {
		std::cout << "Unable to open file!\n";
		exit(1);
	}

	//ispis za provjeru
	/*for (auto tuple : unitigs_and_transitions) {
			printf ("%s\n", std::get<0>(tuple).c_str());
			for (auto tuple2 : std::get<1>(tuple)) {
				printf ("%s %llu %s ", std::get<0>(tuple2).c_str(), std::get<1>(tuple2), std::get<2>(tuple2).c_str());
			}
			std::cout << std::endl;
	}*/

//-----------------------------------------------ZADNJI_ZADATAK-----------------------------------------------------------------------------------------------------------------------

	//std::unordered_map <unsigned long long int, unsigned long long int> reference_genome_minimizers_positions_in_vector;
        //std::vector < std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>>> reference_genome_minimizers;

	/*get_reference_genome_minimizers(
		reference_genome[0],
		k,
		window_size,
		f,
		reference_genome_minimizers_positions_in_vector,
		reference_genome_minimizers
	);*/

	//-------------

	//stvorimo minimizatore unitigova
	//std::vector <std::tuple <std::unordered_map <unsigned long long int, unsigned long long int>, std::vector <std::vector <std::tuple <unsigned long long int, unsigned long long int, bool>>>>> unitig_minimizers;

	//unitig_minimizers = create_unitig_minimizers(unitigs_and_transitions, k, window_size);

	//-------------

	//timestamp_t t0 = get_timestamp();

	//stvorimo regije
	std::unordered_map <unsigned long long int, std::tuple <std::vector<unsigned long long int>, unsigned long long int, unsigned long long int, std::vector<unsigned long long int>, unsigned long long int>> regions;
	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> unitig_u_kojim_regijama;

	create_regions_outer_function(
	         unitigs_and_transitions,
	         REGION_SIZE,
                 regions,
		 unitig_u_kojim_regijama
        );

	//-------------

	//brisanje premalih regija
	std::vector <unsigned long long int> regions_to_delete;

	for (auto it = regions.begin(); it != regions.end(); it++)
		if (std::get<4>(it->second) < REGION_SIZE / 2)
			regions_to_delete.emplace_back(it->first);

	for (auto region : regions_to_delete) {
		for (auto unitig : std::get<0>(regions[region]))
			unitig_u_kojim_regijama[unitig].erase(region);

		regions.erase(region);
	}


	//-------------

	//BRISANJE MALIH REGIJA KOJE NEMAJU SUSJEDA ---- OVO IDE UZ SUSJEDNOST
	/*std::vector <unsigned long long int> regions_to_delete;

	for (auto it = regions.begin(); it != regions.end(); it++)
		if (std::get<4>(it->second) < REGION_SIZE / 2 && std::get<3>(it->second).size() == 0)
			regions_to_delete.emplace_back(it->first);

	for (auto region : regions_to_delete)
		regions.erase(region);
*/
	//-------------

	//stvaramo hash za matchanje

	std::unordered_map <unsigned long long int, std::unordered_map <unsigned long long int, unsigned long long int>> hash_minimizer_region_number_to_add;
	hash_minimizer_region_number_to_add = create_hash_for_matching(unitigs_and_transitions, unitig_u_kojim_regijama, k, window_size);

	//timestamp_t t1 = get_timestamp();

    	//time_spent_on_creating_regions = (t1 - t0) / 1000000.0L;

	//ispis za provjeru

	/*for (auto it : hash_minimizer_region_number_to_add) {
		printf("\nMinimizator: %llu\n", it.first);

		for (auto it2 : it.second) {
			printf("\nRegija: %llu | Dodati: %llu", it2.first, it2.second);
		}

		printf("\n");
	 }*/

	//exit(0);

	/*for (auto it = unitig_u_kojim_regijama.begin(); it != unitig_u_kojim_regijama.end(); it++) {
		printf("%llu |", it->first);

		for (auto regija : it->second)
			printf(" %llu", regija);

		printf("\n");

	}

	//ispis za provjeru
	for (auto it = regions.begin(); it != regions.end(); it++) {
		printf("%llu |", it->first);

		for (auto unitig : std::get<0>(it->second)) {
			printf("%llu ", unitig);
		}
		printf(" |");

		printf("%llu %llu |", std::get<1>(it->second), std::get<2>(it->second));

		for (auto neighbour : std::get<3>(it->second)) {
			printf("%llu ", neighbour);
		}

		printf(" |");
		printf("%llu ", std::get<4>(it->second));
		printf("\n");
	}*/

	//-------------


	//stvorimo thread_pool kako bi paralelizirali izvodenje koda
	std::shared_ptr <thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool();

	std::vector<std::future<void>> thread_futures;

	unsigned int size = sequences.size();

	unsigned int sequences_per_thread = size / number_of_threads;
	unsigned int residue_sequences = size % number_of_threads;

	unsigned int start = 0;
	unsigned int end = 0;

	//-------------

	//vektor za pohraniti koja sekvence se mora mapirati na koje regije grafa

        std::vector <std::unordered_map <unsigned long long int, std::vector <unsigned long long int>>> for_mapping;

        for (unsigned int i = 0; i < number_of_threads; i++) {
		std::unordered_map <unsigned long long int, std::vector <unsigned long long int>> single_thread_mapping;
                for_mapping.emplace_back(single_thread_mapping);
        }

	for (unsigned int i = 0; i < number_of_threads; i++) {
		if (residue_sequences > 0) {
			end = start + sequences_per_thread + 1;
			residue_sequences--;
		}
		else
			end = start + sequences_per_thread;

		thread_futures.emplace_back(thread_pool->submit(
			map_sequence_to_reference,
			//std::ref(reference_genome_minimizers_positions_in_vector),
			//std::ref(reference_genome_minimizers),
			std::ref(sequences),
			//std::ref(reference_genome),
			//std::ref(unitig_minimizers),
			hash_minimizer_region_number_to_add,
			regions,
			unitig_u_kojim_regijama,
			k,
			window_size,
			start,
			end,
			std::ref(for_mapping[i])
		));

		start = end;
	}

	for (auto &it : thread_futures)
		it.wait();


	//-------------

	//zapis svih regija/unitiga na koje ce se mapirati sekvence u zasebne fileove
	create_fasta_graph_files(regions, unitig_strings_for_creating_fasta_file, for_mapping);


	//TODO: pretvori fileove iz .fasta u .gfa
	//fa_to_gfa(path_for_bcalm_convert_to_gfa, sequences.size());


	//razdvoji sekvence koje treba mapirati u zasebne fileove

	//std::vector <std::string> sequences_file_paths = seperate_sequences_for_mapping_into_seperate_files (argv[optind]);


	//TODO: GraphAligner
	//run_GraphAligner(path_for_GraphAligner, sequences_file_paths, sequences.size());

	//printf("\n%.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n", time_spent_on_creating_regions, time_spent_on_matching, time_spent_on_creating_minimizers, time_spent_on_match, time_spent_on_sorting, time_spent_on_dfs);

	return 0;
}
