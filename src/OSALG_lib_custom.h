#include <string>
#include <vector>

namespace OSALG {
	int long_gaps_alignment(std::string const &seq1, std::string const &seq2, int L, std::vector<int> const &u, std::vector<int> const &v, std::vector<std::string> &cigars,
				int match_score, int mismatch_score);
}
