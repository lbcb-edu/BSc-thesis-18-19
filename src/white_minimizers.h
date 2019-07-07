#pragma once

namespace white {

	std::vector<std::tuple<unsigned long long int, unsigned long long int, bool>> minimizers(
										const char* sequence,
										unsigned long long int sequence_length,
                                                                     		unsigned long long int k,
                                                                     		unsigned long long int window_length);
}
