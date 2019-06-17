#ifndef BLUE_ALIGNMENT_H_INCLUDED
#define BLUE_ALIGNMENT_H_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <stdio.h>

namespace blue
{
    enum AlignmentType {global, local, semi_global};    

    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match, int mismatch, int gap,
                           std::string &cigar, unsigned int &target_begin);

    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match, int mismatch, int gap);

    AlignmentType getType(std::string type);
}

#endif // BLUE_ALIGNMENT_H_INCLUDED
