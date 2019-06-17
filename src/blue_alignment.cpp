#include "blue_alignment.hpp"
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <stdio.h>

namespace blue
{
    typedef struct {
        int value;
        char trace;
    } Cell;

    void create_cigar_string(std::vector<std::vector<Cell>> &matrix, int i, int j, std::string &cigar, unsigned int &target_begin, const char* query, unsigned int query_length, const char* target);

    void initialize_matrix(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           int match, int mismatch, int gap, std::vector<std::vector<Cell> > &matrix, AlignmentType type);

    int semi_global_alignment(const char* query, unsigned int query_length,
                                    const char* target, unsigned int target_length,
                                    int match, int mismatch, int gap, std::string &cigar, unsigned int &target_begin);

    int local_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap, std::string &cigar, unsigned int &target_begin);

    int global_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap, std::string &cigar, unsigned int &target_begin);

    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match, int mismatch, int gap) {

                           std::string overloaded;
                           unsigned int broj;
                    return pairwise_alignment(query, query_length, target, target_length, type, match, mismatch, gap, overloaded, broj);
    }

    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match, int mismatch, int gap,
                           std::string &cigar, unsigned int &target_begin) {

        switch (type) {
            case global:
                return global_alignment(query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);
            case local:
               return local_alignment(query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);

            case semi_global:
                return semi_global_alignment(query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);

            default:
                std::cout << "Wrong Alignment type!" << std::endl;
                break;
        }
        return -1;
    }
    //************************************************************************************************************************
    int local_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap,
                                 std::string &cigar, unsigned int &target_begin) {

        std::vector<std::vector<Cell>> matrix;
        matrix.resize(query_length+1, std::vector<Cell>(target_length+1));

        std::pair<int, int> max_cell;
        max_cell.first = 0;
        max_cell.second = 0;
        int max_val = 0;

        AlignmentType type = local;

        initialize_matrix(query, query_length, target, target_length, match, mismatch, gap, matrix, type);

        for(int i = 1; i < query_length+1; i++) {
            for(int j = 1; j<target_length+1; j++ ) {

                if (matrix[i][j].value > max_val) {
                    max_val = matrix[i][j].value;
                    max_cell.first = i;
                    max_cell.second = j;
                }
            }
        }

        create_cigar_string(matrix, max_cell.first, max_cell.second, cigar, target_begin, query, query_length, target);
        std::vector<std::vector<Cell>>().swap(matrix);
        return max_val;
    }

    int global_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap,
                                 std::string &cigar, unsigned int &target_begin) {
        std::vector<std::vector<Cell> > matrix;
        matrix.resize(query_length+1, std::vector<Cell>(target_length+1));

        AlignmentType type = global;
        initialize_matrix(query, query_length, target, target_length, match, mismatch, gap, matrix, type);

        create_cigar_string(matrix, query_length, target_length, cigar, target_begin, query, query_length, target);
        return matrix[query_length][target_length].value;
    }

    int semi_global_alignment(const char* query, unsigned int query_length,
                                    const char* target, unsigned int target_length,
                                    int match, int mismatch, int gap,
                                    std::string &cigar, unsigned int &target_begin) {

        std::vector<std::vector<Cell> > matrix;
        matrix.resize(query_length+1, std::vector<Cell>(target_length+1));

        AlignmentType type = semi_global;

        initialize_matrix(query, query_length, target, target_length, match, mismatch, gap, matrix, type);

        std::pair<int, int> max_cell;
        max_cell.first = 0;
        max_cell.second = 0;
        int max_val = matrix[0][0].value;

        for (int i = 1; i < query_length + 1; i++) {
            if (matrix[i][target_length].value > max_val) {
                max_val = matrix[i][target_length].value;
                max_cell.first = i;
                max_cell.second = target_length;
            }
        }

        for (int i = 1; i < target_length + 1; i++) {
            if (matrix[query_length][i].value > max_val) {
                max_val = matrix[query_length][i].value;
                max_cell.first = query_length;
                max_cell.second = i;
            }
        }
        create_cigar_string(matrix, max_cell.first, max_cell.second, cigar, target_begin, query, query_length, target);
        return max_val;
    }
    //************************************************************************************************************************

    void initialize_matrix(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           int match, int mismatch, int gap, std::vector<std::vector<Cell> > &matrix, AlignmentType type) {

        matrix[0][0].value = 0;
        matrix[0][0].trace = 'n';

        if (type == global) {
            for(int i = 1; i < query_length+1; i++) {
                matrix[i][0].value = gap*i;
                matrix[i][0].trace = 'u';
            }

            for(int i = 1; i<target_length + 1; i++) {
                matrix[0][i].value = gap*i;
                matrix[0][i].trace = 'l';
            }
        } else {
            for (int i=1; i < query_length + 1; i++){
                matrix[i][0].value = 0;
                matrix[i][0].trace = 'n';
            }

            for (int i=1; i < target_length + 1; i++){
                matrix[0][i].value = 0;
                matrix[0][i].trace = 'n';
            }
        }

        for(int i = 1; i < query_length+1; i++) {
            for(int j = 1; j<target_length+1; j++ ) {
                int match_mismatch = query[i-1] == target[j-1] ?
                                        matrix[i-1][j-1].value + match : matrix[i-1][j-1].value + mismatch;
                int insertion = matrix[i][j-1].value + gap;
                int deletion = matrix[i-1][j].value + gap;

                if (match_mismatch < 0 && insertion < 0 && deletion < 0 && type == local) {
                    matrix[i][j].value = 0;
                    matrix[i][j].trace = 'n';
                } else {
                    if(insertion > deletion) {
                        if(match_mismatch > insertion) {
                            matrix[i][j].value = match_mismatch;
                            matrix[i][j].trace = 'd';

                    } else if (match_mismatch == insertion) {
                        if (matrix[i-1][j-1].value > matrix[i][j-1].value) {
                            matrix[i][j].trace = 'd';
                        } else {
                            matrix[i][j].trace = 'l';
                        }
                        matrix[i][j].value = match_mismatch;

                    } else {
                        matrix[i][j].value = insertion;
                        matrix[i][j].trace = 'l';
                    }

    //***********************************************************************************
                    } else {  //deletion > insertion
                        if(match_mismatch > deletion) {
                            matrix[i][j].value = match_mismatch;
                            matrix[i][j].trace = 'd';

                        } else if (match_mismatch == deletion) {
                            if (matrix[i-1][j-1].value > matrix[i-1][j].value) {
                                matrix[i][j].trace = 'd';
                            } else {
                                matrix[i][j].trace = 'u';
                            }
                            matrix[i][j].value = match_mismatch;

                        } else {
                            matrix[i][j].value = deletion;
                            matrix[i][j].trace = 'u';
                        }
                    }
                }
            }
        }
    }

    void create_cigar_string(std::vector<std::vector<Cell> > &matrix, int start_row, int start_column,
                             std::string &adress, unsigned int &target_begin,
                               const char* query, unsigned int query_length, const char* target) {

        int mism_counter, mis_counter, del_counter, ins_counter;
        mism_counter = mis_counter = del_counter = ins_counter = 0;

        std::string cigar_reverse = "";
        Cell current = matrix[start_row][start_column];
        int currentRow = start_row;
        int currentColumn = start_column;

        if (start_row < query_length) {
            cigar_reverse+='S';
            cigar_reverse+= std::to_string(query_length - start_row);
        }

        std::pair<int, int> pos;
        pos.first = start_row;
        pos.second = start_column;

        while(current.trace != 'n') {
            //provjera za match/mismatch
            if(current.trace == 'd') {
                if(del_counter != 0) {
                    cigar_reverse+=std::to_string(del_counter);
                    del_counter = 0;
                } else if (ins_counter != 0) {
                    cigar_reverse +=std::to_string(ins_counter);
                    ins_counter = 0;
                }

                if(query[pos.first-1] == target[pos.second-1]){
                    if(mism_counter != 0) {
                        cigar_reverse+=std::to_string(mism_counter);
                        mism_counter = 0;
                    }
                    if(mis_counter == 0) {
                        cigar_reverse+='=';
                    }
                    ++mis_counter;
                } else {
                    if(mis_counter != 0) {
                        cigar_reverse+=std::to_string(mis_counter);
                        mis_counter = 0;
                    }
                    if(mism_counter == 0) {
                        cigar_reverse+='X';
                    }
                    ++mism_counter;
                }

                // provjera za insertion
            } else if(current.trace == 'l') {
                if(del_counter != 0) {
                    cigar_reverse+=std::to_string(del_counter);
                    del_counter = 0;
                } else if (mis_counter != 0) {
                    cigar_reverse +=std::to_string(mis_counter);
                    mis_counter = 0;
                } else if (mism_counter != 0) {
                    cigar_reverse +=std::to_string(mism_counter);
                    mism_counter = 0;
                }
                if(ins_counter == 0) {
                    cigar_reverse+='I';
                }
                ++ins_counter;

                //deletion
            } else {
                if(ins_counter != 0) {
                    cigar_reverse+=std::to_string(ins_counter);
                    ins_counter = 0;
                } else if (mis_counter != 0) {
                    cigar_reverse +=std::to_string(mis_counter);
                    mis_counter = 0;
                } else if (mism_counter != 0) {
                    cigar_reverse +=std::to_string(mism_counter);
                    mism_counter = 0;
                }
                if(del_counter == 0) {
                    cigar_reverse+='D';
                }
                ++del_counter;
            }

            if (current.trace == 'd'){
                pos.first = currentRow - 1;
                pos.second = currentColumn - 1;
            } else if(current.trace == 'u'){
                pos.first = currentRow - 1;
                pos.second = currentColumn;
            } else if(current.trace == 'l'){
                pos.first = currentRow;
                pos.second = currentColumn - 1;
            }

            if (current.trace == 'd'){
                current = matrix[--currentRow][--currentColumn];
            } else if(current.trace == 'u'){
                current = matrix[--currentRow][currentColumn];
            } else if(current.trace == 'l'){
                current = matrix[currentRow][--currentColumn];
            }
        }

       target_begin = (unsigned int) pos.second;

        if(mis_counter != 0) {
            cigar_reverse+=std::to_string(mis_counter);
        } else if(mism_counter != 0){
            cigar_reverse+=std::to_string(mism_counter);
        } else if(ins_counter != 0){
            cigar_reverse+=std::to_string(ins_counter);
        } else if(del_counter != 0) {
            cigar_reverse+=std::to_string(del_counter);
        }

        if (pos.first != 0) {
            cigar_reverse+='S';
            cigar_reverse+=std::to_string(pos.first);
        }

        std::string number = "";
        std::string cigar = "";
        for(int i=cigar_reverse.length()-1; i>=0; i--) {
            char c = cigar_reverse.at(i);
            if(isdigit(c)) {
                number.insert(0, 1, c);
            } else {
                cigar += number;
                cigar += c;
                number = "";
            }
        }

        adress = cigar;
        return;
    }

    AlignmentType getType(std::string type) {
        if(type.compare("global") == 0) return global;
        if(type.compare("local") == 0) return local;
        if(type.compare("semi_global")==0) return semi_global;
    }
}
