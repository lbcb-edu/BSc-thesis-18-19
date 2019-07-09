# Introduction

MapToGraph is a C++ implementation of a sequence to de Brujin graph mapper. For each sequence it finds a subgraph in the original de Brujin graph to which the sequence maps. The mapper is still in early stages of development. Final version should map sequences to subgraphs using [GraphAligner][gl]. The implementation was made as a part of my Final BSc thesis.

# Dependencies

1. gcc 4.8+
2. cmake 3.2+

# Installation

In order to build MapToGraph from source run the following commands:

```bash
git clone --recursive --single-branch -b dbatic https://github.com/lbcb-edu/BSc-thesis-18-19.git GraphMap
cd MapToGraph
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running the commands, the executable file will be created in `build/bin`.

# Input

GraphAligner outputs mappings in the [FASTA format][fasta].

# Usage

In order to map sequences use the following command:
```bash 
graphmap graph.fa reads.fa / reads.fq
```

# Input

GraphMap takes two imputs:
- a de Brujin graph in FASTA format
- reads in FASTA or FASTQ format

# Output

GraphAligner outputs a subgraph for each sequence as a seperate file in FASTA format.

# Input data test examples

Escherichia coli str. K-12 substr. MG1655
  - [reference genome][ref]

# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.

[gl]: https://github.com/maickrau/GraphAligner
[fasta]: https://en.wikipedia.org/wiki/FASTA_format
[ref]: https://www.ncbi.nlm.nih.gov/genome/167

