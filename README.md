# LORMAP

LORMAP is C++ implementation for mapping long reads. 

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+ or clang 3.4+
2. cmake 3.2+

## Installation

To install LORMAP run the following commands:

```bash
git clone --recursive --single-branch -b kpongracic https://github.com/lbcb-edu/BSc-thesis-18-19.git lormap
cd lormap
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running make, an executable named `lormap` will appear in the `build` directory.

## Usage

Usage of LORMAP is as following:

```bash
./lormap [OPTIONS] <reference> <sequences>
	
    <reference>
        FASTA file containing reference genome (can be compressed with gzip)
    <sequences>
        FASTA/FASTQ file containing a set of fragments (can be compressed with gzip)

    OPTIONS:
        -K  or  --kvalue <int>
            default: 1000
            number of bases to take from start and end
        -w  or  --window_length <int>
            default: 5
            length of window
        -k  or  --kmers <int>
            default: 15
            constraints: largest supported is 16
            number of letters in substrings
        -r  or  --bandwidth <int>
            default: 500
            bandwidth used in chaining and DP-based alignment
        -n  or  --min_hits <int>
            default: 5
            minimum hits when finding region
        -A  or  --match <int>
            default: 2
            match number
        -B  or  --mismatch <int>
            default: 4
            mismatch number
        -O  or  --gap_open <int>
            default: 4
            gap opening penalty
        -E  or  --gap_extension <int>
            default: 2	
            gap extension penalty
        -z  or  --z_drop <int>	
            default: 400	
            Z-drop score
        -t  or  --threads <int>
            default: 3
            number of threads
        -c  or  --cigar 
            output CIGAR in PAF
        -a  or  --sam
            output in the SAM format (PAF by default)
        -v or --version
            prints the version number
        -h or --help
            prints the usage			
```

## Example

The first version of the mapper was tested on an Oxford Nanopore Technologies data set obtained by sequencing the Escherichia coli K-12 substr. MG1655 genome. The data set is freely available from Loman Labs [here](https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta), while the reference genome is freely available from NCBI [here](https://bit.ly/2PCYHWr).

# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.
