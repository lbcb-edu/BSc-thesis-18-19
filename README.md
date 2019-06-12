# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Installation

Tool can be installed and build by running:

```bash
git clone --recursive --single-branch -b sbakic https://github.com/lbcb-edu/BSc-thesis-18-19
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Usage

After building, executable file is located in /bin folder as detection. For detecting chimeric and repeating reads and create .fasta files with separated and annotated reads run:
```bash
./bin/detection file_with_mappings.paf reads.fasta
```

## Data

The data set for testing is freely available from Loman Labs [here](https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta), while the reference genome is freely available from NCBI [here](https://bit.ly/2PCYHWr).

You use data set with reads and reference data set for mapping reads to reference with minimap2. The result of mapping should then be forwarded to detection.

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.
