# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Installation

Clone repository, update submodules and build files with:
```bash
git submodule update --init --recursive
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Usage

Program is consisted of two parts. To make .msa files position yourself into python-scripts folder.
Run: 

```python
python make_msa2.py <folder_with_fastq_files>
```
You will find new spoa_folder with generated .msa files.
Position yourself within /build folder and run:

```bash
./bin/zavrsni <path_to_.msa_file> <path_to_appropriate_fastqFile>
```

## Usage

After building, executable file is located in /bin folder as detection. For detecting chimeric and repeating reads and create .fasta files with separated and annotated reads run:
```bash
./bin/detection file_with_mappings.paf reads.fasta
```

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.
