# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Installation

To install OSALG run commands:

```bash
git clone --recursive --single-branch -b rjpenic https://github.com/lbcb-edu/BSc-thesis-18-19
cd BSc-thesis-18-19
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## Usage

After installation, executable file will be located in /build/bin folder. To use OSALG run:
```bash
cd bin
./OptSeqAlignmentLongGaps reads.fasta references.fasta
```

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.
