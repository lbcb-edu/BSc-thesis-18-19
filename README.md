# A Benchmark of Tools for Metagenomic Species Identification

Python 3.6.5. scripts for simulating metagenomic data sets and evaluating outputs of classification tools [Kraken](https://github.com/DerrickWood/kraken) and [CLARK](http://clark.cs.ucr.edu/Tool/).

## Installation

To build scripts run the following commands:

```bash
git clone --recursive --single-branch -b jlipovac https://github.com/lbcb-edu/BSc-thesis-18-19
cd src
```

## Usage

For simulating test data set you can use data from folder [test_data_for_simulator](https://drive.google.com/file/d/1P8K1lcpy46FMYSz3-A8q1YnC1zQ9mcIx/view) and run simulator script like this:

```bash
python simulator.py -f ../test_data_for_simulator/klebsiella_pneumoniae_INF125_hiseq_2000.fasta -s ../test_data_for_simulator/pantonea_agglomerans_miseq.fasta -t ../test_data_for_simulator/pseudomonas_koreensis_miseq_2_300.fasta -n ../test_data_for_simulator/human_chr.fasta 
```
For advanced simulator options check help:

```bash
python simulator.py -h
```
or :

```bash
python simulator.py --help
```

In order to see how evaluator works, use test outputs from [test_data_for_evaluator](https://drive.google.com/file/d/1t90iMjv2Rt8Sv5GD5g4E2lE6GmGoVCMT/view):

```bash
python evaluator.py -k ../test_data_for_evaluator/testdataset.kra.txt -c ../test_data_for_evaluator/clarktestdataset.csv -p ../test_data_for_evaluator/testdatasetPercentage.txt -i ../test_data_for_evaluator/testdatasetIDs.txt
```
For advanced evaluator options check help:

```bash
python evaluator.py -h
```
or :

```bash
python evaluator.py --help
```


# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.
