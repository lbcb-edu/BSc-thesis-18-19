# Long_reads_mapper

Long_reads_mapper is C++ implementation for mapping long reads. 

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+ or clang 3.4+
2. cmake 3.2+

## Installation

To install Long_reads_mapper run to following commands:

```bash
git clone https://github.com/lbcb-edu/BSc-thesis-18-19.git
cd BSc-thesis-18-19
git checkout kpongracic
git submodule update --init --recursive
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

