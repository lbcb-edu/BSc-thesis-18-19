
#!/bin/bash

file=$1

cut -f 1,2,3,4,8,9 $file | sort -n -k 2 -r > ~/zavrsni/BSc-thesis-18-19/cleanedpaf.txt
