#!/bin/sh

rm -rf *.in

cat > file.in << EOF
[system0]
size 10
h 0.9
J 1
gamma  1
epsilon  0.1
pbc 1
[system]
size 10
h   0.5
J   1
gamma  1
epsilon 0.1
pbc 1
EOF

./check1.out file.in 
#> /dev/null 
#2>&1

OUT=$?

if [ $OUT = 0 ]; then
    echo "check 1 (PBC) passed"
else
    echo "check 1 (PBC) not passed, error " $OUT
    exit
fi


cat > file.in << EOF
[system0]
size 10
h 0.9
J 1
gamma  1
epsilon  0.1
pbc 0
[system]
size 10
h   0.5
J   1
gamma  1
epsilon 0.1
pbc 0
EOF

./check1.out file.in 
#> /dev/null 
#2>&1
OUT=$?

if [ $OUT = 0 ]; then
    echo "check 1 (OBC) passed"
else
    echo "check 1 (OBC) not passed, error " $OUT
    exit
fi