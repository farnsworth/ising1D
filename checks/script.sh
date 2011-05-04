#!/bin/sh

function dorun {
    ./check${1}.out file.in
    OUT=$?
    if [ $OUT = 0 ]; then
	echo "check "${1}" ...... OK"
    else
	echo "check "${1}" not passed, error " $OUT
	exit
    fi
}





rm -rf *.in

echo "CHECK 1: checks on static properties of U and V"
echo "CHECK 2: checks on properties of U(t) and V(t)"
echo "CHECK 3: checks of GGE occupations and gamma evolution (that should be constant)"
echo "CHECK 4: normalization of \sum_\alpha c_\alpha^2"
echo " "
echo "PBC"

cat > file.in << EOF
[system0]
size 100
h 0.9
J 1
gamma  1
epsilon  0.1
pbc 1
[system]
size 100
h   0.5
J   1
gamma  1
epsilon 0.1
pbc 1
EOF

dorun 1
dorun 2
dorun 3

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

dorun 4

echo " "
echo "OBC"
cat > file.in << EOF
[system0]
size 50
h 0.9
J 1
gamma  1
epsilon  0.0
pbc 0
[system]
size 50
h   0.5
J   1
gamma  1
epsilon 0.0
pbc 0
EOF

dorun 1
dorun 2
dorun 3

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

dorun 4