if [[ -z "$1" ]]; then
    echo "Usage: sh ./runBmatchProgram.sh <test case>"
    exit 1
fi

make -j32 libabc.a
gcc -Wall -g -c bmatch.c -o bmatch.o
g++ -g -o bmatch bmatch.o libabc.a -lm -ldl -lreadline

./bmatch $1
