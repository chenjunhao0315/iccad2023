if [[ -z "$1" ]]; then
    echo "Usage: sh ./runBmatch.sh <test case>"
    exit 1
fi

make -j8

./abc -c "bmatchgroup -v -r benchmark/case$1/input result.txt"
