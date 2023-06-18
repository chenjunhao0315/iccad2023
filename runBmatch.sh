if [[ -z "$1" ]]; then
    echo "Usage: sh ./runBmatch.sh <test case>"
    exit 1
fi

make -j8

./abc -c "bmatch -v -r benchmark/case$1/cir1.v benchmark/case$1/cir2.v result.txt"
