bits=(14) # half the size of N
dim=(100 200 300)
blks=(50)

for size in "${bits[@]}"
do
    for n in "${dim[@]}"
    do
        for bz in "${blks[@]}"
        do
            ./test $n 10000 $bz $size 0.714
            python plot.py
        done
    done
done