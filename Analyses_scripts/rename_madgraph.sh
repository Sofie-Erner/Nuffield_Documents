#!/bin/bash
# this script will rename madgraph run_X directories and files depending on variables used
#

input_file=params.dat
g=0
m=0
IFS=' = '

for path in run_*; do
    echo $path
    cd $path

    iter=0
    while IFS= read -r line; do
        #echo $line
        read -ra ValArr <<< "$line"

        if [ $iter == 0 ]; then
            g=${ValArr[1]}
            echo $g
            if [ "${g%.*}" -gt "100" ]; then
                g_temp=$(printf '%.1e' $g)
                g=$g_temp
            fi
        elif [ $iter == 1 ]; then
            m=${ValArr[1]}
            echo $m
        fi 

        iter=$((iter+1))
    done < $input_file

    bind="_"
    name="$m$bind$g"
    echo $name

    #if [[ -f "*.yoda" ]]; then
        #mv *.yoda $name.yoda
        gunzip *.lhe.gz
        mv *.lhe $name.lhe
        #mv *.hepmc.gz $name.hepmc.gz
        #mv *.txt $name.txt
        #mv un*.dat $name.dat
        mv E_theta.dat ../LHE_files/e_theta_rivet/$name.dat 
        #echo "done"
    #fi

    cd ..
done
