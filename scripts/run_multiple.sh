command=$1
n_trees=$2
n_tries=$3
input_tree=$4
scripts_dir=$5
region_sizes_file=$6

infer() {
    $command --postfix=$1 --tree_file=$2 > $1_logs.out
}

run_parallel() {
    for i in `seq 0 $(($2-1))`; do
        infer "$1-$i" "$3" "$i" &
    done
    wait
}

evaluate_robustness() {
    thres=5

    trees=($( ls -1 $1*_tree_inferred.txt | cut -d "_" -f 1))
    scores=($( cat $1*_tree_inferred.txt | grep "Tree posterior" | cut -d " " -f 3 ))
    num=${#trees[@]}

    max=$(printf '%s\n' "${scores[@]}" | awk '$1 > m || NR == 1 { m = $1 } END { print m }')

    max_tree=${trees[0]}
    n_maxs=0
    for i in ${!scores[@]};
    do
       score_diff=$(awk -v a=${scores[$i]} -v b=$max "BEGIN {print b-a}"  OFMT="%.20f")
       is_robust=$(awk -v a=${score_diff} -v b=$thres "BEGIN {if ( a < b ) print 1; else print 0}")
       if [ ${is_robust} == "1" ]
       then
           n_maxs=$((n_maxs+1))
           max_tree=${trees[$i]}
       fi

       echo "$i\t${trees[$i]}\t${scores[$i]}\t${score_diff}\t${is_robust}"
    done
    robustness_score=$(awk -v a=$num -v b=${n_maxs} "BEGIN {print b/a}")

    echo ""
    echo "Maximum: ${max_tree} ${max}"
    echo "Robustness: ${robustness_score} ($n_maxs/$num)"
}


evaluate_robustness_2() {
    scores=($( cat $1*_tree_inferred.txt | grep "Tree posterior" | cut -d " " -f 3 | tr '\n' ',' | sed 's/.$//'))

    # Compute cell-cell distance for each tree
    for i in `seq 0 $((${n_trees}-1))`; do
        python ${scripts_dir}/compute_pdist.py ${region_sizes_file} "$1"-"$i"_tree_inferred.txt "$1"-"$i"_cell_node_ids.tsv
    done

    files=($(for i in `seq 0 $((${n_trees}-1))`; do echo $1-"$i"_pdist.txt; done | tr '\n' ',' | sed 's/.$//'))
    # Compute tree-tree distances
    python ${scripts_dir}/compute_tree_distances.py ${files} --scores ${scores}
}

for i in `seq 0 $((${n_tries}-1))`; do
    block="$i"
    echo "Running block $i from ${input_tree}..."

    # Run trees
    run_parallel "$i" "${n_trees}" "${input_tree}"

    # Check robustness
    rscore=$(evaluate_robustness "$block" | grep Robustness | cut -d " " -f 2)
    input_tree="$(evaluate_robustness "$block" | grep Maximum | cut -d " " -f 2)_tree_inferred.txt"

    evaluate_robustness "$block"

    # Check robustness
    score=$(evaluate_robustness_2 "$block")
    echo "Tree-distance-based robustness: ${score}" 
    echo ${score} >> "tree_distances.txt"

    if [[ $rscore == "1" ]]; then
        break
    fi
done


