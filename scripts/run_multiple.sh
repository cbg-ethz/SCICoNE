command=$1
n_trees=$2
n_tries=$3
input_tree=$4
scripts_dir=$5
region_sizes_file=$6

maximum_target=$7
initial_target=5
step_size=4
initial_beta=1
beta_rate=0.001

infer() {
    $command --postfix=$1 --tree_file=$2 --beta_rate=${beta_rate} --n_genotypes_target=$3 --beta=$4 > $1_logs.out
}

run_parallel() {
    for i in `seq 0 $(($2-1))`; do
        infer "$4-$1-$i" "$3" "$4" "$5" &
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

# Achieve convergence for each target until maximum_target is reached
target=${initial_target}
converged_tree=$input_tree
while true; do
    echo "Searching convergence on trees with maximum $target genotypes..."
    t=0
    while true; do
        block="$t"
        prefix="$target-$t"
        echo "Running block $t/$(($n_tries-1)) from ${input_tree} with maximum ${target} genotypes and beta=${initial_beta}..."

        # Run trees
        run_parallel "$t" "${n_trees}" "${input_tree}" "${target}" "${initial_beta}"

        # Check robustness
        rscore=$(evaluate_robustness "$prefix" | grep Robustness | cut -d " " -f 2)
        input_tree="$(evaluate_robustness "$prefix" | grep Maximum | cut -d " " -f 2)_tree_inferred.txt"
        best_prefix="$(echo $input_tree | cut -d_ -f1)"
        initial_beta="$(cat ${best_prefix}_beta_values.csv | rev | cut -d, -f1 | rev)"

        evaluate_robustness "$prefix"

        # # Check robustness
        # score=$(evaluate_robustness_2 "$block")
        # echo "Tree-distance-based robustness: ${score}" 
        # echo ${score} >> "tree_distances.txt"

        if (( $(echo "$rscore > 0.5" | bc) )); then
            converged_tree=$input_tree
            echo "Achieved convergence ($rscore)."
            break
        fi

        if [[ $t == $(($n_tries-1)) ]]; then
            echo "Could not converge in $n_tries tries for maximum $target genotypes."
            break 2
        fi

        t=$((t+1))

    done

    target=$((target+${step_size}))
    
    if [[ $target -gt ${maximum_target} ]]; then
        echo "Finished increasing the maximum number of genotypes."
        break
    fi

done

echo "\nLast converged tree: $converged_tree"
