The following command solves a tree-ensemble optimization instance using three formulation (TEOC, TEOR, TEOM) described in Section 5.1 in the paper:
```
julia main_solver_single_instance.jl
```

Users can modify `main()` function to solve different problem instances. 
```
function main()
    dataset_dir = "dataset/concrete_bt/" # "concrete_bt", "concrete_rf", "redwine_bt", "redwine_rf"
    num_trees = 10
    is_LP = false
    time_limit_sec = 300
    for formulation in ["TEOM", "TEOR", "TEOC"]
        solve_single_instance(dataset_dir, num_trees, formulation, is_LP, time_limit_sec)
    end
end
```

* `dataset_dir`: Tree ensemble model dataset to use. 
  * There are four datasets in `dataset/`.
    * `dataset/concrete_bt`: Boosted trees using [concrete strength data](https://archive.ics.uci.edu/dataset/165/concrete+compressive+strength)
    * `dataset/concrete_rf`: Random forest using [concrete strength data](https://archive.ics.uci.edu/dataset/165/concrete+compressive+strength)
    * `dataset/redwine_bt`: Boosted trees using [redwine quality data](https://archive.ics.uci.edu/dataset/186/wine+quality)
    * `dataset/redwine_rf`: Random forest using [redwine quality data](https://archive.ics.uci.edu/dataset/186/wine+quality)
* `num_trees`: The number of decision trees
* `time_limit_sec`: Time limit for a MIP solve in seconds
  