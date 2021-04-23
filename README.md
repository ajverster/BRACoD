# BRACoD

### Installation

Installation in python: 

    pip install BRACoD

There is also an R interface, which depends on the python version being installed. There is a helper function that will do it for you, but it might be easier to do it with pip.

    devtools::install_github("ajverster/BRACoD/BRACoD.R")

### Python Walkthrough

1. Simulate some data and normalize it

    ```python
    import BRACoD
    sim_counts, sim_y, contributions = BRACoD.simulate_microbiome_counts(BRACoD.example_otu_data)
    sim_relab = BRACoD.scale_counts(sim_counts)
    ```

2. Run BRACoD

    ```python
    trace = BRACoD.run_bracod(sim_relab, sim_y, n_sample = 1000, n_burn=1000, njobs=4)
    ```
    
3. Examine the diagnostics

    ```python
    BRACoD.convergence_tests(trace, sim_relab)
    ```

4. Examine the results

    ```python
    df_results = BRACoD.summarize_trace(trace, sim_counts.columns, 0.3)
    ```

5. Compare the results to the simulated truth

    ```python
    bugs_identified = df_results["bugs"].values
    bugs_actual = np.where(contributions != 0)[0]

    precision, recall, f1 = BRACoD.score(bugs_identified, bugs_actual)
    print("Precision: {}, Recall: {}, F1: {}".format(precision, recall, f1))
    ```

6. Try with your real data. We have included some functions to help you threshold and process your data
    
    ```python
    df_counts = BRACoD.threshold_count_data(df_counts)
    df_rel = BRACoD.scale_counts(df_counts)
    df_rel, Y = remove_null(df_rel, Y)
    trace = BRACoD.run_bracod(df_rel, Y, n_sample = 1000, n_burn=1000, njobs=4)
    df_results = BRACoD.summarize_trace(trace, sim_counts.columns, 0.3)
    ```
    
### R Walkthrough

1. Simulate some data and normalize it

    ```R
    library('BRACoD.R')
    data(obesity)
    r <- simulate_microbiome_counts(obesity)

    sim_counts <- r[[1]]
    sim_y <- r[[2]]
    contributions <- r[[3]]
    sim_relab <- scale_counts(sim_counts)
    ```

2. Run BRACoD

    ```R
    trace <- run_bracod(sim_relab, sim_y, n_sample = 1000, n_burn=1000, njobs=4)
    ```
    
3. Examine the diagnostics

    ```R
    convergence_tests(trace, sim_relab)
    ```

4. Examine the results

    ```R
    df_results <- summarize_trace(trace, colnames(sim_counts))
    ```

5. Compare the results to the simulated truth

    ```R
    bugs_identified <- df_results$bugs
    bugs_actual <- which(contributions != 0)

    r <- score(bugs_identified, bugs_actual)
    
    precision <- r[[1]]
    recall <- r[[2]]
    f1 <- r[[3]]

    print(sprintf("Precision: %.2f, Recall: %.2f, F1: %.2f",precision, recall, f1))
    ```

6. Try with your real data. We have included some functions to help you threshold and process your data
    
    ```R
    df_counts <- threshold_count_data(df_counts)
    df_rel <- scale_counts(df_counts)
    r <- remove_null(df_rel, Y)
    df_rel <- r[[1]]
    Y <- r[[2]]
    
    trace <- run_bracod(df_rel, Y, n_sample = 1000, n_burn=1000, njobs=4)
    df_results <- summarize_trace(trace, sim_counts.columns, 0.3)
    ```

