# BRACoD: Bayesian Regression Analysis of Compositional Data

### Installation

Installation in python: 

    pip install BRACoD

There is also an R interface, which depends on the python version being installed. There is a helper function that will do it for you, but it might be easier to do it with pip. You can install via CRAN.

    install.packages("BRACoD.R")

Or you can install via github.

    devtools::install_github("ajverster/BRACoD/BRACoD.R")

### Python Walkthrough

1. Simulate some data and normalize it

    ```python
    import BRACoD
    import numpy as np
    sim_counts, sim_y, contributions = BRACoD.simulate_microbiome_counts(BRACoD.df_counts_obesity)
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
    taxon_identified = df_results["taxon_num"].values
    taxon_actual = np.where(contributions != 0)[0]

    precision, recall, f1 = BRACoD.score(taxon_identified, taxon_actual)
    print("Precision: {}, Recall: {}, F1: {}".format(precision, recall, f1))
    ```

6. Try with your real data. We have included some functions to help you threshold and process your data
    
    ```python
    df_counts = BRACoD.threshold_count_data(BRACoD.df_counts_obesity)
    df_rel = BRACoD.scale_counts(df_counts)
    df_rel, Y = BRACoD.remove_null(df_rel, BRACoD.df_scfa_obesity["butyric"].values)
    trace = BRACoD.run_bracod(df_rel, Y, n_sample = 1000, n_burn=1000, njobs=4)
    df_results = BRACoD.summarize_trace(trace, df_rel.columns, 0.3)
    ```

The taxonomy information for these OTUs is available at ```BRACoD.df_taxonomy```
    
### R Walkthrough

1. Simulate some data and normalize it

    ```R
    library('BRACoD.R')
    data(obesity)
    r <- simulate_microbiome_counts(df_counts_obesity)

    sim_counts <- r$sim_counts
    sim_y <- r$sim_y
    contributions <- r$contributions
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
    taxon_identified <- df_results$taxon_num
    taxon_actual <- which(contributions != 0)

    r <- score(taxon_identified, taxon_actual)
    
    print(sprintf("Precision: %.2f, Recall: %.2f, F1: %.2f", r$precision, r$recall, r$f1))
    ```

6. Try with your real data. We have included some functions to help you threshold and process your data
    
    ```R
    df_counts_obesity_sub <- threshold_count_data(df_counts_obesity)
    df_rel <- scale_counts(df_counts_obesity_sub)
    r <- remove_null(df_rel, df_scfa$butyric)
    df_rel <- r$df_rel
    Y <- r$Y
    
    trace <- run_bracod(df_rel, Y, n_sample = 1000, n_burn=1000, njobs=4)
    df_results <- summarize_trace(trace, colnames(df_counts_obesity_sub), 0.3)
    ```

