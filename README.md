# BRACoD

Installation: 

    pip install BRACoD

If you want to use the R interface, install reticulate in R


Walkthrough

1. Simulate some data and normalize it

    sim_counts, sim_y, contributions = BRACoD.simulate_microbiome_counts(BRACoD.example_otu_data)
    sim_relab = BRACoD.scale_counts(sim_counts)


2. Run BRACoD

    trace = BRACoD.run_bracod(sim_relab, sim_y, n_sample = 1000, n_burn=1000, njobs=4)

3. Examine the diagnostics

    BRACoD.convergence_tests(trace, sim_relab)

4. Examine the results

    df_results = BRACoD.summarize_trace(trace, sim_counts.columns, 0.3)

5. Compare the results to the simulated truth

    bugs_identified = df_results["bugs"].values
    bugs_actual = np.where(contributions != 0)[0]

    precision, recall, f1 = BRACoD.score(bugs_identified, bugs_actual)
    print("Precision: {}, Recall: {}, F1: {}".format(precision, recall, f1))

6. Try with your real data. We have included some functions to help you threshold and process your data
    df_counts = BRACoD.threshold_count_data(df_counts)
    df_rel = BRACoD.scale_counts(df_counts)
    df_rel, Y = remove_null(df_rel, Y)
    trace = BRACoD.run_bracod(df_rel, Y, n_sample = 1000, n_burn=1000, njobs=4)
    df_results = BRACoD.summarize_trace(trace, sim_counts.columns, 0.3)
