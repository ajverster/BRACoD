
import numpy as np
import pandas as pd
import argparse
import BRACoD

def simulate_microbiome_counts(df_counts, n_contributors = 20, coeff_contributor = 0.0, min_ab_contributor = -9, sd_Y = 1.0, n_reads = 100000, var_contributor = 5.0, use_uniform = True, n_samples_use = None, corr_value = None, return_absolute = False, seed = None):
    if seed is None:
        np.random.seed(seed)

    df_relab = df_counts.apply(lambda x: x / np.sum(x), 1)
    if n_samples_use is None:
        n_samples_use = df_counts.shape[0]
    # lowest value
    pseudo = np.min(np.min(df_relab[df_relab > 0])) / 10.0

    # Randomize the order of bacteria
    order = np.array(range(df_counts.shape[1]))
    np.random.shuffle(order)
    df_relab = df_relab.iloc[:,order]

    # Threshold to only those high abundant microbes
    mu = df_relab.apply(lambda x: np.mean(np.log(x + pseudo)), 0)
    df_relab = df_relab.loc[:,mu > min_ab_contributor]
    df_relab.reset_index(inplace = True, drop=True)
    df_relab.columns = df_counts.columns[mu > min_ab_contributor]
    df_relab.index = df_counts.index

    # Average and sd of each bug
    mu = df_relab.apply(lambda x: np.mean(np.log(x + pseudo)), 0)
    var = df_relab.apply(lambda x: np.var(np.log(x + pseudo)), 0)

    Sigma = np.diag(var)
    # Add in some inter-microbe correlations if you want
    if corr_value is not None:
        corr_value = np.float(corr_value)
        for i in range(n_contributors):
            # First 10 bugs are contributors, next 10 bugs are correlated
            j = i + n_contributors
            Sigma[i,j] = corr_value * np.sqrt(var[i] * var[j])
            Sigma[j,i] = corr_value * np.sqrt(var[j] * var[i])

    # Simulate log absolute abundances
    sim_absolute_log = np.random.multivariate_normal(mean=mu, cov=Sigma, size=int(n_samples_use))
    n_decoys = sim_absolute_log.shape[1] - n_contributors

    X = sim_absolute_log

    # Simulate environmental values
    if var_contributor is None:
        var_contributor = coeff_contributor * 0.25
    if not use_uniform:
        beta_contributor = list( np.random.normal(coeff_contributor, scale = var_contributor, size = int(n_contributors) ))
    else:
        beta_contributor = list( np.random.uniform(coeff_contributor - var_contributor, coeff_contributor + var_contributor, size = int(n_contributors) ))

    beta = np.array(beta_contributor + [0 for i in range(int(n_decoys))])

    Y = np.random.normal(X.dot(beta), scale=sd_Y)
    Y = Y - np.mean(Y)

    if return_absolute:
        return X, Y, beta

    X_nonlog = np.exp(X)
    X_rel = np.apply_along_axis(lambda x: x / np.sum(x), 1, X_nonlog)
    X_counts = np.apply_along_axis(lambda x: np.random.multinomial(int(n_reads), x), 1, X_rel)
    df_counts = pd.DataFrame(X_counts)
    df_counts.columns = df_relab.columns
    df_counts.index = df_relab.index[:df_counts.shape[0],]
    return df_counts, Y, beta


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile_otus", help="csv of 16S sequencing counts that will be used as a model for data simulation")
    args = parser.parse_args()

    df_counts = pd.read_csv(args.infile_otus)
    print(df_counts)
    sim_counts, sim_y, contributions =  simulate_microbiome_counts(df_counts)
    print(sim_counts)
    sim_relab = BRACoD.scale_counts(sim_counts)

    trace = BRACoD.run_bracod(sim_relab, sim_y)

