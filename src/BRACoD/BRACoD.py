
import pymc3 as pm
from arviz.stats import diagnostics
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale

import typing
import logging

import theano
#theano.config.blas.ldflags = '-lf77blas -latlas -lgfortran'

def scale_counts(df):
    """
    Takes a DataFrame of OTU counts and normalizes it for use with run_bracod()
    :param df: DataFrame of OTU counts
    :return: A DataFrame of relative abundance data
    """
    df = pd.DataFrame(df)
    assert 'int' in str(df.iloc[:,0].dtype), "This is not counts data"
    # Convert to relative abundance
    df = df.apply(lambda x: x / np.sum(x),1)

    # Add a pseudo count
    df[df == 0] = np.min(df[df > 0]).min() / 10.0
    # Renorm
    df = df.apply(lambda x: x / np.sum(x),1)

    return df


def check_chains_equal(trace):
    """
    Checks whether there are any OTUs whose inclusion metric is radically different between chains
    :param trace:
    :return:
    """
    n_chains = trace.nchains
    inclusion_perchain = np.zeros((trace.get_values('p').shape[1], n_chains))

    for i in range(n_chains):
        inclusion_perchain[:,i] = trace.get_values('p', chains=[i]).sum(0) / trace.get_values('p', chains=[i]).shape[0]
        
    # Your chains are radically different
    radically_different = np.where(np.apply_along_axis(lambda x: np.max(x) - np.min(x),1,inclusion_perchain) > 0.5)[0]
    if len(radically_different) > 0:
        return radically_different
    return None


def get_positives(trace, inclusion_cutoff=0.30):
    """
    Return locations of contributing OTUs
    :param trace:
    :param inclusion_cutoff:
    :return:
    """
    inclusion_full = trace.get_values('p').sum(0) / trace.get_values('p').shape[0]
    found_full = np.where(inclusion_full >= inclusion_cutoff)[0]
    return inclusion_full, found_full


def convergence_tests(trace, df_otus = None, inclusion_cutoff=0.30):
    """
    This function runs a series of convergence tests for the important variables in Bannoc
    Meant to replace the less focused warnings that come from pm.sample
    :param trace:
    :param inclusion_cutoff:
    :return:
    """
    # Calculate the inclusion probabilities
    if trace.nchains >= 2:
        inclusion_full, pos_values = get_positives(trace)
        all_chains = list(range(trace.nchains))
        halfway = int(trace.nchains/2)
        inclusion_1 = trace.get_values('p', chains=all_chains[:halfway]).sum(0) / trace.get_values('p', chains=all_chains[:halfway]).shape[0]
        inclusion_2 = trace.get_values('p', chains=all_chains[halfway:]).sum(0) / trace.get_values('p', chains=all_chains[halfway:]).shape[0]
        found1 = set(np.where(inclusion_1 >= inclusion_cutoff)[0])
        found2 = set(np.where(inclusion_2 >= inclusion_cutoff)[0])
        found_uncertain = np.array(list((found1 - found2).union(found2 - found1)))
        # Only accept the uncertain if the difference in inclusion probabilities is substantial
        if len(found_uncertain) > 0:
            difference = np.absolute(inclusion_1[found_uncertain] - inclusion_2[found_uncertain])
        
            if len(found_uncertain[difference > 0.1]) > 0:
                 logging.warning("Warning! The following bacteria were found in only one of the chains: {}".format(
                    " ".join([str(x) for x in found_uncertain])))

        diff = check_chains_equal(trace)
        if diff is not None:
            logging.warning("Warning! the following bacteria are radically different between chains {}".format(" ".join([str(x) for x in diff])))


        # Check the effective number of samples for the positive bugs
        ess_data = diagnostics.ess(trace)
        rhat_values = diagnostics.rhat(trace)
        effn_p = ess_data["p"] 
        effn_betas_one = ess_data["betas_one"] 
        gr_p = rhat_values["p"]
        gr_betas_one = rhat_values["betas_one"]
        if (any(effn_p[pos_values] < 200) | any(effn_betas_one[pos_values] < 200)):
            logging.warning("Warning! Some parameters have an effective sample size less than 200.")
            logging.warning("Either Rerun with more steps or be wary of including them in your interpretation.")
            if any(effn_p[pos_values] <= 200):
                problem_bugs = pos_values[np.where(effn_p[pos_values] < 200)[0]]
                if df_otus is not None:
                    problem_bugs = df_otus.columns[problem_bugs]
                logging.warning("The effective n for the p variable is less than 200 for the following bacteria {}".format(problem_bugs))
            if any(effn_betas_one[pos_values] <= 200):
                problem_bugs = pos_values[np.where(effn_betas_one[pos_values] < 200)[0]]
                if df_otus is not None:
                    problem_bugs = df_otus.columns[problem_bugs]
                logging.warning("The effective n for the beta variable is less than 200 for the following bacteria {}".format(problem_bugs))

        if (any(gr_p[pos_values] >= 1.2) | any(gr_betas_one[pos_values] >= 1.2)):
            logging.warning("Warning! Some parameters have a Gelman-Rubin statistic greater than 1.2.")
    else:
        print("You need at least 2 chains to do convergence tests.")


def summarize_trace(trace, taxon_names = None, inclusion_cutoff=0.30):
    """
    Summarizes the trace object from run_svss()
    :param trace: trace object from the pymc3 run
    :param taxon_names: Optional list of taxon names to include in the spreadsheet
    :param inclusion_cutoff: fraction of samples a bug must be selected to consider it positive
    :return: dataframe with the inclusion probabilities, and regression coefficients for the included bacteria
    """
    inclusion_full, found_full = get_positives(trace, inclusion_cutoff)

    Df = pd.DataFrame({"taxon_num": found_full, "inclusion_p": inclusion_full[found_full],"beta_coeff": trace.get_values('beta_slab').mean(0)[found_full]})
    if taxon_names is not None:
        assert len(taxon_names) == len(inclusion_full)
        taxon_names = np.array(taxon_names)
        Df["taxon_name"] = taxon_names[found_full]

    return Df


def remove_null(X, Y):
    """
    If you have null values in Y, this removes them and the corresponding rows of X
    :param X: microbiome data
    :param Y: metabolite data
    :return: X, Y (subset)
    """
    Y = pd.Series(Y)
    locs_keep = ~pd.isnull(Y)
    locs_keep = np.where(locs_keep)[0]
    Y = Y[locs_keep]
    X = X.iloc[locs_keep, :]
    return X, Y


def threshold_count_data(df, min_counts = 1000, min_ab=1e-4):
    assert 'int' in str(df.iloc[:,0].dtype), "This is not counts data"

    if min_counts is not None:
        df = df.loc[df.sum(1) >= min_counts,:]

    # Convert to relative abundance
    df_rel = df.apply(lambda x: x / np.sum(x),1)
    pseudo = np.min(df_rel[df_rel > 0]) / 10.0

    mu = df_rel.apply(lambda x: np.mean(x), 0)
    df = df.loc[:,mu > min_ab]
    return df


def score(identified, actual_pos):
    tp = np.sum(np.isin(identified,actual_pos))
    fp = np.sum(~np.isin(identified,actual_pos))
    precision = tp / (tp + fp)
    recall = tp / len(actual_pos)
    f1 = 2 * precision * recall / (precision + recall)
    return precision, recall, f1


def run_bracod(X_prop: np.array, Y: np.array, n_sample: int = 1000, n_burn: int = 5000, g: float = 500.0, njobs: int = 2, tau_fixed: typing.Union[float] = None, mu_t: float = 5.0,sigma_t: float = 1.0, inclusion_prior: float = 0.25, init_method="auto") -> object:
    """
    Initilizes the model and uses MCMC sampling on the posterior
    :param X_prop: Relative abundance values for the microbiome
    :param Y: Metabolite concentration values
    :param n_sample: number of MCMC samples after the burn-in period
    :param n_burn: number of burn-in samples
    :param g: The scaling factor for how much larger the included coefficient distribution is
    :param njobs: Number of chains in the pymc model
    :param tau_fixed: if a float then this value if used for the coefficient variance, if None this is sampled from the model
    :param mu_t: prior mean of the total abundance distribution (log-normal)
    :param sigma_t: prior variance of the total abundance distribution (log-normal)
    :param inclusion_prior: prior parameter of the inclusion probability (bernoulli)
    :return: pymc3 trace object
    """
    n_samples = X_prop.shape[0]
    n_bugs = X_prop.shape[1]

    # Sometimes there can be issues with the R object not being a numpy array
    Y = np.array(Y)

    # warning about too many bugs
    if n_bugs >= 300:
        logging.warning("Warning! you have a lot of bacteria in here, did you threshold to the most abundant?")
    assert np.isnan(Y).sum() == 0, "You have nan values in your environmental variable"
    assert np.sum(X_prop == 0).sum() == 0, "You have 0 values in your OTU abundance, you need a pseudo count"
    assert np.allclose(X_prop.sum(1), 1, atol = 0.0001), "This is not relative abundance data with bacteria as rows and microbiomes as columns"
    assert Y.shape[0] == X_prop.shape[0], "Environmental variable must have the same number of samples as the OTU data"

    # Normalize
    Y = scale(Y)

    linear_model = pm.Model()
    with linear_model:
        # Define as a deterministic value in case we want to compare values after the fact
        X_Ab = pm.Deterministic('abs_ab', pm.math.log(X_prop))

        # Inclusion value
        p = pm.Bernoulli("p", inclusion_prior, shape=n_bugs)

        # Either set tau, or estimate it from the model
        if tau_fixed is None:
            tau = pm.HalfNormal('tau', sd=1)
        else:
            tau = tau_fixed

        # select the spike or slab, depending on the value of b
        # betas * b (zero or 1) is variable selection. It zeros out variables that aren't included
        betas_one = pm.Normal('betas_one', mu=0, sd=tau * g, shape=n_bugs, testval=0.)
        betas_zero = pm.Normal('betas_zero', mu=0, sd=tau, shape=n_bugs, testval=0.)
        betas_slab = pm.Deterministic('beta_slab', (1 - p) * betas_zero + p * betas_one)

        # Regression intercept
        alpha = pm.Normal("alpha", mu=Y.mean(), sd=100)

        # Mean of the metabolite values
        Y_hat = pm.Deterministic('predicted', alpha + pm.math.dot(X_Ab, betas_slab))
        # Variance of the metabolite values
        sigma_squared = pm.HalfNormal('sigma', 5)

        def logp(y):
            """
            Cutsom likelihood function
            See derivation in supplementary
            :param y: metabolite values
            :return: log likelihood
            """
            beta_sum = pm.math.sum(betas_slab)
            a = (pm.math.sqr(beta_sum)) / sigma_squared + 1 / sigma_t
            b = beta_sum * (Y - alpha - pm.math.dot(X_Ab, betas_slab)) / sigma_squared + mu_t / sigma_t
            c = pm.math.sum(pm.math.sqr(y - Y_hat)) / sigma_squared
            return 1 / 2 * (-n_samples * pm.math.log(sigma_squared) - n_samples * pm.math.log(a) + pm.math.sum(pm.math.sqr(
                b) / a) - c)

        # Include the custom likelihood function
        loglik = pm.DensityDist('loglik', logp, observed=dict(y=Y))

        # Sample
        trace = pm.sample(int(n_sample), tune=int(n_burn), chains=int(njobs), cores=int(njobs)*2, init=init_method)

    return trace


