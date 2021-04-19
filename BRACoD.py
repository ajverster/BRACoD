
import pymc3 as pm
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from scipy.stats import pearsonr

import typing
import logging

import theano
theano.config.blas.ldflags = '-lf77blas -latlas -lgfortran'

def scale_counts(df, min_counts = 1000, top_n_bugs = None):
    """
    Takes a DataFrame of OTU counts and normalizes it for use with run_bracod()
    :param df: DataFrame of OTU counts
    :param top_n_bugs: most abundant n bugs will be returned
    :return: A DataFrame of relative abundance data
    """
    df = pd.DataFrame(df)
    assert 'int' in str(df.iloc[:,0].dtype), "This is not counts data"
    if min_counts is not None:
        df = df.loc[:,df.sum(0) >= min_counts]
    # Convert to relative abundance
    df = df.apply(lambda x: x / np.sum(x),1)
    # Add a pseudo count
    df[df == 0] = np.min(df[df > 0]).min() / 2
    # Renorm
    df = df.apply(lambda x: x / np.sum(x),1)
    # Sort by descending abundance
    cols_sort = df.mean(0).sort_values(ascending=False).index
    assert all([x in df.columns for x in cols_sort])
    df = df.loc[:,cols_sort]

    if top_n_bugs is not None:
        # Subset to most abundance bugs
        top_n_bugs = int(top_n_bugs)
        df = df.iloc[:,:top_n_bugs]
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
    inclusion_perchain = np.zeros((trace.get_values('b').shape[1], n_chains))

    for i in range(n_chains):
        inclusion_perchain[:,i] = trace.get_values('b', chains=[i]).sum(0) / trace.get_values('b', chains=[i]).shape[0]
        
    # Your chains are radically different
    radically_different = np.where(np.apply_along_axis(lambda x: np.max(x) - np.min(x),1,inclusion_perchain) > 0.5)[0]
    if len(radically_different) > 0:
        return radically_different
    return None


def get_positives(trace, inclusion_cutoff=0.50):
    """
    Return locations of contributing OTUs
    :param trace:
    :param inclusion_cutoff:
    :return:
    """
    inclusion_full = trace.get_values('b').sum(0) / trace.get_values('b').shape[0]
    found_full = np.where(inclusion_full >= inclusion_cutoff)[0]
    return inclusion_full, found_full


def convergence_tests(trace, inclusion_cutoff=0.50):
    """
    This function runs a series of convergence tests for the important variables in Bannoc
    Meant to replace the less focused warnings that come from pm.sample
    :param trace:
    :param inclusion_cutoff:
    :return:
    """
    # Calculate the inclusion probabilities
    if trace.nchains >= 2:
        all_chains = list(range(trace.nchains))
        halfway = int(trace.nchains/2)
        inclusion_1 = trace.get_values('b', chains=all_chains[:halfway]).sum(0) / trace.get_values('b', chains=all_chains[:halfway]).shape[0]
        inclusion_2 = trace.get_values('b', chains=all_chains[halfway:]).sum(0) / trace.get_values('b', chains=all_chains[halfway:]).shape[0]
        found1 = set(np.where(inclusion_1 >= inclusion_cutoff)[0])
        found2 = set(np.where(inclusion_2 >= inclusion_cutoff)[0])
        found_uncertain = (found1 - found2).union(found2 - found1)

        
        if len(found_uncertain) > 0:
             logging.warning("Warning! The following bugs were found in only one of the chains: {}".format(
                " ".join([str(x) for x in found_uncertain])))

        diff = check_chains_equal(trace)
        if diff is not None:
            logging.warning("Warning! the following bugs are radically different between chains {}".format(" ".join([str(x) for x in diff])))


        # Check the effective number of samples for the positive bugs
        effn_b = pm.diagnostics.effective_n(trace, ["b"])["b"] 
        effn_betas_one = pm.diagnostics.effective_n(trace, ["betas_one"])["betas_one"] 
        gr_b = pm.diagnostics.gelman_rubin(trace, ["b"])["b"]
        gr_betas_one = pm.diagnostics.gelman_rubin(trace, ["betas_one"])["betas_one"]
        inclusion_full, pos_values = get_positives(trace)
        if (any(effn_b[pos_values] <= 200) | any(effn_betas_one[pos_values] <= 200)):
            logging.warning("Warning! Some parameters have an effective sample size less than 200.")

        if (any(gr_b[pos_values] >= 1.2) | any(gr_betas_one[pos_values] >= 1.2)):
            logging.warning("Warning! Some parameters have a Gelman-Rubin statistic greater than 1.2.")
    else:
        print("You need at least 2 chains to do convergence tests.")


def summarize_trace(trace, bug_names = None, inclusion_cutoff=0.50):
    """
    Summarizes the trace object from run_svss()
    :param trace: trace object from the pymc3 run
    :param outfile: Optional savefile for a graph of the consistency of the run
    :param inclusion_cutoff: fraction of samples a bug must be selected to consider it positive
    :return: dataframe with the inclusion probabilities, and regression coefficients for the included bugs
    """
    inclusion_full, found_full = get_positives(trace, inclusion_cutoff)

    Df = pd.DataFrame({"bugs": found_full, "inclusion_p": inclusion_full[found_full],"coefficients": trace.get_values('beta_slab').mean(0)[found_full]})
    if bug_names is not None:
        assert len(bug_names) == len(inclusion_full)
        bug_names = np.array(bug_names)
        Df["bug_names"] = bug_names[found_full]

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

    # warning about too many bugs
    if n_bugs >= 300:
        logging.warning("Warning! you have a lot of bugs in here, did you threshold to the most abundant?")
    assert np.isnan(Y).sum() == 0, "You have nan values in your environmental variable"
    assert np.allclose(X_prop.sum(1), 1, atol = 0.0001), "This is not relative abundance data with bugs as rows and microbiomes as columns"
    assert Y.shape[0] == X_prop.shape[0], "Environmental variable must have the same number of samples as the OTU data"

    # Normalize
    Y = scale(Y)

    linear_model = pm.Model()
    with linear_model:
        # Define as a deterministic value in case we want to compare values after the fact
        X_Ab = pm.Deterministic('abs_ab', pm.math.log(X_prop))

        # Inclusion value
        b = pm.Bernoulli("b", inclusion_prior, shape=n_bugs)

        # Either set tau, or estimate it from the model
        if tau_fixed is None:
            tau = pm.HalfNormal('tau', sd=1)
        else:
            tau = tau_fixed

        # select the spike or slab, depending on the value of b
        # betas * b (zero or 1) is variable selection. It zeros out variables that aren't included
        betas_one = pm.Normal('betas_one', mu=0, sd=tau * g, shape=n_bugs, testval=0.)
        betas_zero = pm.Normal('betas_zero', mu=0, sd=tau, shape=n_bugs, testval=0.)
        betas_slab = pm.Deterministic('beta_slab', (1 - b) * betas_zero + b * betas_one)

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


