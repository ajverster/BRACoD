from BRACoD import *
from BRACoD.BRACoD import *
from BRACoD.simulations import *
import pandas as pd


import pkg_resources

# Could be any dot-separated package/module name or a "Requirement"
resource_package = __name__
try:
    resource_path_otu = '/'.join(('data', 'OTUCounts_obesitystudy.csv'))  # Do not use os.path.join()
    infile_otu_data = pkg_resources.resource_stream(resource_package, resource_path_otu)
    resource_path_scfa = '/'.join(('data', 'SCFA_obesitystudy.csv'))  # Do not use os.path.join()
    infile_scfa_data = pkg_resources.resource_stream(resource_package, resource_path_scfa)
except:
    resource_path_otu = '/'.join(('../data', 'OTUCounts_obesitystudy.csv'))  # Do not use os.path.join()
    infile_otu_data = pkg_resources.resource_stream(resource_package, resource_path_otu)
    resource_path_scfa = '/'.join(('../data', 'SCFA_obesitystudy.csv'))  # Do not use os.path.join()
    infile_scfa_data = pkg_resources.resource_stream(resource_package, resource_path_scfa)

df_counts_obesity = pd.read_csv(infile_otu_data).T
df_scfa_obesity = pd.read_csv(infile_scfa_data)
