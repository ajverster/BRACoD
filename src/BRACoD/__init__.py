from BRACoD import *
from BRACoD.BRACoD import *
from BRACoD.simulations import *
import pandas as pd


import pkg_resources

# Could be any dot-separated package/module name or a "Requirement"
resource_package = __name__
try:
    resource_path = '/'.join(('data', 'OTUCounts_obesitystudy.csv'))  # Do not use os.path.join()
    template = pkg_resources.resource_stream(resource_package, resource_path)
except:
    resource_path = '/'.join(('../data', 'OTUCounts_obesitystudy.csv'))  # Do not use os.path.join()
    template = pkg_resources.resource_stream(resource_package, resource_path)

#try:
#    import importlib.resources as pkg_resources
#except ImportError:
#    # Try backported to PY<37 `importlib_resources`.
#    import importlib_resources as pkg_resources


#from . import Data
#f = pkg_resources.open_text(Data,"OTUCounts_obesitystudy.csv")
example_otu_data = pd.read_csv(template).T
