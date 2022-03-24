import setuptools

__version__ = '0.3.6'
__author__ = ['Adrian Verster']
__email__ = 'adrian.verster@hc-sc.gc.ca'

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BRACoD",
    install_requires=['pymc3==3.9.0',
                      'pandas>=0.24.0',
                      'numpy>=1.15,<=1.20',
                      'scikit-learn>=0.20',
                      'arviz<=0.10',
                      'Theano>=1.0.5'
                      ],
    python_requires='>3.6',
    description="BRACoD is a method to identify associations between bacteria and physiological variables in Microbiome data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ajverster/BRACoD/tree/main",
    package_dir={"": "src"},
    package_data= {"BRACoD": ["data/OTUCounts_obesitystudy.csv","data/SCFA_obesitystudy.csv", "data/Taxonomy.csv"]},
    packages=setuptools.find_packages(where="src"),
    version=__version__,
    author=__author__,
    author_email=__email__,
)
