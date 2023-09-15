from setuptools import Command, find_packages, setup
import sys

DESCRIPTION = "Fast and easy analysis of netCDF data in Python"
LONG_DESCRIPTION = """

**ecoval** is an automated ERSEM evaluation toolkit 


"""

PROJECT_URLS = {
    "Bug Tracker": "https://github.com/pmlmodelling/ecoval/issues",
    "Source Code": "https://github.com/pmlmodelling/ecoval",
}

extras_require: dict() = dict()


extras_require["complete"] = ["geoviews", "rioxarray", "cfchecker", "geocube", "geopandas"]

REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]

setup(name='ecoval',
      version='0.1.0',
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      python_requires='>=3.6.1',
      classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],

      project_urls=PROJECT_URLS,
      url = "https://github.com/pmlmodelling/nctoolkit",
      author='Robert Wilson',
      maintainer='Robert Wilson',
      author_email='rwi@pml.ac.uk',
      include_package_data=True,
      package_data={
      'ecoval': ['data/*'] },

      packages = ["ecoval"],
      setup_requires=[
        'setuptools',
        'setuptools-git',
        'wheel',
    ],
      install_requires = REQUIREMENTS,
      extras_require = extras_require,
      zip_safe=False)




