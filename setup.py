from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), 'r', encoding='utf-8') as f:
    long_description = f.read()

DESCRIPTION = "Perform permutation-based pathway enrichment analysis"
LONG_DESCRIPTION = long_description,
LONG_DESCRIPTION_CONTENT_TYPE = 'text/markdown'
NAME = "pathwayenrichment"
AUTHOR = "Semidán Robaina Estévez"
AUTHOR_EMAIL = "srobaina@ull.edu.es"
MAINTAINER = "Semidán Robaina Estévez"
MAINTAINER_EMAIL = "srobaina@gmail.com"
DOWNLOAD_URL = 'http://github.com/robaina/pathwayEnrichment'
LICENSE = 'Creative Commons Attribution 4.0 International'
VERSION = '0.0.3'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=['pathwayenrichment'],
      install_requires=['numpy', 'pandas', 'matplotlib', 'requests', 'scikit_learn']
      )
