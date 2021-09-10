from setuptools import setup, find_packages
import sys
import os

if sys.version_info.major != 3:
    raise RuntimeError("My_GSEA requires Python 3")

dir = os.path.dirname(os.path.abspath(__file__))
req_file = os.path.join(dir, 'requirements.txt')
with open(req_file) as f:
    required = f.read().splitlines()

setup(
    name="My_GSEA",
    description="Gene Set Enrichment Analysis",
    version="0.0",
    author="John-William Sidhom",
    author_email="johnwilliamsidhom@gmail.com",
    packages=find_packages(),
    install_requires = required,
    url="https://github.com/sidhomj/my_gsea",
    license="LICENSE",
    long_description=open(os.path.join(dir,"README.md")).read(),
    long_description_content_type='text/markdown',
    package_data={'my_gsea':[os.path.join('gene_sets','*')]}
)