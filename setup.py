from os import name
from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name      = 'PyMM',
    version   = '1.0',
    long_description=readme,
    license=license,
    url='https://github.com/ChenGiuseppe/PyMM',
    packages  = find_packages(),
    entry_points= {
        'console_scripts': ['pymm=pymm.pymm:main']
    }
)
