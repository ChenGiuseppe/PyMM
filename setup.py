from os import name
from setuptools import setup, find_packages

setup(
    name      = 'pmm',
    version   = '0.0.2a1',
    packages  = find_packages(),
    entry_points= {
        'console_scripts': ['pmm=pmm.pmm:main']
    }
)