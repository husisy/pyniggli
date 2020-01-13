from setuptools import setup, find_packages

setup(
    name='pyniggli',
    version='0.1',
    packages=find_packages(),
    description='niggli algorithm to reduce lattice',
    requires=['numpy'],
)
