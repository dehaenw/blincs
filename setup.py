from setuptools import setup

setup(
    name='blincs',
    version='0.0.1alpha',    
    description='breadth first chemical line notation',
    url='https://github.com/dehaenw/blincs',
    author='Wim Dehaen',
    packages=['blincs'],
    install_requires=['rdkit>=2022.03'],
)
