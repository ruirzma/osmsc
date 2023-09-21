"""
OSMsc setup script.

See license in LICENSE.txt.
"""

from setuptools import setup

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

with open("requirements.txt") as f:
    requirements = [line.strip() for line in f.readlines()]

setup(name='osmsc',
    version='0.1.20',
    author='Rui Ma',
    author_email='rui.rz.ma@gmail.com',
    install_requires= requirements,
    description='Construct semantic city models from OpenStreetMap',
    long_description=readme,
    url='https://github.com/ruirzma/osmsc',
    license='MIT',
    packages=['osmsc'],
    zip_safe=False)

