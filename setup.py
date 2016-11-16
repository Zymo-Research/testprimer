import os
import re
import codecs

from setuptools import setup, find_packages


def read(*parts):
    # from pip setup.py
    here = os.path.abspath(os.path.dirname(__file__))
    return codecs.open(os.path.join(here, *parts), 'r').read()


def find_version(*file_paths):
    # from pip setup.py
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name='testprimer',
    version=find_version('testprimer', '__init__.py'),
    description='The tool for evaluating primer pools performance.',
    long_description=read('README.md'),
    author='Zymo Research',
    author_email='mjin@zymoresearch.com',
    url='http://www.zymoresearch.com/',
    packages=find_packages(exclude=['tests']),
    entry_points={
        'console_scripts': [
            'testprimer=testprimer.cli:main',
        ],
    },
    install_requires=[
        'pandas',
        'biopython',
        'openpyxl',
    ],
)
