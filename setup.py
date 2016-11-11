import os
import re
import codecs

from setuptools import setup, find_packages
# from pip.req import parse_requirements


# requirements = parse_requirements(
    # os.path.join(os.path.dirname(__file__), 'requirements.txt'),
    # session='hack'
# )
# INSTALL_REQUIRES = [str(req.req) for req in requirements] 

def find_version(*file_paths):
    # inspired by pip repository
    def read(*parts):
        path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            *parts)
        return codecs.open(path, 'r').read()

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
    long_description='',
    author='Zymo Research',
    author_email='mjin@zymoresearch.com',
    url='https://github.com/Zymo-Research/testprimer/tree/master',
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

