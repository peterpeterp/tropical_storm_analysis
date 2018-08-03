
from setuptools import setup

setup(
    name='tropical_storm_analysis',
    version='0.1.0',
    description=(
        'Library .'
    ),
    long_description=open('README.md').read(),
    author='Peter Pfleiderer',
    author_email='',
    license='BSD 3-Clause',
    url='None',
    packages=['feature_tracking'],
    package_data={'': ['*.rst', '*.txt']},
    test_suite='tests',
    classifiers=[
        'Programming Language :: Python',
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        'Programming Language :: Python :: Implementation :: PyPy',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
    ],
)
