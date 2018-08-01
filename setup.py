#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='nasp',
    version='1.1.2',
    description='Northern Arizona SNP Pipeline',
    long_description=readme + '\n\n' + history,
    author='Darrin Lemmer',
    author_email='dlemmer@tgen.org',
    url='https://github.com/TGenNorth/nasp',
    packages=[
        'nasp'
        #, 'nasp.vtm'
    ],
    package_dir={
        'nasp': 'nasp',
    },
    include_package_data=True,
    install_requires=requirements,
    license="Academic and Research License",
    zip_safe=False,
    keywords='nasp',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    test_suite='tests',
    tests_require=test_requirements,
    entry_points = {
        'console_scripts': [
            'nasp = nasp.nasp:main',
            'format_fasta = nasp.format_fasta:main',
            'find_duplicates = nasp.find_duplicates:main',
            'convert_external_genome = nasp.convert_external_genome:main',
            'vcf_to_matrix = nasp.vcf_to_matrix:main',
        ]
    },
    scripts = [
        'scripts/filter_matrix_by_coord.py',
        'scripts/filter_matrix_by_distance.py',
        'scripts/filter_matrix_by_genome.py',
        'scripts/matrix_to_fasta.py',
        'scripts/merge_matrices.py',
        'scripts/report_single_snps_single_isolate.py'
    ]
)
