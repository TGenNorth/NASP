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
    version='0.9.7',
    description='Northern Arizona SNP Pipeline',
    long_description=readme + '\n\n' + history,
    author='Darrin Lemmer',
    author_email='dlemmer@tgen.org',
    url='https://github.com/TGenNorth/nasp',
    packages=[
        'nasp',
    ],
    package_dir={'nasp':
                 'nasp'},
    include_package_data=True,
    install_requires=requirements,
    license="Academic and Research License",
    zip_safe=False,
    keywords='nasp',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    entry_points = {
        'console_scripts': [
            'nasp = nasp.nasp:main',
            'format_fasta = nasp.format_fasta:main',
            'find_duplicates = nasp.find_duplicates:main',
            'convert_external_genome = nasp.convert_external_genome:main',
            'vcf_to_matrix = nasp.vcf_to_matrix:main'
        ]
    },
    scripts = [
        'scripts/filter_matrix_by_coord.py',
        'scripts/filter_matrix_by_genome.py',
        'scripts/filter_matrix_remove_genomes.py',
        'scripts/matrix_to_fasta.py',
        'scripts/report_single_snps_single_isolate.py'
    ]
)
