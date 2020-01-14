from setuptools import setup, find_packages


setup(
    name='kegg_query',
    version='0.0.1',
    author='Yue Zhang',
    author_email='yjzhang@cs.washington.edu',
    url='https://github.com/yjzhang/kegg_query',
    license='MIT',
    install_requires=['backports.functools_lru_cache'],
    packages=find_packages("."),
    # this is for including the data dir in the package.
    zip_safe=False,
    package_data={'kegg_query': ['data/kegg.db']},
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
    ],

)
