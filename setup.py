from setuptools import setup

ENTRY_POINTS={
        'console_scripts': [
            'MTWSPy_main = MTWSPy.MTWSPy_main:main',
            'match_catalog = MTWSPy.match_catalog:main',
            'find_twin = MTWSPy.find_twin:main',
            'match_twin_files = MTWSPy.match_twin_files:main',
            'phase_association = MTWSPy.phase_association:main',
            'correlate_twin = MTWSPy.correlate_twin:main',
            'process_tdl_files = MTWSPy.post_processing.process_tdl_files:main',
        ]}

setup(
    name='MTWSPy',
    version='1.0.0',
    author='Alistair Boyce',
    maintainer='Alistair Boyce',
    author_email='alistair.boyce@univ-lyon1.fr',
    description='Python implementation of the Morphological Time Window Selection (MTWS)',
    license='MIT',
    url='https://github.com/alistairboyce11/MTWSPy',
    packages=['MTWSPy', 'MTWSPy.post_processing',
                'utils', 'utils.gcmt', 'utils.EarthScope_FetchData', 'utils.specfem'],
    package_dir={'MTWSPy': 'MTWSPy'},
    package_data={'MTWSPy': ['data/gcmt_catalog.csv', 
                             'data/obs/e2008/py_formatted/20080101063232/20080101063232_II_*T',
                             'data/syn/e2008/py_formatted/20080101063232/20080101063232_II_*T']},
    python_requires='>=3.5',
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'geographiclib',
        'obspy',
        'PyYAML',
        'pytest',
        'ellipticipy',
    ],
    entry_points=ENTRY_POINTS,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords='example documentation tutorial',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)

