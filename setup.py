from setuptools import setup

ENTRY_POINTS={
        'console_scripts': [
            'MTWSPy_main = MTWSPy.MTWSPy_main:main',
            'mk_events_csv = MTWSPy.mk_events_csv:main',
            'find_twin_obs = find_twin_obs.find_twin_obs:main',
            'find_twin_syn = find_twin_syn.find_twin_syn:main',
            'match_twin_files = match_twin_files.match_twin_files:main',
            'phase_association_obs = phase_association_obs.phase_association_obs:main',
            'phase_association_syn = phase_association_syn.phase_association_syn:main',
            'correlate_twin = correlate_twin.correlate_twin:main',
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
    package_data={'MTWSPy': ['data/gcmt_catalog.csv']},
    python_requires='>=3.5',
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'obspy',
        'PyYAML',
        'pytest',
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

