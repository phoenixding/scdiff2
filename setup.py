#!/usr/bin/env python
from setuptools import setup

setup(  name='scdiff2',
		version='1.0.0',
		description='Single Cell Differentiation Model package 2.0',
		author='Jun Ding',
		author_email='jund@andrew.cmu.edu',
		url="https://github.com/phoenixding/scdiff2",
		license='MIT',
		packages=['scdiff2'],
		package_data={'scdiff2':['tfdata/HumanTFList.txt']},
		entry_points={'console_scripts':['scdiff2=scdiff2.scdiff2:main','prerun=scdiff2.prerun:main']},
		install_requires=['scipy>0.13.3','numpy>1.8.2','scikit-learn>=0.20,<=0.22','matplotlib>=3.1.2','imbalanced_learn<0.5.0','anndata>=0.7','scanpy>=1.5','pandas>=0.23','h5py>=2.10','python-igraph>=0.8.0','leidenalg>=0.8.0'],
		classifiers=[
			'License :: OSI Approved :: MIT License',
			'Programming Language :: Python :: 3.6',
		],
		)
		
