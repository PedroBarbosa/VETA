import setuptools, shutil

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name='veta',
    version='0.5.1',
    licence='MIT',
    author='Pedro Barbosa',
    author_email='pedro.barbosa@medicina.ulisboa',
    description="Package that benchmarks variant prediction tools",
    url="https://github.com/PedroBarbosa/VETA",
    download_url="https://pypi.python.org/pypi/veta",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
    ],

    python_requires='>=3.6',    
    install_requires=[
        'matplotlib>=3.3.3',
        'numpy>=1.19.5',
        'pandas',
        'seaborn>=0.11.1',
        'cyvcf2>=0.30.4',
        'pysam>=0.16',
        'scikit-learn>=0.24.0',
        'imbalanced-learn>=0.7.0',
        'hgvs>=1.5.1',
        'fastcluster>=1.1.26',
        'statannot>=0.2.3',
        'gplearn>=0.4.1',
        'sklearn-pandas>=2.0.4',
        'dtreeviz>=1.1.3'],
    include_package_data=True,
    package_data={'src': ['config/tools_config.txt']},
    entry_points={
        'console_scripts': ['veta=src.veta:main']
    }
)
