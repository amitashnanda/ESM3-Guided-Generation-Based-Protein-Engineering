# In setup.py

from setuptools import setup, find_packages

setup(
    name="esm_foldx_guidedgeneration",
    version="0.1.0",
    author="Amitash Nanda",
    author_email="amitashnanda01@gmail.com",
    description="A tool for Guided Generation based Protein Design and Engineering using FoldX, or any other properties adding to custom-scoring function.",
    long_description=open('README.md').read(), 
    long_description_content_type="text/markdown",
    url="https://github.com/amitashnanda/ESM3-Guided-Generation-Based-Protein-Engineering.git", 
    
    package_dir={"": "src"}, 
    packages=find_packages(where="src"), 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        'torch',
        'pandas',
        'seaborn',
        'matplotlib',
        'esm', 
        'biopython',
        'huggingface-hub',
        'tqdm'
    ],

    entry_points={
        'console_scripts': [
            'esm_foldx_guidedgeneration = esm_foldx_guidedgeneration.main:main',
        ],
    },
)