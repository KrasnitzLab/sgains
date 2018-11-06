import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sgains",
    version="1.0.0",
    author="Lubomir Chorbadjiev",
    author_email="lubomir.chorbadjiev@gmail.com",
    description="Sparse Genomic Analysis of Individual Nuclei by "
    "Sequencing (s-GAINS)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KrasnitzLab/sgains",
    # packages=setuptools.find_packages(
    #     # 'sgains',
    #     exclude=[
    #         'docs', 'tests', 
    #         'commands.tests',
    #         'pipeline.tests',
    #     ]
    # ),
    packages=['sgains', 'sgains.commands', 'sgains.pipelines'],
    package_dir={
        'sgains': 'sgains',
    },
    package_data={
        'sgains': ['scripts/*.R'],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'sgains-tools=sgains.cli:main',
        ]
    },
    classifiers=(
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    python_requires='>=3.6',
    
)
