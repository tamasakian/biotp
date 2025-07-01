from setuptools import setup, find_packages

setup(
    name="biotp",
    version="1.4.2",
    packages=find_packages(),

    author="tamasakian",
    description="A Bioinformatics Textfile Processor",

    install_requires=[
        "pandas",
        "biopython",
    ],
)
