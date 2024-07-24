from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="qDNA",
    version="0.1.0",
    author="Dennis Herb",
    author_email="dennis.herb@uni-ulm.de",
    description="A package to calculate lifetimes, average charge separation and dipole moments of excited states along DNA within the formalism of open quantum systems.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/dehe1011/QuantumDNA",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=required,
)
