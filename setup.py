from setuptools import setup, find_packages

with open("configs/requirements/requirements.txt") as f:
    required = f.read().splitlines()

with open("configs/dev-requirements.txt") as f:
    dev = f.read().splitlines()

with open("configs/docs-requirements.txt") as f:
    docs = f.read().splitlines()

extras = dev + docs

setup(
    name="qDNA",
    version="0.1.1",
    author="Dennis Herb",
    author_email="dennis.herb@uni-ulm.de",
    description="A package to calculate lifetimes, average charge separation and dipole moments of excited states along DNA within the formalism of open quantum systems.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/dehe1011/QuantumDNA",
    license="BSD-3-Clause",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    install_requires=required,
    extras_require=extras,
)
