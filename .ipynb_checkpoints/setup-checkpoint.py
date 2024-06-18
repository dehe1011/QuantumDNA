from setuptools import setup, find_packages

setup(
    name='Quantum_DNA', 
    version='1.0', 
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Dennis Herb',
    author_email='dennis.herb@uni-ulm.de',
    url='https://github.com/dehe1011/quantum_DNA', 
    packages=find_packages(),
    install_requires=[
        'numpy', 'matplotlib', 'qutip', 'pandas'
    ],
    tests_require=[
        'unittest',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

