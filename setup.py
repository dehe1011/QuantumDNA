from setuptools import setup, find_packages

setup(
    name='QuantumDNA', 
    version='1.0.0', 
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Dennis Herb',
    author_email='dennis.herb@uni-ulm.de',
    url='https://github.com/dehe1011/QuantumDNA', 
    packages=find_packages(),
    install_requires=[
        'numpy', 'matplotlib', 'qutip', 'pandas'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        # 'Operating System :: 'Windows', 'Linux', 'macOS'
    ],
    python_requires='>=3.6',
)

