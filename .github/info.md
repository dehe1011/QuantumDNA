# Workflows

This file describes some workflows that are used frequently. 

# Testing with pytest

1. Install pytest package:
```pip install pytest```

2. Run pytest:
Windows Powershell: ```python -m pytest```
Conda Powershell: ```pytest```
All files of the form ```test_*.py``` in all folders are automatically detected. 

Alternative (to specify a directory and disable warnings):
```python -m pytest -vv tests/ --disable-pytest-warnings```


# Linting with pylint

1. Install pylint package:
```pip install pylint```

2. Run pylint:
Windows Powershell: ```python -m pylint qDNA/**/*.py```
Conda Powershell: ```pylint qDNA/**/*.py```


# Formatting with black 

1. Install black package:
```pip install black```

2. Run black:
```cd .. ```
Windows Powershell: ```python -m black QuantumDNA```
Conda Powershell: ````black QuantumDNA```
reformats all .py files inside the whole project. Alternatively you can apply black to selected folders and files. 


# Coverage with coveralls.io

1. Install coverage, coveralls and pytest packages:
```pip install coverage coveralls pytest```

2. Run coverage:
Windows Powershell: ```python -m coverage run -m pytest```
Conda Powershell: ```coverage run -m pytest```

The command creates a .coverage file in the working directory.
You can either print the report to the console (```(python -m) coverage report -m```) or save a .html file in the working directory (```(python -m) coverage html```).

3. Upload report to coveralls.io:

```$env:COVERALLS_REPO_TOKEN = "secret" ```
Windows Powershell: ```python -m coveralls```
Conda Powershell: ```coveralls```


# Documentation on readthedocs.org

1. Install sphinx package:
```pip install sphinx```

2. Start sphinx to automatically create all relevant files:
```sphinx-quickstart```
Collect all files in the ```docs/``` directory. 

3. Customize documentation:

* adopt ```docs/index.rst```
* adopt ```docs/conf.py```
* if you want to use the sphinx_rtd_theme in conf.py you have to install it first (```pip install sphinx_rtd_theme```)

4. Create a local html file:
```docs/make.bat html```
The file is located at ```docs/_build/html/index.html```.

5. Upload documentation to readthedocs.org:

* create ```readthedocs.yaml```
* create ```docs/requirements.txt``` (contains all packages needed to build the documentation)
* make a new commit and push changes to GitHub 
* login to readthedocs.org
* refresh and select your GitHub project (repository must be public)
* start "Build-Version" to create the webpage


# Publish package on pypi.org

1. Install setuptools, wheel and twine packages:
```pip install setuptools wheel twine```

2. Create source distribution and binary distribution from the setup.py file
```python setup.py sdist bdist_wheel ```

Creates a ````dist/``` and ```.egg-info``` directory. 

The package is saved locally in the dist/ directory and can be installed via
```pip install dist/qDNA-0.1.1-py3-none-any.whl``` (from the binary distribution)
```pip install dist/qDNA-0.1.1.tar.gz``` (from the source distribution)

3. Upload to pypi.org:

```$env:TWINE_USERNAME="__token__"```
```$env:TWINE_PASSWORD="secret"```
Windows Powershell ```python -m twine upload dist/*```
Conda Powershell ```twine upload dist/*```

Uploads sdist and bdist to pypi.