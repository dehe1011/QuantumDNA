# Create a documentation on readthedocs.org

Install and launch sphinx:

```pip install sphinx```
```sphinx-quickstart```

* collect all files in the ```docs/``` directory
* adopt ```docs/index.rst```
* adopt ```docs/conf.py```
* if you want to use the sphinx_rtd_theme in conf.py you have to install it first:

```pip install sphinx_rtd_theme```

Create a local html file ```docs/_build/html/index.html```
```docs/make.bat html```

* create ```readthedocs.yaml```
* create ```docs/requirements.txt```

* make a new commit and push changes to GitHub 
* login to readthedocs.org
* refresh and select your GitHub project (repository must be public)
* start "Build-Version" to create the webpage