TST_NOTEBOOK:=example-notebook.Rmd
TST_IPYNB:=example-notebook.ipynb

VERSION:=$(shell cat VERSION)


.PHONY: doc test package tst-upload 

zip:
	zip lca-algebraic.zip algebraic_model.ipynb lca_algebraic.py doc/*

doc:
	pdoc3 --html lca_algebraic --force
	mv html/lca_algebraic/* doc/
	rm -r html

clean:
	rm -r dist

package-pip:
	python setup.py sdist bdist_wheel --universal

package-conda :
	conda build --py 3.9 conda-recipe

package : package-pip package-conda

tst-upload-pip:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/lca_algebraic*

upload-pip:
	twine upload dist/lca_algebraic*

upload-conda:
	anaconda upload --force ~/mambaforge/conda-bld/noarch/lca_algebraic-$(VERSION)-py_0.tar.bz2

upload : upload-pip upload-conda

pytest:
	pytest test

test_notebook:
	# ipyrmd -y --from Rmd --to ipynb $(TST_NOTEBOOK) -o $(TST_IPYNB)
	pytest --nbmake *.ipynb

lint:
	pre-commit run

lint-all:
	pre-commit run --all-files

test: pytest test_notebook

tox:
	tox run
