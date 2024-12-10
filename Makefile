TST_NOTEBOOK:=example-notebook.Rmd
TST_IPYNB:=example-notebook.ipynb

VERSION:=$(shell cat VERSION)


.PHONY: doc test package tst-upload 

zip:
	zip lca-algebraic.zip algebraic_model.ipynb lca_algebraic.py doc/*

doc:
	cd doc && $(MAKE) html

clean:
	rm -r dist

package:
	python setup.py sdist bdist_wheel --universal

tst-upload:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/lca_algebraic*

upload:
	twine upload dist/lca_algebraic*


pytest:
	pytest test

test_notebook:
	# ipyrmd -y --from Rmd --to ipynb $(TST_NOTEBOOK) -o $(TST_IPYNB)
	pytest --nbmake notebooks/*.ipynb

lint:
	pre-commit run

lint-all:
	pre-commit run --all-files

test: pytest test_notebook

tox:
	tox run
