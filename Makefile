TST_NOTEBOOK:=example-notebook.ipynb
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
	conda build conda-recipe/meta.yaml

package : package-pip package-conda

tst-upload-pip:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/lca_algebraic*

upload-pip:
	twine upload -u oie-minesparistech dist/lca_algebraic*

upload-conda:
	anaconda upload --force ~/anaconda3/conda-bld/noarch/lca_algebraic-$(VERSION)-py_0.tar.bz2

upload : upload-pip upload-conda

pytest:
	pytest test/tests.py

notebook_test:
	# ipyrmd -y --from Rmd --to ipynb $(TST_NOTEBOOK) -o $(TST_IPNB)
	./test/notebook_runner.py $(TST_NOTEBOOK)

test: pytext notebook_test
