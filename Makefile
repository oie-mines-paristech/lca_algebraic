TST_NOTEBOOK:=example-notebook.Rmd
TST_IPYNB:=example-notebook.ipynb
DIST:=dist

# `twine check` in `package` needs twine>=6.1: isolated builds pull a recent
# setuptools (>=77, PEP 639) that emits Metadata-Version 2.4, which older twine
# rejects with the misleading "Metadata is missing required fields: Name, Version".
.PHONY: doc test package check clean tst-upload upload upload-tst

zip:
	zip lca-algebraic.zip algebraic_model.ipynb lca_algebraic.py doc/*

doc:
	cd doc && $(MAKE) html

clean:
	rm -rf $(DIST)

package: clean
	python -m build
	$(MAKE) check

check:
	twine check --strict $(DIST)/*

tst-upload: check
	twine upload --repository-url https://test.pypi.org/legacy/ $(DIST)/*

upload: check
	twine upload --repository-url https://upload.pypi.org/legacy/ $(DIST)/*

upload-tst: check
	twine upload --repository testpypi $(DIST)/*

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
