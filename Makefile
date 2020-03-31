TST_NOTEBOOK:=algebraic_model.ipynb
#TST_IPNB:=algebraic_model_tst.ipynb

.PHONY: doc test package tst-upload 

zip:
	zip lca-algebraic.zip algebraic_model.ipynb lca_algebraic.py doc/*

doc:
	pdoc3 --html lca_algebraic --force
	mv html/lca_algebraic/* doc/
	rm -r html

clean:
	rm -r dist

package:
	python setup.py sdist bdist_wheel --universal

tst-upload:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/lca_algebraic*

upload:
	twine upload -u oie-minesparistech dist/lca_algebraic*

test:
	# ipyrmd -y --from Rmd --to ipynb $(TST_NOTEBOOK) -o $(TST_IPNB)
	./test/notebook_runner.py $(TST_NOTEBOOK)
