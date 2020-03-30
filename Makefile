TST_NOTEBOOK:=algebraic_model.ipynb
#TST_IPNB:=algebraic_model_tst.ipynb

.PHONY: doc test

zip:
	zip lca-algebraic.zip algebraic_model.ipynb lca_algebraic.py doc/*

doc:
	pdoc3 --html lca_algebraic --force
	mv html/lca_algebraic/* doc/
	rm -r html

test:
	# ipyrmd -y --from Rmd --to ipynb $(TST_NOTEBOOK) -o $(TST_IPNB)
	./test/notebook_runner.py $(TST_NOTEBOOK)
