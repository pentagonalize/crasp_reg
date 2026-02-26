all: install test

install:
	sage -pip install --upgrade --no-index -v .
	pip install --upgrade --no-index -v .
test:
	sage -tp --force-lib pysemigroup/*py 

coverage:
	sage -coverage pysemigroup/*

diste: 
	python3 setup.py sdist

register: dist
	twine register dist/pysemigroup-$$VERSION.tar.gz

upload: diste
	 twine upload dist/pysemigroup-$$VERSION.tar.gz

