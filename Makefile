# Makefile for plastid
#
# Run in project folder, not source folder.
date := $(shell date +%Y-%m-%d)

help:
	@echo "plastid make help"
	@echo " "
	@echo "Please use \`make <target>\`, choosing <target> from the following:"
	@echo "    dist        to make HTML documentation and eggs for distribution"
	@echo "    docs        to make HTML documentation"
	@echo "    cleandoc    to remove previous generated documentation components"
	@echo "    clean       to remove everything previously built"
	@echo " "

docs/source/class_substitutions.txt :
	mkdir -p docs/source
	docs/bin/get_class_substitutions plastid plastid
	mv plastid_substitutions.txt docs/source/class_substitutions.txt

docs/build/html : docs/source/class_substitutions.txt | docs/source/generated
	cp CHANGES.rst docs/source
	$(MAKE) html -C docs

docs/source/generated :
	sphinx-apidoc -M -e -o docs/source/generated plastid
	#rm docs/source/generated/plastid.test*rst

docs : | docs/build/html docs/source/generated

cleandoc : 
	rm -rf docs/build
	rm -rf docs/source/generated
	rm -rf docs/class_substitutions.txt

clean : cleandoc
	rm -rf dist
	rm -rf build
	
.PHONY : docs dist clean cleandoc help
