PACKAGE=shmm
VERSION=0.1
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

all:
	make doc-update
	make install
	make pdf
	make test

doc-update:
	echo "roxygen2::roxygenize('shmm/', roclets=c('rd', 'collate', 'namespace'))" | R --slave
#echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | R --slave	

build-package:
	R CMD build --resave-data=no $(PACKAGE)

install:
	make build-package
	R CMD INSTALL --preclean --no-multiarch $(TARBALL)
	date

unexport TEXINPUTS
pdf:
	rm -f $(PACKAGE).pdf
	R CMD Rd2pdf --no-preview $(PACKAGE)

check:
	R CMD check $(PACKAGE)

#test:
#	echo 'source("../work/production_model/spict/testing/src/maketest.R")' | R --vanilla

