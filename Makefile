# prepare the package for release
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: docs build install clean 

full: deps docs vignettes check install clean

deps:
	R -e 'for(pkg in c("Matrix","fields","deSolve","devtools","stats")) if (!require(pkg,character.only=TRUE)) install.packages(pkg, repos="http://cran.rstudio.com")'

docs:
	R -e 'library(devtools); document()'

build:
	cd ..;\
	R CMD build --no-manual --no-build-vignettes $(PKGSRC)

build-cran:
	cd ..;\
	R CMD build $(PKGSRC)

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build-cran
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

vignettes: build vignettes/*Rmd
	R -e 'require(devtools); build_vignettes()'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

