## Unix Makefile to automate many IDL/MAS tasks.


WWW_SERVER=marecilab.mbi.ufl.edu
WWW_USER=marecilab
WWW_PATH=/home/marecilab/www/software/MAS/

INSTALLDIR=/export/faraday/thmlab/idl/mas/
INSTALLSRC=/export/faraday/thmlab/idl/mas_source/

ifndef $(IDL)
	IDL=idl71
endif

.PHONY: all
all: mas.sav

.PHONY: install
install: #all upload
	@test -e mas.sav && echo cp -av mas.sav $(INSTALLDIR)
	echo $(shell pwd)

.PHONY: upload
upload: mas.sav current_version.txt ChangeLog
	scp mas.sav ChangeLog current_version.txt \
	    $(WWW_USER)@$(WWW_SERVER):$(WWW_PATH)

.PHONY: mas.sav
mas.sav: mas.pro
	$(IDL) -e "@BUILD_SCRIPTS/EXPORT_MAS"

.PHONY: ChangeLog
ChangeLog: 
	git log --tags --pretty --decorate \
	      --summary --stat --name-only \
	      --no-merges --abbrev-commit > ChangeLog

.PHONY: current_version.txt
current_version.txt: mas.pro
	git tag | tail -1 | tr [a-zA-Z] 0 > current_version.txt

.PHONY: export
export: mas.sav ChangeLog

.PHONY: clean
clean:
	rm -f mas.sav ChangeLog current_version.txt

