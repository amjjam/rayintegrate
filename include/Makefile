includedir=$(prefix)/include

build:

install: $(includedir)/rayintegrate.H

$(includedir)/rayintegrate.H: rayintegrate.H
	cp rayintegrate.H $(includedir)

uninstall:
	- rm -f $(includedir)/rayintegrate.H

clean:
