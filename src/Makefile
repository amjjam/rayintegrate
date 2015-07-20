libdir=$(prefix)/lib

CPPFLAGS=-Wall -g -fPIC
CPP=g++

RAYINTEGRATE_SRC=rayintegrate.C
RAYINTEGRATE_OBJ=$(RAYINTEGRATE_SRC:.C=.o)

build: libRayintegrate.a libRayintegrate.so

libRayintegrate.a: $(RAYINTEGRATE_OBJ)
	ar rc $@ $^

libRayintegrate.so: $(RAYINTEGRATE_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^

install: $(libdir)/libRayintegrate.a $(libdir)/libRayintegrate.so

$(libdir)/libRayintegrate.a $(libdir)/libRayintegrate.so: libRayintegrate.a libRayintegrate.so
	cp libRayintegrate.a libRayintegrate.so $(libdir)

uninstall:
	- rm -f $(libdir)/libRayintegrate.a $(libdir)/libRayintegrate.so

clean:
	- rm -f libRayintegrate.a libRayintegrate.so $(RAYINTEGRATE_OBJ)
