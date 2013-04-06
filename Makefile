CFLAGS = -std=c99 -DPD -O3 -Wall -W -Wshadow -Wstrict-prototypes -Wno-unused -Wno-parentheses -Wno-switch -funroll-loops -pipe -fomit-frame-pointer -ffast-math


UNIVERSAL=-arch i386 -arch x86_64
DARWINCFLAGS = $(CFLAGS) -DDARWIN $(UNIVERSAL)
DARWIN_LIBS=$(UNIVERSAL)
NTCFLAGS = $(CFLAGS) -DNT
WIN_CC=gcc.exe
WIN_STRIP=strip.exe

linux: bandlimited~.c
	gcc $(CFLAGS) -o bandlimited~.o -c bandlimited~.c
	gcc $(CFLAGS) -o bandlimited_util.o -c bandlimited_util.c
	ld -export_dynamics -shared -o bandlimited~.pd_linux bandlimited_util.o bandlimited~.o
	strip --strip-unneeded bandlimited~.pd_linux

darwin: bandlimited~.c
	  cc  $(DARWINCFLAGS) -pedantic -o bandlimited~.o -c bandlimited~.c
	  cc  $(DARWINCFLAGS) -pedantic -o bandlimited_util.o -c bandlimited_util.c
		cc -bundle -undefined suppress -flat_namespace $(DARWIN_LIBS) -o bandlimited~.pd_darwin bandlimited_util.o bandlimited~.o 

win32: bandlimited~.c
	${WIN_CC} $(NTCFLAGS) -o bandlimited~.o  -c  bandlimited~.c   
	${WIN_CC} $(NTCFLAGS) -o bandlimited_util.o  -c  bandlimited_util.c   
	${WIN_CC} $(NTCFLAGS) -LC:/Program\ Files/pd/bin -lpd  -shared -o bandlimited~.dll  bandlimited_util.o bandlimited~.o -W1  
 	#${WIN_STRIP} --strip-unneeded bandlimited~.dll

clean:
	rm *.o
	rm bandlimited~.pd*
