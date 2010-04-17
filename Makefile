CFLAGS = -std=c99 -DPD -O3 -Wall -W -Wshadow -Wstrict-prototypes -Wno-unused -Wno-parentheses -Wno-switch -funroll-loops -pipe -fomit-frame-pointer -ffast-math


UNIVERSAL=-arch i386 -arch ppc
DARWINCFLAGS = $(CFLAGS) -DDARWIN $(UNIVERSAL)
DARWIN_LIBS=$(UNIVERSAL)
NTCFLAGS = $(CFLAGS) -DNT
WIN_CC=gcc.exe
WIN_STRIP=strip.exe

linux: bandlimited~.c
	gcc $(CFLAGS) -o bandlimited~.o -c bandlimited~.c
	ld -export_dynamics -shared -o bandlimited~.pd_linux bandlimited~.o
	strip --strip-unneeded bandlimited~.pd_linux

darwin: bandlimited~.c
	  cc  $(DARWINCFLAGS) -pedantic -o bandlimited~.o -c bandlimited~.c
		cc -bundle -undefined suppress -flat_namespace $(DARWIN_LIBS) -o bandlimited~.pd_darwin bandlimited~.o

win32: bandlimited~.c
	${WIN_CC} $(NTCFLAGS) -o bandlimited~.o  -c  bandlimited~.c   
	${WIN_CC} $(NTCFLAGS) -Lc:/Arquivos\ de\ programas/pd/bin/ -lpd  -shared -o bandlimited~.dll   bandlimited~.o -W1  
 	#${WIN_STRIP} --strip-unneeded bandlimited~.dll

clean:
	rm *.o
	rm bandlimited~.pd*
