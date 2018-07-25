# Makefile

include Make_include

LIBS    += -L$L -lgenfft -lm -lfftw3f $(LIBSM)
#OPTC += -g -O0 -Wall 

ALL: marchen 

SRCH	= newmarchenko.c \
		getFileInfo.c  \
		readData.c \
		new.readShotData.c \
		readTinvData.c \
		applyMute.c \
		writeData.c \
		writeDataIter.c \
		wallclock_time.c \
		name_ext.c  \
		verbosepkg.c  \
		atopkge.c \
		getpars.c \
		docpkge.c \
		pml.c

OBJH	= $(SRCH:%.c=%.o)

marchen:	$(OBJH) 
	$(CC) $(LDFLAGS) $(OPTC) $(CFLAGS) -o marchen $(OBJH) $(LIBS)

install: marchen 
	cp marchen $B


clean:
		rm -f core $(OBJH)

realclean: clean
		rm -f $B/marchen




