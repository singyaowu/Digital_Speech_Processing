.PHONY: all clean

CFLAGS+= -g
LDFLAGS+=-lm -g    # link to math library

TARGET=train test

all: $(TARGET)
# type make/make all to compile test_hmm
	echo $(CFLAGS)
clean:
	$(RM) $(TARGET)   # type make clean to remove the compiled file
