CC = gcc
CFLAGS = -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c
LFLAGS = -static -pthread -L. -o 
SRCS = Turing.c CompTuring.c 
LIBS = -l:libhgt.a -l:libmpfr.a -l:libgmp.a
OBJS = $(SRCS:.c=.o)
DEPS = hgt.h turing.h
TARGET = turing

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(TARGET) $(OBJS) $(LIBS)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
