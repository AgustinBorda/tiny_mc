# Compilers
CC = gcc

# Flags
CFLAGS = -std=c11 -Wall -Wextra -Ofast -march=native -fopenmp -ftree-vectorize #-fopt-info-vec-missed -ffast-math
LDFLAGS = -lm

# Binary file
TARGET = tiny_mc

# Files
C_SOURCES = tiny_mc.c wtime.c
C_OBJS = $(patsubst %.c, %.o, $(C_SOURCES))

# Rules
all: $(TARGET)

$(TARGET): $(C_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) 


clean:
	rm -f $(TARGET) *.o

