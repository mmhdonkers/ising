FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O3
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LDFLAGS = 
LIBS = $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS = 
OBJS += ising.o

all: run

run: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

clean: 
	$(RM) run $(OBJS) *.mod
