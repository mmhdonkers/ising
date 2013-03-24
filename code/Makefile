
############################
#### Compiler and Flags ####
############################

# Compiler
FC = gfortran

FFLAGS = -Wall -Wextra -march=native -O3

# Compiler flags
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LDFLAGS = 
LIBS = $(shell pkg-config --libs plplotd-f95)



#####################
### General Rules ###
#####################

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)


# Files (modules) to link to main program
OBJS = 
OBJS += metrop.o
OBJS += plot.o
OBJS += ising.o
OBJS += wolff.o

all: run

run: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

clean: 
	$(RM) run $(OBJS) *.mod
