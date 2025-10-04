# Makefile for ode
#
# JM and DP, Ljubljana January 2018
#
FFLAGS = -fdefault-real-8 -g -Wall -fcheck=all
FC = gfortran
#
allobj =    \
rk4.o       \
main.o     
#
ode: $(allobj)
	$(FC) $(FFLAGS) -o ode $(allobj)

main.o  : sizes.ode
rk4     : sizes.ode
