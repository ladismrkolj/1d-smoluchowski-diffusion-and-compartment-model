# Makefile for diff_cyl 
# JM, Ljubljana , March 2022 
#
FFLAGS = -O2
FC = gfortran
#
diffuseobj = main.o   \
             fd1d.o   \
             sd1d.o   \
            fft1d.o   \
            fftcc.o   \
            fftrc.o   \
            fft2c.o
#
diff_cyl: $(diffuseobj)
	$(FC)  $(FFLAGS) -o diff_cyl $(diffuseobj)

main.o   : sizes.diff_cyl
fd1d.o   : sizes.diff_cyl
sd1d.o   : sizes.diff_cyl
fft1d.o  : sizes.diff_cyl





