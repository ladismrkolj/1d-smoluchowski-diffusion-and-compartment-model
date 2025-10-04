      Subroutine SD1D(Npoi,x,sd,L)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calculates real second derivatives sd() of the real function given at Npoi points x(),
c     interval is L long      
c     
c     INPUT :
c     Npoi = # of points
c     L = length of the interval
c     x() = function values 
c     OUTPUT :
c     sd() second derivative at dicrete points
c
c     Required routines:
c     FFT1D
c     FFT2C     
c     FFTRC
c     FFTCC
c
c     JM, Ljubljana, March 2022
c     Ladi Smrkolj, Ljubljana March 2022 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit None
      Include 'sizes.diff_cyl'
      Real*8 Zero, x(MaxBin), sd(MaxBin)
      Complex*16 A(MaxBin)
      Integer Npoi

      Real*8 Pii,L
      Complex*16 Wave
      Integer k,n,i

      Pii = 4.0d0*Atan(1.0d0)
      Zero = 0.0d0

      Do i=1,Npoi
         A(i)= Dcmplx(x(i), Zero)
      End Do

      Call FFT1D(A,Npoi,1)
      
      Do k=1,Npoi
         n = k - 1
         If (n.lt.Npoi/2) Then
            Wave = Cmplx(0.0d0, -2.0d0*Pii*Dble(n)/L)
         Else
            Wave = Cmplx(0.0d0, -2.0d0*Pii*Dble(-Npoi+n)/L)
         End If
         A(k) = A(k)*(Wave**2)
      End Do

      Call FFT1D(A,Npoi,-1)

      Do i=1,Npoi
         sd(i)= Dreal(A(i))
      End Do

      
      Return
      End




      
