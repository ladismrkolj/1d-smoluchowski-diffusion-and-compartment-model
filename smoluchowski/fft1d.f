      Subroutine FFT1D(A,N,IJob)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Performs FFT in 1D Complex->Complex (IJob=1)
c
c     or Inverse FFT  in 1D               (IJob=-1)
c
c     Input/Output : A() initial complex array in which also 
c                 the results are stored
c
c             N = number of points: 
c             it is recommended to be power of 2 
c             Ijob  (sse above)
c     Routines required : FFTCC (IMSL)
c    
c     JM, Ljubljana, January 1999
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit None
      Include 'sizes.diff_cyl'
      Complex*16 A(max)
      Integer N,IJob
      Logical Debug
c     Working arrays for FFTCC
      Real*8 Wk(max)
      Integer Iwk(max),i
      Data Debug /.False./

      If (Debug) Then
         Write(6,*) 'A() in fft1d before FFT'
         Do i=1,N
            Write(6,*) i, A(i)
         End Do
      End If
      
      If (IJob.eq.1) Then
         Call FFTCC(A,N,Iwk,Wk)
         Do i=1,N
            A(i) = A(i)/N
         End Do
      Else
c     Inverse FFT
         Do i=1,N
            A(i) = Conjg(A(i))
         End Do
         Call FFTCC(A,N,Iwk,Wk)
         Do i=1,N
            A(i) = Conjg(A(i))
         End Do
      End If

      If (Debug) Then
         Write(6,*) 'A() in fft1d after FFT'
         Do i=1,N
            Write(6,*) i, A(i)
         End Do
      End if
      
      Return
      End

      
