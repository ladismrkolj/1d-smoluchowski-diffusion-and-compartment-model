      Program Mju_Diffuse_Cyl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Time and coordinate dependent concentrations of local anesthetics u(x,t)
c     penetrating the outer neuron membrane.
c     Finally implemented coordinate dependent chemical potential amju(x)
c      
c           
c     Diffusion equation by chemical potential dependending on coordinate
c     by A. I.Livshits, Phys. Lett. 380 (2016) 1891-1894.
c     On the left hand side is always von Neumann BC applied
c     On the right hand side can be either von Neumann BC (for debugging) or Dirichlet BC
c     J.Mavri, D.Pregeljc, Ljubljana,London, September 2020
c     Cylindrical coordinates revisited
c     Origin is at the center of neuron
c     x() is radius
c
c     du(x,t)/dt = D* { d2u(x,t)/dx2 +(1/x)*du(x,t)/dx +
c                       (1/akT)*(du(x,t)/dx)*(dmju(x)/dx) +
c                       (1/akT)*u(x,t)*d2mju(x)/dx2  }   
c                   
c     Exp. diffusion coefficients F. Brouneus et al., Int. J. Pharmaceutics, 218 (2001) 57-62
c     D = 74.9d0               lidocaine  A^2/ns
c     D = 67.1d0               bupivacain A^2/ns
c     time unit is 1ns=10^-9s      
c     
c     J.Mavri, Sv. Lenart, January 2022
c     First and second derivatives can be calculated by FFT
c     V. Smrkolj, D. Pregeljc and J.Mavri, Ljubljana, March-May, 2022      
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit None
      Include 'sizes.diff_cyl'

      Double Precision x(MaxBin),U(MaxBin),U_New(MaxBin),amju(MaxBin),
     &     damju(MaxBin), d2amju(MaxBin),f(MaxBin),Pop(MaxBin),h1,h0,g1
      Double Precision a,b,dx,t,dt,alength, Pop_min, dU(MaxBin),
     & d2U(MaxBIn)
      Integer i,l,j,iend
      Integer(kind=8) :: Nstlim
      Double Precision Zero,akT,D,very_small,total,T1,T2,T3,T4
      Double Precision Pii, One, Third, a_membrane
      Double Precision a_membrane_start,a_membrane_end,pmf1,pmf2,pmf3
      Double Precision c0
      Integer N_wraps
      Double Precision water_slab, a_membrane_slab,r1,r2,r3
      Integer NBin,Ndruck, Nsmooth
      Logical Debug, Dirichlet, Fourier, Boltzmann, Lego_c, Lego_pmf
      Logical Axo_full
      Logical Ranvier, Schwann
      Debug=.true.
      Read(5,*) a,b,dx     ! left boundary, right boundary, dx
      Read(5,*) Nstlim,dt  ! number of steps, applied time step for integration
      Read(5,*) Dirichlet  ! if .true. then we will have Dirchlet (=open) BC on the right
      Read(5,*) Fourier    ! if .true. then first and second derivatives are calculated by FFT
      Read(5,*) Boltzmann  ! if .true. then initial populations are according to Boltzmann distribution
      Read(5,*) Lego_c     ! if .true. then initial concentration profile is read from file  'initial_c'
      Read(5,*) Axo_full,c0 ! if .true. then the entire axoplasm region is filled with concentration c0
      Read(5,*) g1         ! Dirichlet BC, steady concentration on the right
      Read(5,*) Ndruck     ! how often intermediate results are printed
      Read(5,*) Nsmooth    ! how many times smoothing of the potential of mean force is perfomed
      Read(5,*) D          ! diffusion coefficient A^2/ns; 74.9d0 for lidocaine and 67.1d0  for bupivacaine
      Read(5,*) a_membrane_start,a_membrane_end,pmf1,pmf2,pmf3  !radius of the membrane start, membrane end, and the corresponding pmfs
      Read(5,*) Schwann    ! if .true. then we will generate pmf for membrane at the Schwann cell
      If (Schwann) Then
c        in addition read  number of wraps, thickness of water_slab and a_membrane_slab
         Read(5,*) N_wraps, water_slab, a_membrane_slab
      End If

      Write(6,*)'On the left side we will have von Neumann BC'
      If (Dirichlet) Then
         Write(6,*) 'On the right hand side we will have Dirichlet BC'
      Else
         Write(6,*)'On the right side we will have von Neumann BC'
      End If

      If (Fourier) Then
         Write(6,*) 'First and second derivatives by FFT'
      End If

      If (Boltzmann) Then
         Write(6,*) 'Inital populations Boltzmann distributed'
      Else If (Lego_c) Then
         Write(6,*) 'Inital populations imported from initial_c'
      Else
         Write(6,*) 'Inital populations generated automatically'
      End If

      If (Schwann) Then
         Write(6,*) 'We will preform calculations
     & for membrane at Schwann cell'
      End If

      If (Axo_full) Then
         Write(6,*) 'Entire axoplasm region
     & is filled with concentraion c0'
      End If
      
      Write(6,*) 'a,b,dx ',a,b,dx
      Write(6,*) 'Nstlim,dt ', Nstlim,dt
      Nbin= int((b-a)/dx) + 1
      alength = b - a
      Write(6,*) 'Nbin', Nbin
      Write(6,*) 'alength', alength
      Write(6,*) 'Ndruck', Ndruck
      Write(6,*) 'Nsmooth', Nsmooth
      Write(6,*) 'D', D
      Write(6,*) 'a_membrane_start,a_membrane_end',
     & a_membrane_start,a_membrane_end
      Write(6,*) 'pmf1 - axoplasm ,pmf2 - membrane, pmf3 - ECF',
     & pmf1,pmf2,pmf3
      If (Schwann) Then
         Write(6,*) 'N_wraps, water_slab, a_membrane_slab',
     & N_wraps, water_slab, a_membrane_slab
      End If

c     Security
      If (NBin.gt.MaxBin) Then
         Write(6,*) 'NBin>MaxBin'
         Stop
      End If

      Zero = 0.0d0
c      very_small = 1.0d-10
      very_small = 0.0d0
      One=1.0d0
      Pii=4.0d0*atan(One)
      Third=One/3.0d0

      h0=Zero  ! von Neumann BC, at left impermeable wall
      h1=Zero  ! von Neumann BC, at right impermeable wall
      akT = 0.593d0             !k_B*T in kcal/mol at room temperature
c      D = 74.9d0               !Exp. lidocaine diffusion coefficient A^2/ns
c      D = 67.1d0               !Exp. bupivacaine diffusion coefficient A^2/ns

  
c     Sources and sinks along the coordinate
      Do i=1,NBin
         f(i) = Zero
      End Do

c     Set up the grid
      Do i=1,NBin
           x(i) = a + (i-1)*dx
      End Do     

      If  (Lego_pmf) Then
c        Read the potential of mean force to be applied from a file
         Open(Unit=10, File='PMF',Status='Old')
         Do i=1,NBin
            Read(10,*) amju(i)
         End Do
         Close(Unit=10)
      Else
c     Generate it manually    
c     Potential of mean force; data for lidocaine at physiological pH
c     set pmf of the extracellular fluid to zero       
      
c        axoplasm      
         Do i=1,int(a_membrane_start/dx)
            amju(i) = pmf1
         End Do
      
c        membrane
         Do i=int(a_membrane_start/dx)+1, int(a_membrane_end/dx)
            amju(i) = pmf2
         End Do

c        extracellular liquid
         Do i=int(a_membrane_end/dx)+1, int(b/dx)+1
            amju(i) = pmf3
         End Do
         If (Schwann) Then
c         Generate in addition pmf for Schwann cell layers
c         Loop over 40 Schwann cell layers            
            Do i=1, N_wraps
               r1 = a_membrane_end +(i-1)*(water_slab + a_membrane_slab)
               r2 = r1 + water_slab
               r3 = a_membrane_end + i*(water_slab + a_membrane_slab)
               Do j=int(r1/dx)+1, int(r2/dx)
                  amju(j) = pmf3 !extracellular fluid
               End Do
               Do j=int(r2/dx)+1, int(r3/dx)
                  amju(j) = pmf2 !the same pmf as in the membrane pmf2
               End Do
            End Do
         End If
                 
c        Smooth potential of mean force in order to make integration stable
         Do l=1, Nsmooth
            Do i=2,Nbin-1
               amju(i) = Third*amju(i)+
     &                   Third*amju(i-1) + Third*amju(i+1)
            End Do
         End Do
      End If
c
      Open (Unit=8, File='GENERATED_PMF',Status='Unknown')
      Do i=1,NBin
         write(8,*) x(i),amju(i)
      End Do
      Close(Unit=8)
      Write(6,*) 'Applied pmf is in file GENERATED_PMF'

      If (Lego_c) Then
c        Read the initial concentration profile
         Open(Unit=10, File='INITIAL_C', Status='Old')
         Do i=1,NBin
            Read(10,*) U(i)
         End Do
         Close(Unit=10)
      Else If (Boltzmann) Then
         Write(6,*) 'Equilibrium Pop(i) according to amju(i) calculated'
         Do i=1,Nbin
            Pop(i) = exp(-amju(i)/akT)
         End Do

         Pop_min = Pop(1)
         Do i=1,Nbin
            If (Pop(i).lt. Pop_min) Then
               Pop_min = Pop(i)
            End If
         End Do

         Write(6,*) 'Scaled eq.populations are calculated'
      
         Write(6,*) 'Pop_min', Pop_min
         Do i=1,Nbin
            Pop(i)=Pop(i)/Pop_min
         End Do

         Write(6,*)'u(x,t=0) corresponds to equilibrium'
         Do i=1,Nbin
            U(i) = Pop(i)
         End Do
      Else
c     We will generate intial concentration profile manually
         Do i=1,Nbin
            U(i)=very_small
         End Do
         If (Axo_full) Then
c            we will fill only the axoplasm  with concentration c0 in order to study transfer to the membrane/Schwann cell
            Do i=int(a_membrane_start/dx), int(a_membrane_end/dx)
               U(i) = c0
            End Do
         Else
            If (Schwann) Then
               iend = int((a_membrane_end +
     &              N_wraps*(water_slab + a_membrane_slab) ) /dx)
            Else
               iend = int((a_membrane_end)/dx)
            End If
            Do i=iend, Nbin
               U(i) = c0
            End Do
         End If
c        Smooth initial concentration profile in order to make integration stable
         Do l=1, Nsmooth
            Do i=2,Nbin-1
               U(i) = Third*U(i)+
     &                   Third*U(i-1) + Third*U(i+1)
            End Do
         End Do
      End If
      
      Write(6,*) 'Inital profile is stored in initial'
      Open(Unit=10, File='INITIAL', Status='Unknown')
      Do i=1,Nbin
         Write(10,*) x(i), U(i)
      End Do
      Close(Unit=10)
      
c     First derivative of amju(i) with respect to the coordinate
      If (Fourier) Then
            Call FD1D(Nbin,amju,damju,alength)
      Else 
         Do i=2, Nbin-1
            damju(i) = (amju(i+1)-amju(i-1))/(dx)
         End Do
         damju(1) = (amju(2)-amju(1))/dx
         damju(Nbin) = (amju(Nbin)-amju(Nbin-1))/dx
      End If

      Open(Unit=16,File='damju',Status='Unknown')
      Do i=1,Nbin
         Write(16,*) x(i), damju(i)
      End Do
      Close(Unit=16)
      
c     Second Derivative of amju(i) with respect to the coordinate
      If (Fourier) Then
         Call SD1D(Nbin,amju,d2amju,alength)
      Else
         Do i=2, Nbin-1
            d2amju(i) = (amju(i+1)-2.0d0*amju(i)+amju(i-1))/(dx*dx)
         End Do
         d2amju(1) = d2amju(2)
         d2amju(Nbin) = d2amju(Nbin-1)
      End If

      Open(Unit=17,File='d2amju',Status='Unknown')
      Do i=1,Nbin
         Write(17,*) x(i), d2amju(i)
      End Do
      Close(Unit=17)

      Write(6,*) 'Diffusion equation with included chemical potential'
      Write(6,*) 'J. Mavri,D. Pregeljc, L. Smrkolj LJ 2020-2022'
      write(6,2010) dt,nstlim
 2010 format(/' dt ',f10.5' fs.  nstlim',i10/)
      write(6,2011) dt*nstlim
 2011 format('Total simulation time',e15.5,'ns')


      Open(Unit=13,File='c_t',Status='Unknown')
      Open(Unit=15,File='MEMBRANE',Status='Unknown')

c     Propagation loop ******************************
      
      Do l=1,Nstlim
         t = (l-1)*dt
         Do i=2,Nbin-1
            If (Fourier) Then
c              FFT first and second derivatives of U(x,t)              
               Call FD1D(Nbin,U,dU,alength)
               Call SD1D(Nbin,U,d2U,alength)
               T1 = d2U(i)
               T2 = (One/x(i))*dU(i)
               T3 = (One/akT)*dU(i)*damju(i)
               T4 = (One/akT)*U(i)*d2amju(i)
            Else
               T1 = (U(i+1)-2.0*U(i)+U(i-1))/(dx*dx)
c               T2 = (One/x(i))*((U(i+1)-U(i-1))/(dx))
               T2 = 0
               T3 = (One/akT)*(U(i+1)-U(i-1))/(dx)*damju(i)
               T4 = (One/akT)*U(i)*d2amju(i)
             End If
            U_new(i) = U(i) + D*dt*(T1 + T2 + T3 + T4)
            If (U_New(i).lt.Zero) Then
               U_New(i) = very_small
c               Write(6,*) 'CONCENTRATION NEGATIVE! Timestep:',l, 'x:',i
            End If
         End Do


c     Here we will apply von Neumann BC: right boundary is impermeable,
c     derivative of the concentration is fixed to zero
        
         If (Dirichlet) Then
            U_New(Nbin) = g1
         Else
         U_New(Nbin) = U(Nbin) + D*(dt/(dx*dx)*(U(Nbin-1)- U(Nbin)
     &           +h1*dx) + f(Nbin))      
         End If
c    left boundary is by definition impermeable since we have cylindrical coordinates
c    this is center of the neuron         

        U_New(1) = U(1) + D*(dt/(dx*dx)*(U(2)-U(1)+h0*dx)+f(1))

c        Copy new versions of U to old        
         Do i=1,Nbin
            U(i) = U_New(i)
         End Do

c     Write the concentration profile
         If (Mod(l,Ndruck).eq.0) Then
              Do i=1,NBin
                   Write(13,'(2e20.10)') x(i),U(i) 
              End Do
              Write(13,*) '  '  
              Write(13,*) '  '  
         End If

c        population in the membrane or in the Schwann cell
         a_membrane = Zero
         If (Schwann) Then
            iend = int((a_membrane_end +
     &           N_wraps*(water_slab + a_membrane_slab) ) /dx)
         Else
            iend = int((a_membrane_end)/dx)
         End If
         
         Do i=int(a_membrane_start/dx),iend
              a_membrane = a_membrane + Pii*U(i)*x(i)*x(i)*dx
         End Do
         Write(15,'(2e15.5)') t, a_membrane

      End Do        
c     End of propagation loop *******************************
      Close(Unit=13)
c      Close(Unit=14)
      Close(Unit=15)

c     Write the final concentrations to a file final
      Open(Unit=12,File='FINAL',Status='Unknown')
      Do i=1,NBin
         Write(12,'(2e20.10)') x(i),U(i)
      End Do
      Close(Unit=12)

      Write(6,*) 'Real time simulation duration ', Nstlim*dt, 'ns'
      Write(6,*) 'What corresponds to ', Nstlim*dt*1.0d-9, 's'

      Write(6,*) 'Initial profile is written in initial'
      Write(6,*) 'Final profile is written in final'
      Write(6,*) 'membrane population(t) in ns is in file MEMBRANE '

      Write(6,*) 'Job done.'


      End


