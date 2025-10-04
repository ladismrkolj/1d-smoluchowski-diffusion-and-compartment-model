      Program Ode
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Solves a system of ordinary differential equations for a system
c     of coupled compartments. This will be used for phenomoenological
c     simulations of a synapse
c   
c     Units: SI
c
c     JM,DP Ljubljana, February 2018
c     SKETCH OF THE MODEL:
c                 ____________     ak(8)     ______________
c                |  MEMBRANE  | ----------> |   AXOPLASM   |
c                |     (3)    | <---------- |     (4)      |
c                |____________|    ak(7)    |______________|
c                 A  /                      A  /
c          ak(3) /  / ak(4)          ak(6) /  / ak(5)
c               /  /                      /  /
c  bkw(1)  ____/__V   ak(1)   ___________/  / 
c  -----> |  ECF  | -------> |  SCHWANN  | /
c  <----- |  (1)  | <------- |    (2)    |V
c  akw(1) |_______|   ak(2)  |___________|
c           |   A
c     ak(10)|   |ak(9)
c           V   |
c       ________________     
c      |     empty    |
c      |       (5)      |
c      |________________|
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit None
      Include 'sizes.ode'
      INTEGER(8) Npoi, Ndruck, i
      Integer j,nthree,n_ode,n_ak
      Double Precision dt, dydx(NMAX), zeit, one, zero
      Double Precision funk(NMAX), yout(NMAX), derivative(NMAX), 
     &     Alevel(NMAX),Avolume
      
      double precision ak,akw,bkw,tau1,blevel  ! coupling coefficiens for ODE, ak-between compartments, akw is washout coefficient
      Logical empty, washout, exponential
      common/coeff/ak(NMAX*2),akw(NMAX),bkw(NMAX),Avolume(NMAX),
     &             tau1,blevel
      common/logic/empty,exponential
      data zero /0.0d0/
      data one  /1.0d0/
      data nthree /3/
c      character(len=2) stri
           
      
      Write(6,*) 'Welcome in the program Ode'
      Write(6,*) 'Solving ODE for a system of compartments by rk4'
      Write(6,*) 'JM,DP Ljubljana, February 2018'
      Write(6,*) ' '

      Read(5,*) zeit
      Read(5,*) Npoi
      Read(5,*) dt
      Read(5,*) Ndruck
      Read(5,*) empty
      n_ode = 3
      n_ak = 6
      Read(5,*) washout, akw(1)
      If (.NOT. washout) Then
c     If we do not simulate washout from compartments, set all akw to 0, else read them from INP file
            Do i=1,n_ode
                  akw(i) = zero
            End Do
      End If
      Read(5,*) exponential, bkw(1), blevel, tau1
      If (.NOT. exponential) Then
            Do i=1,n_ode
                  bkw(i) = zero
            End Do
      End If
      Read(5,*) (Alevel(i),i=1,n_ode)    ! read initial concentrations
      Read(5,*) (Avolume(i),i=1,n_ode)   ! read volumes of the compartments
      Read(5,*) (Ak(i),i=1,n_ak)         ! read coupling coefficient between compartments
      
      Do i=1,n_ak
         If (Ak(i).ne.0) Then
            Ak(i) = 1/Ak(i)
         End If
      End Do

      If (akw(1).ne.0) Then
          akw(1)=1/akw(1)
      End If

      If (bkw(1).ne.0) Then
          bkw(1)=1/bkw(1)
      End If

      If (empty) Then
         bkw(1) = 0
         If (Alevel(1).eq.0) Then
            Alevel(1) = 1
         End If
         Alevel(2) = Alevel(1) * Ak(1) / Ak(2)
         Alevel(3) = Alevel(2) * Ak(6) / Ak(5)
      End If
      Write(6,'(a,i5)') 'n_ode', n_ode
      write(6,'(a,f10.5)') 'zeit', zeit
      write(6,'(a,i10)') 'npoi',npoi
      Write(6,'(a,1p,g15.2)') 'dt', dt
      Write(6,'(a,1p,g15.5)') 'simulation time', dt*npoi
      Write(6,*) 'empty', empty
      Write(6,'(a,1p,5g12.3)')
     &        '(Alevel(i),i=1,n_ode)', (Alevel(i),i=1,n_ode)
      write(6,'(a,1p,5e12.3)')
     &        '(Avolume(i),i=1,n_ode)', (Avolume(i),i=1,n_ode)
c      write(6,'(a,5f10.5)') '(Ak(i),i=1,n_ode)',(Ak(i),i=1,n_ak)

      Do i=1,n_ode
c           set funk(i) to initial concentrations provided by ALevel(i)
            funk(i)=ALevel(i)
c            write(stri,'(I5)') i
c            Open(Unit=(6+i),file=('A'//stri),Status='Unknown')
      End Do
      Open(Unit=7,file='A1',Status='Unknown')
      Open(Unit=8,file='A2',Status='Unknown')
      Open(Unit=9,file='A3',Status='Unknown')

c      Write(7,*) zeit, funk(1)
c      Write(8,*) zeit, funk(2)
c      Write(9,*) zeit, funk(3)
c      Write(10,*) zeit, funk(4)
c      If (empty) Then
c            Write(11,*) zeit, funk(5)
c      End If

      Do i=1,Npoi
c     Do runge kutta for each timestep
          zeit=(i-1)*dt
          Call Derivs(zeit,funk,derivative)
          Call RK4(funk,derivative,n_ode,zeit,dt,yout,Derivs)
c         Bellow we write the concentrations to file
          If (Mod(i,Ndruck).eq.0) Then
            Write(7,'(1p,e17.7,e15.5)') zeit, yout(1)
            Write(8,'(1p,e17.7,e15.5)') zeit, yout(2)
            Write(9,'(1p,e17.7,e15.5)') zeit, yout(3)
          End If
          Do j=1,n_ode
c            write((6+j),*) zeit, yout(j)
            funk(j)=yout(j)
          End Do
      End Do
c      Do i=1,n_ode
c          Close(Unit=(6+i))
c      End Do

      Close(Unit=7)
      Close(Unit=8)
      Close(Unit=9)


      Write(6,*) 'Time course of the levels is stored in A files.'
     
      Write(6,*) 'Bye.'
      End



      Subroutine Derivs(zeit,funk,derivative)
      implicit none
      include 'sizes.ode'
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Provides derivatives for integration using RK4
c     JM, DP, February 2018 
c
c     funk = concentration
c     dc(t)/dt = (dn(t)/dt) / V = SUM(-ak*funk) / Avolume
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision zeit
      double precision funk(NMAX)
      double precision derivative(NMAX)

      double precision ak,akw,bkw,Avolume,tau1,blevel
      logical empty, exponential
      common/coeff/ak(NMAX*2),akw(NMAX),bkw(NMAX),Avolume(NMAX),
     &             tau1,blevel
      common/logic/empty,exponential
          derivative(1) = (-akw(1)*funk(1))+
     &         (bkw(1) * blevel * EXP(-zeit/tau1))*(blevel-funk(1))+
     &                   ((-ak(1)*funk(1)+ak(2)*funk(2)) * Avolume(2)+
     &        (-ak(3)*funk(1)+ak(4)*funk(3)) * Avolume(3)) / Avolume(1)
          derivative(2) = ((ak(1)*funk(1)-ak(2)*funk(2))+
     &                    (-ak(6)*funk(2)+ak(5)*funk(3)))
          derivative(3) = ((ak(6)*funk(2)-ak(5)*funk(3)) * Avolume(2)+
     &        (-ak(4)*funk(3)+ak(3)*funk(1)) * Avolume(3)) / Avolume(3)

      return
      End







