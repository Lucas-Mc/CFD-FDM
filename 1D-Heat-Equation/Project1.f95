!!! Lucas McCullum !!!
!!! ENME 498       !!!
!!! Dr. Meilin Yu  !!!

program plotting

    implicit none

    ! Define all of the parameters
    character (len = 3), parameter :: scheme = "22"! Define the scheme to be used
    ! 111-Forward Euler (done)
    ! 112-Two_stage TVD RK (done)
    ! 113-Backward Euler (done)
    ! 114-Crank_Nicolson (done)
    ! 12-Two_stage TVD R-K (done)
    ! 13-Crank_Nicolson (done)
    ! 21-Two_stage TVD R-K (done)
    ! 22-Crank_Nicolson (done)
    real, parameter :: alpha = 1.0                  ! Thermal diffusivity
    real, parameter :: minS = -1.0                  ! Space start point
    real, parameter :: maxS = 1.0                   ! Space end point
    real, parameter :: minT = 0.0                   ! Time start point
    real, parameter :: maxT = 1.5                   ! Time end point
    real, parameter :: h = 0.025                    ! Space step
    real, parameter :: k = 0.00025                 	! Time step
    real, parameter :: TminS = 5.0                  ! T(x at left bound)
    integer, parameter :: tmax = (maxT-minT)/k      ! Number of time steps
    integer, parameter :: smax = ((maxS-minS)/h)-1  ! Number of space steps (nodes)

    ! Define all of the input and output variables
    real :: r                       ! Find the stability of the scheme... also known as the CFL number
    real :: t                       ! The value of the current time step
    real :: Res                     ! TVD R-K variable
    real :: Res1                    ! TVD R-K variable
    real :: Res_mat(0:smax+1)       ! TVD R-K variable
    real :: u_1                     ! TVD R-K variable
    real :: FD_mat(0:smax+1)        ! Array of final temperature values
    real :: FD_mat500(0:smax+1)     ! Array of temperature values at t = 0.5
    real :: FD_mat1500(0:smax+1)    ! Array of temperature values at t = 1.0
    real :: S_mat(0:smax+1)         ! Array of temporary spatial values
    real :: x_mat(0:smax+1)         ! Array of x values
    real, dimension (smax,3) :: At  ! Array of temp values for Backward Euler method (Al=B)
    real, dimension (smax,1) :: Bt  ! Array of temp values for Backward Euler method (Al=B)
    real :: a(1:smax)               ! Array of temp Crank-Nicolson values
    real :: b(1:smax)               ! Array of temp Crank-Nicolson values
    real :: c(1:smax)               ! Array of temp Crank-Nicolson values
    real :: d(1:smax)               ! Array of temp Crank-Nicolson values
    integer :: i                    ! Variable to be looped in do loop
    integer :: j                    ! Variable to be looped in do loop
    integer :: z                    ! Variable to be looped in do loop
    real :: current_val             ! Temp variable to determine spatial value at current spatial step
    real :: t_val                   ! Temp variable to determine temporal value at current spatial step

    !!! Start the main code !!!

    ! Determine the stability where r < 1/2 is stable (also lambda)
    r = (alpha*k)/(h**2)

    ! Create an array for the initial condition
    x_mat(0) = minS
    S_mat(0) = (0.5)*(TminS+F(scheme,minS))

    do i = 1,smax
        current_val = x_mat(i-1)+h
        x_mat(i) = current_val
        !print*,"Current: ",F(scheme,current_val)
        S_mat(i) = F(scheme,current_val)
    end do
    if (scheme == "13") then
        S_mat(smax+1) = (0.5)*(S_mat(smax)+F(scheme,maxS))
    else
        S_mat(smax+1) = (0.5)*(TmaxS(scheme,0.0)+F(scheme,maxS))
    end if
    x_mat(smax+1) = maxS
    do i = 0,smax+1
        Res_mat(i) = S_mat(i)
        !print*,Res_mat(i)
    end do


    ! Move through time and update spatial components
    do i = 1,tmax

        t_val = minT+((i-1)*k)

        if (scheme == "113") then

            print*,"Backward Euler"
            ! Backward Euler
            do j = 1,smax
                At(j,2) = -1*(2*r+1)
                if (j == 1) then
                    Bt(1,1) = -S_mat(j)-r*TminS
                else
                    if (j /= smax) then
                        Bt(j,1) = -S_mat(j)
                    end if
                end if
                if (j == smax) then
                    Bt(smax,1) = -S_mat(j)-r*TmaxS(scheme,t_val)
                    At(j,1) = r
                else
                    if (j /= 1) then
                        At(j,1) = r
                    end if
                    At(j,3) = r
                end if
            end do

            do j = 2,smax
                At(j,2) = At(j,2)-At(j-1,3)*At(j,1)/At(j-1,2)
                Bt(j,1) = Bt(j,1)-Bt(j-1,1)*At(j,1)/At(j-1,2)
            end do

            FD_mat(0) = TminS
            FD_mat(smax) = Bt(smax,1)/At(smax,2)
            do j = smax-1,1,-1
                FD_mat(j) = (Bt(j,1)-At(j,3)*FD_mat(j+1))/At(j,2)
            end do
            FD_mat(smax+1) = TmaxS(scheme,t_val)

        else if (scheme == "114" .or. scheme == "13" .or. scheme == "22") then

            print*,"Crank-Nicolson"
            ! Crank-Nicolson
            do j = 1,smax
                a(j) = -r/2
                b(j) = 1+r
                c(j) = -r/2
                d(j) = (r/2)*S_mat(j-1)+(1-r)*S_mat(j)+(r/2)*S_mat(j+1)+(k/2)*(S(scheme,x_mat(j))+S(scheme,x_mat(j)))
            end do

            if (scheme == "13") then
                !d(1) = d(1) - 2*h*r*TminS/2
                d(smax) = d(smax) + 2*h*r*TmaxS(scheme,t_val)/2
                c(1) = -r
                a(smax) = -r
            else
                d(1) = d(1) + r*TminS/2
                d(smax) = d(smax) + r*TmaxS(scheme,t_val)/2
            end if

            call TRIDI(smax,a,b,c,d)

            do j = 1,smax
                FD_mat(j) = d(j)
            end do

            FD_mat(0) = TminS
            if (scheme == "13") then
                FD_mat(smax+1) = FD_mat(smax)
            else
                FD_mat(smax+1) = TmaxS(scheme,t_val)
            end if

        else

            do j = 1,smax
                if (scheme == "112" .or. scheme == "12" .or. scheme == "21") then
                    ! Two-stage TVD Runge-Kutta
                    print*,"Two-stage TVD Runge-Kutta"
                    Res1 = (alpha/(h**2))*(Res_mat(j+1)-2*Res_mat(j)+Res_mat(j-1))+S(scheme,x_mat(j))
                    Res = (alpha/(h**2))*(S_mat(j+1)-2*S_mat(j)+S_mat(j-1))+S(scheme,x_mat(j))
                    u_1 = S_mat(j)+k*Res
                    FD_mat(j) = (0.5)*(S_mat(j)+u_1)!+k*Res1)
                else
                    ! Forward Euler
                    print*,"Forward Euler"
                    FD_mat(j) = r*S_mat(j-1)+(1-(2*r))*S_mat(j)+r*S_mat(j+1)+k*S(scheme,x_mat(j))
                end if
            end do
			
            FD_mat(0) = TminS
            if (scheme == "13") then
                FD_mat(smax+1) = FD_mat(smax)
            else
                FD_mat(smax+1) = TmaxS(scheme,t_val)
            end if

        end if

        do z = 0,smax+1
            S_mat(z) = FD_mat(z)
            if (t_val == 0.5) then
                FD_mat500(z) = FD_mat(z)
            end if
            if (t_val == 1.0) then
                FD_mat1500(z) = FD_mat(z)
            end if
        end do

    end do

	! Export the corresponding data to a '.dat' file for future plotting
    open(unit=48, file='Project1_22.dat')
    do i = 0,smax+1
        current_val = x_mat(i)
        print*,"(",current_val,",",F(scheme,current_val),",",FD_mat(i),")"
        write(48,*) current_val,F(scheme,current_val),FD_mat500(i),FD_mat1500(i),FD_mat(i)
    end do
    close(48)
	! I am plotting with Gnuplot so this calls that program for plotting... must have the file named that for it to work
    call system('gnuplot -p Project1_FDM.plt')

	! Sanity check to see if everything is working right through resolution and stability
    print*,"Time steps: ",tmax,"Space steps: ",smax
    print*,"Stability is: ",r

    ! Import all of the required outside functions...

contains

    ! Determine the boundary conditions at the maxS
    function TmaxS(scheme,t) result(mt)

        character (len = 3), intent(in) :: scheme   ! Input = define which scheme to use
        real :: t                                   ! Temp = used to setup the optional parameter tt
        real :: mt                                  ! Output = the temperature at the max x-bound
        real(16), parameter :: PI = 4*atan(1.0_16)  ! Constant = pi ~= 3.14159...

        if (scheme == "12") then
            mt = 1+sin(PI*t)
        else
            mt = 1
        end if

    end function

    ! Initial function at t = 0
    function F(scheme,x) result(y)

        character (len = 3), intent(in) :: scheme   ! Input = define which scheme to use
        real, intent(in) :: x                       ! Input = position along space mesh
        real :: y                                   ! Output = value of initial function at this position

        if (scheme /= "12") then
            y = -4*(x**2)-2*x+7
        else
            y = -2*x+3
        end if

    end function F

    ! Forcing / source function
    function S(scheme,x) result(v)

        character (len = 3), intent(in) :: scheme   ! Input = define which scheme to use
        real, intent(in) :: x                       ! Input = position along space mesh
        real :: v                                   ! Output = value of forcing function at this position

        if (scheme == "21" .or. scheme == "22") then
            if (x>=-0.5 .and. x<=0.5) then
                v = exp(-40*(x**2))+1
            else
                v = 0
            end if
        else
            v = 0
        end if

    end function S

	! Solve the tridiagonal matrix through this clever technique
    subroutine TRIDI(N,A,B,C,D)

        implicit none

        real, dimension (1:N) :: A,B,C,D    ! Input = Polynomial matrix coefficients
        integer, intent(in) :: N            ! Input = Size of the matrix
        integer :: I                        ! Temp = Loop variable
        real :: RATIO                       ! Temp = Temporary variable to store the ratio of values

        do I = 2,N
            RATIO = A(I)/B(I-1)
            B(I) = B(I) - RATIO*C(I-1)
            D(I) = D(I) - RATIO*D(I-1)
        end do

        D(N) = D(N)/B(N)
        do I = N-1,1,-1
            D(I) = (D(I) - C(I)*D(I+1))/B(I)
        end do

    return
	
    end

end program plotting
