!!! Lucas McCullum !!!
!!! ENME 489       !!!
!!! Dr. Meilin Yu  !!!
!!! Project 2U     !!!

program plotting

    implicit none

    ! Define all of the parameters
    real, parameter :: rho = 1.0                    ! Thermal diffusivity
    real, parameter :: c = 1.0                      ! Space start point
    real, parameter :: minS = -8.0                  ! Space start point
    real, parameter :: maxS = 8.0                   ! Space end point
    real, parameter :: minT = 0.0                   ! Time start point
    real, parameter :: maxT = 2.0                   ! Time end point
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

        ! Two-stage TVD Runge-Kutta
        print*,"Two-stage TVD Runge-Kutta"
        do j = 1,smax
            Res1 = (alpha/(h**2))*(Res_mat(j+1)-2*Res_mat(j)+Res_mat(j-1))+S(scheme,x_mat(j))
            Res = (alpha/(h**2))*(S_mat(j+1)-2*S_mat(j)+S_mat(j-1))+S(scheme,x_mat(j))
            u_1 = S_mat(j)+k*Res
            FD_mat(j) = (0.5)*(S_mat(j)+u_1)!+k*Res1)
        end do

        FD_mat(0) = TminS
        if (scheme == "13") then
            FD_mat(smax+1) = FD_mat(smax)
        else
            FD_mat(smax+1) = TmaxS(scheme,t_val)
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
    open(unit=48, file='Project2_U.dat')
    do i = 0,smax+1
        current_val = x_mat(i)
        print*,"(",current_val,",",F(scheme,current_val),",",FD_mat(i),")"
        write(48,*) current_val,F(scheme,current_val),FD_mat500(i),FD_mat1500(i),FD_mat(i)
    end do
    close(48)
	! I am plotting with Gnuplot so this calls that program for plotting... must have the file named that for it to work
    call system('gnuplot -p Project2_U.plt')

	! Sanity check to see if everything is working right through resolution and stability
    print*,"Time steps: ",tmax,"Space steps: ",smax
    print*,"Stability is: ",r

    ! Import all of the required outside functions...

contains

    ! Determine the boundary conditions for the fluid pressure
    function P_Bound(edge) result(v)

        implicit none

        character (len = 3), intent(in) :: edge   ! Input = define which scheme to use
        real :: v                                  ! Output = the temperature at the max x-bound

        if (edge == "min") then
            v = 0
        else if (edge == "max") then
            v = 0
        end if

    end function P_Bound

    ! Determine the boundary conditions for the fluid velocity
    function U_Bound(edge) result(v)

        implicit none

        character (len = 3), intent(in) :: edge   ! Input = define which scheme to use
        real :: v                                  ! Output = the temperature at the max x-bound

        if (edge == "min") then
            v = 0
        else if (edge == "max") then
            v = 0
        end if

    end function U_Bound

    ! Initial fluid pressure at t = 0
    function P(x) result(y)

        implicit none

        real, intent(in) :: x                       ! Input = position along space mesh
        real :: y                                   ! Output = value of initial function at this position

        y = exp(-(x**2))

    end function P

    ! Initial fluid velocity at t = 0
    function U(x) result(y)

        implicit none

        real, intent(in) :: x                       ! Input = position along space mesh
        real :: y                                   ! Output = value of initial function at this position

        y = (0.5)*exp((-2)*(x**2))

    end function U

end program plotting
