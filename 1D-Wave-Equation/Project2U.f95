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
    real, parameter :: intervals = 800              ! Number of intervals (mesh size)
    real, parameter :: h = (maxS-minS)/intervals    ! Space step
    real, parameter :: k = (maxT-minT)/intervals   	! Time step
    integer, parameter :: tmax = (maxT-minT)/k      ! Number of time steps
    integer, parameter :: smax = ((maxS-minS)/h)-1  ! Number of space steps (nodes)

    ! Define all of the input and output variables
    real :: r                       ! Find the stability of the scheme... also known as the CFL number
    real :: t                       ! The value of the current time step
    real :: ResP                    ! TVD R-K variable (pressure)
    real :: ResU                    ! TVD R-K variable (velocity)
    real :: ResP1                   ! TVD R-K variable (pressure)
    real :: ResU1                   ! TVD R-K variable (velocity)
    real :: ResP_mat(0:smax+1)      ! TVD R-K variable (pressure)
    real :: ResU_mat(0:smax+1)      ! TVD R-K variable (velocity)
    real :: uP_1                    ! TVD R-K variable (pressure)
    real :: uU_1                    ! TVD R-K variable (velocity)
    real :: P_mat(0:smax+1)         ! Array of final pressure values
    real :: P_mat500(0:smax+1)      ! Array of pressure values at t = 0.5
    real :: P_mat1500(0:smax+1)     ! Array of pressure values at t = 1.0
    real :: U_mat(0:smax+1)         ! Array of final velocity values
    real :: U_mat500(0:smax+1)      ! Array of velocity values at t = 0.5
    real :: U_mat1500(0:smax+1)     ! Array of velocity values at t = 1.0
    real :: SP_mat(0:smax+1)        ! Array of temporary pressure spatial values
    real :: SU_mat(0:smax+1)        ! Array of temporary velocity spatial values
    real :: x_mat(0:smax+1)         ! Array of x values
    integer :: i                    ! Variable to be looped in do loop
    integer :: j                    ! Variable to be looped in do loop
    integer :: z                    ! Variable to be looped in do loop
    real :: current_val             ! Temp variable to determine spatial value at current spatial step
    real :: t_val                   ! Temp variable to determine temporal value at current spatial step

    !!! Start the main code !!!

    ! Create an array for the initial condition
    x_mat(0) = minS
    SP_mat(0) = P_Bound("min")!(0.5)*(TminS+F(scheme,minS))
    SU_mat(0) = U_Bound("min")!(0.5)*(TminS+F(scheme,minS))

    do i = 1,smax
        current_val = x_mat(i-1)+h
        x_mat(i) = current_val
        !print*,"Current: ",F(scheme,current_val)
        SP_mat(i) = P0(current_val)
        SU_mat(i) = U0(current_val)
    end do

    x_mat(smax+1) = maxS
    SP_mat(smax+1) = P_Bound("max")!0.5)*(TmaxS(scheme,0.0)+F(scheme,maxS))
    SU_mat(smax+1) = U_Bound("max")!0.5)*(TmaxS(scheme,0.0)+F(scheme,maxS))

    do i = 0,smax+1
        ResP_mat(i) = SP_mat(i)
        ResU_mat(i) = SU_mat(i)
        !print*,Res_mat(i)
    end do

    ! Move through time and update spatial components
    do i = 1,tmax

        t_val = minT+((i-1)*k)

        ! Two-stage TVD Runge-Kutta
        print*,"Two-stage TVD Runge-Kutta"
        do j = 1,smax
            !ResP1 = -rho*(c**2)*( (1/h)*(ResP_mat(j)-ResP_mat(j-1)) )
            !ResP = -rho*(c**2)*( (1/h)*(SP_mat(j)-SP_mat(j-1)) )!Pressure (dp/dt):-rho*(c**2)*(du/dx) BDM
            !uP_1 = SP_mat(j)+k*ResP
            !P_mat(j) = (0.5)*(SP_mat(j)+uP_1+k*ResP1)

            !ResU1 = -(1/rho)*( (1/h)*(ResU_mat(j)-ResU_mat(j-1)) )
            !ResU = -(1/rho)*( (1/h)*(SU_mat(j)-SU_mat(j-1)) )!Velocity (du/dt):-(1/rho)*(dp/dx) BDM
            !uU_1 = SU_mat(j)+k*ResU
            !U_mat(j) = (0.5)*(SU_mat(j)+uU_1+k*ResU1)

            !-rho*(c**2)*( (1/h)*(SU_mat(j)-SU_mat(j-1)) )!Pressure (dp/dt):-rho*(c**2)*(du/dx) BDM
            !-(1/rho)*( (1/h)*(SP_mat(j)-SP_mat(j-1)) )!Velocity (du/dt):-(1/rho)*(dp/dx) BDM

            ! First-Order
            if ((j == 1) .or. (j == smax)) then
                ResP = -1*( (0.5/h)*( (SP_mat(j)+SU_mat(j))-(SP_mat(j-1)+SU_mat(j-1)) ) &
                    + (0.5/h)*( (-SP_mat(j+1)+SU_mat(j+1))-(-SP_mat(j)+SU_mat(j)) ) )
                ResU = -1*( (0.5/h)*( (SP_mat(j)+SU_mat(j))-(SP_mat(j-1)+SU_mat(j-1)) ) &
                    + (0.5/h)*( (SP_mat(j+1)-SU_mat(j+1))-(SP_mat(j)-SU_mat(j)) ) )

                !ResP1 = ResP
                !ResU1 = ResU
            else
                ! Second-Order
                ResP = -1*( (0.5/h)*( 3*(SP_mat(j)+SU_mat(j))-4*(SP_mat(j-1)+SU_mat(j-1))+1*(SP_mat(j-2)+SU_mat(j-2)) ) &
                    + (0.5/h)*( -1*(-SP_mat(j+2)+SU_mat(j+2))+4*(-SP_mat(j+1)+SU_mat(j+1))-3*(-SP_mat(j)+SU_mat(j)) ) )
                ResU = -1*( (0.5/h)*( 3*(SP_mat(j)+SU_mat(j))-4*(SP_mat(j-1)+SU_mat(j-1))+1*(SP_mat(j-2)+SU_mat(j-2)) ) &
                    + (0.5/h)*( -1*(SP_mat(j+2)-SU_mat(j+2))+4*(SP_mat(j+1)-SU_mat(j+1))-3*(SP_mat(j)-SU_mat(j)) ) )

                !ResP1 = ResP
                !ResU1 = ResU
            end if

            uP_1 = SP_mat(j)+k*ResP
            uU_1 = SU_mat(j)+k*ResU

            !-rho*(c**2)*( (1/h)*(SU_mat(j)-SU_mat(j-1)) )
            !-(1/rho)*( (1/h)*(SP_mat(j)-SP_mat(j-1)) )

            ResP1 = ResP
            ResU1 = ResU

            P_mat(j) = (0.5)*(SP_mat(j)+uP_1+k*ResP1)
            U_mat(j) = (0.5)*(SU_mat(j)+uU_1+k*ResU1)
        end do

        P_mat(0) = P_Bound("min")
        P_mat(smax+1) = P_Bound("max")

        U_mat(0) = U_Bound("max")
        U_mat(smax+1) = U_Bound("max")

        do z = 0,smax+1
            SP_mat(z) = P_mat(z)
            SU_mat(z) = U_mat(z)
            if (t_val == 0.5) then
                P_mat500(z) = P_mat(z)
                U_mat500(z) = U_mat(z)
            end if
            if (t_val == 1.0) then
                P_mat1500(z) = P_mat(z)
                U_mat1500(z) = U_mat(z)
            end if
        end do

    end do

	! Export the corresponding data to a '.dat' file for future plotting
    open(unit=48, file='Project2_Ut_test.dat')
    do i = 0,smax+1
        current_val = x_mat(i)
        print*,"(",current_val,",",P0(current_val),",",P_mat(i),",",U0(current_val),",",U_mat(i),")"
        write(48,*) current_val,P0(current_val),P_mat500(i),P_mat1500(i),P_mat(i),U0(current_val),U_mat500(i),U_mat1500(i),U_mat(i)
    end do
    close(48)

	! I am plotting with Gnuplot so this calls that program for plotting... must have the file named that for it to work
    call system('gnuplot -p Project2_U.plt')

	! Sanity check to see if everything is working right through resolution and stability
    print*,"Time steps: ",tmax,"Space steps: ",smax
    !print*,"Stability is: ",r

    ! Import all of the required outside functions...

contains

    ! Determine the boundary conditions for the fluid pressure
    function P_Bound(edge) result(v)

        implicit none

        character (len = 3), intent(in) :: edge     ! Input = define which scheme to use
        real :: v                                   ! Output = the temperature at the max x-bound

        if (edge == "min") then
            v = 0
        else if (edge == "max") then
            v = 0
        end if

    end function P_Bound

    ! Determine the boundary conditions for the fluid velocity
    function U_Bound(edge) result(v)

        implicit none

        character (len = 3), intent(in) :: edge     ! Input = define which scheme to use
        real :: v                                   ! Output = the temperature at the max x-bound

        if (edge == "min") then
            v = 0
        else if (edge == "max") then
            v = 0
        end if

    end function U_Bound

    ! Initial fluid pressure at t = 0
    function P0(x) result(y)

        implicit none

        real, intent(in) :: x                       ! Input = position along space mesh
        real :: y                                   ! Output = value of initial function at this position

        y = exp(-(x**2))

    end function P0

    ! Initial fluid velocity at t = 0
    function U0(x) result(y)

        implicit none

        real, intent(in) :: x                       ! Input = position along space mesh
        real :: y                                   ! Output = value of initial function at this position

        y = (0.5)*exp((-2)*(x**2))

    end function U0

end program plotting
