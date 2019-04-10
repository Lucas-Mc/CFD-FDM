!!! Lucas McCullum !!!
!!! ENME 489       !!!
!!! Dr. Meilin Yu  !!!
!!! Project 2G     !!!

program plotting

    implicit none

    ! Define all of the parameters
    real, parameter :: minS = 0.0                   ! Space start point
    real, parameter :: maxS = 3.0                   ! Space end point
    real, parameter :: minT = 0.0                   ! Time start point
    real, parameter :: maxT = 2.0                   ! Time end point
    real, parameter :: gam = 1.4                  ! Fluid specific weight

    real, parameter :: intervals = 800              ! Number of intervals (mesh size)
    real, parameter :: h = (maxS-minS)/intervals    ! Space step
    real, parameter :: k = (maxT-minT)/intervals   	! Time step

    real, parameter :: u_out = 0.08                 ! Fluid exit velocity
    real, parameter :: rho_out = 1                  ! Fluid exit density
    real, parameter :: p_out = 1/gam                ! Fluid exit pressure
    real, parameter :: u_in = 0.08                  ! Fluid inlet velocity
    real, parameter :: rho_in = 1                   ! Fluid inlet density
    real, parameter :: p_in = 1/gam                 ! Fluid inlet pressure

    real,parameter :: al = 0.07                     ! Lower bound for Bisection Method
    real,parameter :: bh = 0.99                     ! Upper bound for Bisection Method
    real,parameter :: tol = 0.0001                  ! Tolerence for Bisection Method
    integer,parameter :: max_itter = 10000          ! Max iterations for Bisection Method

    integer, parameter :: tmax = (maxT-minT)/k      ! Number of time steps
    integer, parameter :: smax = ((maxS-minS)/h)-1  ! Number of space steps (nodes)

    ! Define all of the input and output variables
    real :: ResQ1                   ! TVD R-K variable (first split)
    real :: ResQ2                   ! TVD R-K variable (second split)
    real :: ResQ3                   ! TVD R-K variable (third split)
    real :: ResQ1_1                 ! TVD R-K variable (first split)
    real :: ResQ2_1                 ! TVD R-K variable (second split)
    real :: ResQ3_1                 ! TVD R-K variable (third split)
    real :: resQ1_mat(0:smax+1)     ! TVD R-K variable (first split)
    real :: resQ2_mat(0:smax+1)     ! TVD R-K variable (second split)
    real :: resQ3_mat(0:smax+1)     ! TVD R-K variable (third split)

    real :: uQ1_1                   ! TVD R-K variable (first split)
    real :: uQ2_1                   ! TVD R-K variable (second split)
    real :: uQ3_1                   ! TVD R-K variable (third split)
    real :: p0_mat(0:smax+1)        ! Array of initial pressure values
    real :: p_mat(0:smax+1)         ! Array of final pressure values
    real :: p_mat500(0:smax+1)      ! Array of pressure values at t = 0.5
    real :: p_mat1500(0:smax+1)     ! Array of pressure values at t = 1.0
    real :: rho0_mat(0:smax+1)      ! Array of initial density values
    real :: rho_mat(0:smax+1)       ! Array of final density values
    real :: rho_mat500(0:smax+1)    ! Array of density values at t = 0.5
    real :: rho_mat1500(0:smax+1)   ! Array of density values at t = 1.0
    real :: u0_mat(0:smax+1)         ! Array of initial velocity values
    real :: u_mat(0:smax+1)         ! Array of final velocity values
    real :: u_mat500(0:smax+1)      ! Array of velocity values at t = 0.5
    real :: u_mat1500(0:smax+1)     ! Array of velocity values at t = 1.0
    real :: Q1_mat(0:smax+1)        ! Array of temporary first split spatial values
    real :: Q2_mat(0:smax+1)        ! Array of temporary second split spatial values
    real :: Q3_mat(0:smax+1)        ! Array of temporary third split spatial values

    real :: pf_mat(0:smax+1)
    real :: rhof_mat(0:smax+1)
    real :: uf_mat(0:smax+1)
    real :: pf_mat500(0:smax+1)
    real :: rhof_mat500(0:smax+1)
    real :: uf_mat500(0:smax+1)
    real :: pf_mat1500(0:smax+1)
    real :: rhof_mat1500(0:smax+1)
    real :: uf_mat1500(0:smax+1)

    real :: x_mat(0:smax+1)         ! Array of x values
    real :: current_val             ! Temp variable to determine spatial value at current spatial step
    real :: t_val                   ! Temp variable to determine temporal value at current spatial step
    integer :: i                    ! Variable to be looped in do loop
    integer :: j                    ! Variable to be looped in do loop
    integer :: z                    ! Variable to be looped in do loop

    real :: c                       ! Fluid wave velocity and eigenvalue
    real :: M                       ! Fluid Mach number
    real :: A_star                  ! Reference fluid cross-sectional area (at max Mach number)
    real :: m_dot                   ! Fluid mass flow rate

    !!! Start the main code !!!

    ! Find some inital values first
    c = sqrt((gam*p_in)/rho_in)
    print*,"C = ",c
    M = u_in/c
    print*,"M = ",M
    A_star = (M*A(0.0)) / ( ((1+0.2*(M**2))/1.2)**3 )
    print*,"A_star = ",A_star
    m_dot = rho_in*u_in*A(0.0)
    print*,"m_dot = ",m_dot

    ! Create an array for the initial condition
    x_mat(0) = minS

    p0_mat(0) = p_in
    rho0_mat(0) = rho_in
    u0_mat(0) = u_in

    Q1_mat(0) = rho0_mat(0)
    Q2_mat(0) = rho0_mat(0)*u0_mat(0)
    Q3_mat(0) = rho0_mat(0)*E(p0_mat(0),rho0_mat(0),gam,u0_mat(0))

    do i = 1,smax
        current_val = x_mat(i-1)+h
        x_mat(i) = current_val

        ! Find the velocity based on the initial profile of the nozzle (using Bisection Method)
        M = M_guess(al,bh,current_val,A_star,gam,tol,max_itter)

        u0_mat(i) = M*c
        rho0_mat(i) = m_dot/(u0_mat(i)*A(current_val))
        p0_mat(i) = rho0_mat(i)/gam

        ! Fix these to represent the appropriate Q values
        ! Will probably have to change in the future
        Q1_mat(i) = rho0_mat(i)
        Q2_mat(i) = rho0_mat(i)*u0_mat(i)
        Q3_mat(i) = rho0_mat(i)*E(p0_mat(i),rho0_mat(i),gam,u0_mat(i))
    end do

    x_mat(smax+1) = maxS
    ! From the symmetry of the nozzle
    Q1_mat(smax+1) = Q1_mat(0)
    Q2_mat(smax+1) = Q2_mat(0)
    Q3_mat(smax+1) = Q3_mat(0)

    do i = 0,smax+1
        resQ1_mat(i) = Q1_mat(i)
        resQ2_mat(i) = Q2_mat(i)
        resQ3_mat(i) = Q3_mat(i)
        p_mat(i) = p0_mat(i)
        rho_mat(i) = rho0_mat(i)
        u_mat(i) = u0_mat(i)
    end do

    ! Move through time and update spatial components
    do i = 1,tmax

        t_val = minT+((i-1)*k)

        ! Two-stage TVD Runge-Kutta
        print*,"Two-stage TVD Runge-Kutta"
        do j = 1,smax

            if ((j == 1) .or. (j == smax)) then

                print*,p_mat(j)
                print*,rho_mat(j)
                print*,u_mat(j)

                ! First-Order
                ResQ1 = -1*( (( F(1.0,1.0,c,gam,rho_mat(j),u_mat(j))-F(1.0,1.0,c,gam,rho_mat(j-1),u_mat(j-1)) )/(2*h)) &
                    + (( F(-1.0,1.0,c,gam,rho_mat(j+1),u_mat(j+1))-F(-1.0,1.0,c,gam,rho_mat(j),u_mat(j)) )/(2*h)) &
                    + (-rho_mat(j)*u_mat(j)*(1/A(x_mat(j)))*((A(x_mat(j))-A(x_mat(j-1)))/h)) );

                ResQ2 = -1*( (( F(1.0,2.0,c,gam,rho_mat(j),u_mat(j))-F(1.0,2.0,c,gam,rho_mat(j-1),u_mat(j-1)) )/(2*h)) &
                    + (( F(-1.0,2.0,c,gam,rho_mat(j+1),u_mat(j+1))-F(-1.0,2.0,c,gam,rho_mat(j),u_mat(j)) )/(2*h)) &
                    + (-rho_mat(j)*(u_mat(j)**2)*(1/A(x_mat(j)))*((A(x_mat(j))-A(x_mat(j-1)))/h)) );

                ResQ3 = -1*( (( F(1.0,3.0,c,gam,rho_mat(j),u_mat(j))-F(1.0,3.0,c,gam,rho_mat(j-1),u_mat(j-1)) )/(2*h)) &
                    + (( F(-1.0,3.0,c,gam,rho_mat(j+1),u_mat(j+1))-F(-1.0,3.0,c,gam,rho_mat(j),u_mat(j)) )/(2*h)) &
                    + ((-u_mat(j)*(rho_mat(j)*E(p_mat(j),rho_mat(j),gam,u_mat(j))+p_mat(j)))*((A(x_mat(j))-A(x_mat(j-1)))/h)) );

            else

                ! Second-Order
                ResQ1 = -1*(((1/(2*h))*(3*F(1.0,1.0,c,gam,rho_mat(j),u_mat(j)))-4*F(1.0,1.0,c,gam,rho_mat(j-1),u_mat(j-1)) &
                    +F(1.0,1.0,c,gam,rho_mat(j-2),u_mat(j-2))) + ((1/(2*h))*(-F(-1.0,1.0,c,gam,rho_mat(j+2),u_mat(j+2)) &
                    +4*F(-1.0,1.0,c,gam,rho_mat(j+1),u_mat(j+1))-3*F(-1.0,1.0,c,gam,rho_mat(j),u_mat(j)))) &
                    + (-rho_mat(j)*u_mat(j)*(1/A(x_mat(j)))*((3*A(x_mat(j))-4*A(x_mat(j-1))+A(x_mat(j-2)))/(2*h))));

                ResQ2 = -1*(((1/(2*h))*(3*F(1.0,2.0,c,gam,rho_mat(j),u_mat(j)))-4*F(1.0,2.0,c,gam,rho_mat(j-1),u_mat(j-1)) &
                    +F(1.0,2.0,c,gam,rho_mat(j-2),u_mat(j-2))) + ((1/(2*h))*(-F(-1.0,2.0,c,gam,rho_mat(j+2),u_mat(j+2)) &
                    +4*F(-1.0,2.0,c,gam,rho_mat(j+1),u_mat(j+1))-3*F(-1.0,2.0,c,gam,rho_mat(j),u_mat(j)))) &
                    + (-rho_mat(j)*(u_mat(j)**2)*(1/A(x_mat(j)))*((3*A(x_mat(j))-4*A(x_mat(j-1))+A(x_mat(j-2)))/(2*h))));

                ResQ3 = -1*(((1/(2*h))*(3*F(1.0,3.0,c,gam,rho_mat(j),u_mat(j)))-4*F(1.0,3.0,c,gam,rho_mat(j-1),u_mat(j-1)) &
                    +F(1.0,3.0,c,gam,rho_mat(j-2),u_mat(j-2))) + ((1/(2*h))*(-F(-1.0,3.0,c,gam,rho_mat(j+2),u_mat(j+2)) &
                    +4*F(-1.0,3.0,c,gam,rho_mat(j+1),u_mat(j+1))-3*F(-1.0,3.0,c,gam,rho_mat(j),u_mat(j)))) &
                    + ((-u_mat(j)*(rho_mat(j)*E(p_mat(j),rho_mat(j),gam,u_mat(j))+p_mat(j))) &
                    *((3*A(x_mat(j))-4*A(x_mat(j-1))+A(x_mat(j-2)))/(2*h))));

            end if

            uQ1_1 = Q1_mat(j)+k*ResQ1
            uQ2_1 = Q2_mat(j)+k*ResQ2
            uQ3_1 = Q3_mat(j)+k*ResQ3

            ! Residual is zero for this problem since it doesn't depend on time
            ResQ1_1 = 0!ResQ1
            ResQ2_1 = 0!ResQ2
            ResQ3_1 = 0!ResQ3

            Q1_mat(j) = (0.5)*(Q1_mat(j)+uQ1_1+k*ResQ1_1)
            rho_mat(j) = Q1_mat(j)
            Q2_mat(j) = (0.5)*(Q2_mat(j)+uQ2_1+k*ResQ2_1)
            u_mat(j) = Q2_mat(j)/rho_mat(j)
            Q3_mat(j) = (0.5)*(Q3_mat(j)+uQ3_1+k*ResQ3_1)
            p_mat(j) = ((Q3_mat(j)/rho_mat(j))-((u_mat(j)**2)/2))*rho_mat(j)*(gam-1)

        end do

        ! Should change from Q1,Q2,Q3 to rho,u,p to get the correct output values
        p_mat(0) = p_in
        p_mat(smax+1) = p_out

        rho_mat(0) = rho_in
        rho_mat(smax+1) = rho_out

        u_mat(0) = u_in
        u_mat(smax+1) = u_out

        ! Should change from Q1,Q2,Q3 to rho,u,p to get the correct output values
        do z = 0,smax+1
            !Q1_mat(z) = rho_mat(z)
            !Q2_mat(z) = rho_mat(z)*u_mat(z)
            !Q3_mat(z) = rho_mat(z)*E(p_mat(z),rho_mat(z),gam,u_mat(z))
            pf_mat(z) = p_mat(z)
            rhof_mat(z) = rho_mat(z)
            uf_mat(z) = u_mat(z)
            !print*,rhof_mat(z)*uf_mat(z)*A(x_mat(z))
            if (t_val == 0.5) then
                !Q1_mat500(z) = rho_mat(z)
                !Q2_mat500(z) = rho_mat(z)*u_mat(z)
                !Q3_mat500(z) = rho_mat(z)*E(p_mat(z),rho_mat(z),gam,u_mat(z))
                pf_mat500(z) = p_mat(z)
                rhof_mat500(z) = rho_mat(z)
                uf_mat500(z) = u_mat(z)
            end if
            if (t_val == 1.0) then
                !Q1_mat1500(z) = rho_mat(z)
                !Q2_mat1500(z) = rho_mat(z)*u_mat(z)
                !Q3_mat1500(z) = rho_mat(z)*E(p_mat(z),rho_mat(z),gam,u_mat(z))
                pf_mat1500(z) = p_mat(z)
                rhof_mat1500(z) = rho_mat(z)
                uf_mat1500(z) = u_mat(z)
            end if
        end do

    end do

	! Export the corresponding data to al '.dat' file for future plotting
    open(unit=48, file='Project2_G_test.dat')
    do i = 0,smax+1
        current_val = x_mat(i)
        print*,"(",current_val,",",p0_mat(i),",",pf_mat(i),",",rho0_mat(i),",",rhof_mat(i),",",u0_mat(i),",",uf_mat(i),")"
        write(48,*) current_val,p0_mat(j),pf_mat500(i),pf_mat1500(i),pf_mat(i) &
            ,rho0_mat(j),rhof_mat500(i),rhof_mat1500(i),rhof_mat(i) &
            ,u0_mat(j),uf_mat500(i),uf_mat1500(i),uf_mat(i)
    end do
    close(48)

	! I am plotting with Gnuplot so this calls that program for plotting... must have the file named that for it to work
    call system('gnuplot -p Project2_G.plt')

    ! Import all of the required outside functions...

contains

    ! Determine the specific total energy
    function E(p,rho,gam,u) result(y)

        implicit none

        real, intent(in) :: p       ! Input = fluid pressure
        real, intent(in) :: rho     ! Input = fluid density
        real, intent(in) :: gam     ! Input = fluid specific weight
        real, intent(in) :: u       ! Input = fluid velocity
        real :: y                   ! Output = fluid specific total energy

        y = (p/(rho*(gam-1)))+((u**2)/2)

    end function E

    ! Determine the cross-sectional area of the nozzle
    ! _____     _____
    !      \___/
    ! x->   ___
    ! _____/   \_____
    !
    function A(x) result(y)

        implicit none

        real, intent(in) :: x       ! Input = position along the nozzle
        real :: y                   ! Output = cross-sectional area at the current position

        y = 1+2.2*((x-1.5)**2)

    end function A

    ! Determine the spatial split values for both positive and negative
    function F(sgn,level,c,gam,rho,u) result(y)

        implicit none

        real, intent(in) :: sgn     ! Input = determines whether to calculate F- or F+
        real, intent(in) :: level   ! Input = determines which equations to calculate: Q1,Q2,Q3
        real, intent(in) :: gam     ! Input = fluid specific weight
        real, intent(in) :: rho     ! Input = fluid density
        real, intent(in) :: u       ! Input = fluid velocity
        real, intent(in) :: c       ! Input = wave speed (to be determined by the input parameters)
        real :: y                   ! Output = value of the split spatial values

        ! Should always be 1
        !c = sqrt((gam*p)/rho)
        !print*,c

        ! Negative
        if (sgn == -1.0) then
            if (level == 1.0) then
                y = ((rho*(u-c))/(2*gam))
                !print*,rho
                !print*,u
                !print*,gam
                !print*,y
                !print*,c
            else if (level == 2.0) then
                y = ((rho*(u-c))/(2*gam))*(u-c)
            else
                y = ((rho*(u-c))/(2*gam))*( (((u-c)**2)/2) * (((3-gam)/(gam-1))*((c**2)/2)) )
            end if
        ! Positive
        else
            if (level == 1.0) then
                y = (rho/(2*gam))*((((2*gam)-1)*u)+c)
            else if (level == 2.0) then
                y = (rho/(2*gam))*( (2*(gam-1)*(u**2))+((u+c)**2) )
            else
                y = (rho/(2*gam))*( ((gam-1)*(u**3))+(((u+c)**3)/2) &
                    * (((3-gam)/(gam-1))*((c**2)/2)*(u+c)) )
            end if
        end if
        !print*,y

    end function F

    ! Determine the Mach number at the desired spatial position along the nozzle
    function M_val(valu,x,A_star,gam) result(y)

        implicit none

        real, intent(in) :: valu    ! Input = the guessed Mach number
        real, intent(in) :: x       ! Input = determines the spatial position along the nozzle
        real, intent(in) :: A_star  ! Input = reference area (has the highest Mach number)
        real, intent(in) :: gam     ! Input = fluid specific weight
        real :: y                   ! Output = value of the difference between guessed and real Mach

        y = (A_star/A(x)) - ((((gam+1)/2)**((gam+1)/(2*(gam-1)))) &
        * (valu/((1+((gam-1)/2)*(valu**2))**((gam+1)/(2*(gam-1))))))

    end function M_val

    ! Determine the Mach number at the desired spatial position along the nozzle
    function M_guess(al,bh,x,A_star,gam,tol,max_itter) result(y)

        implicit none

        real, intent(in) :: al              ! Input = the lower bound to start guessing
        real, intent(in) :: bh              ! Input = the upper bound to start guessing
        real, intent(in) :: x               ! Input = the current spatial value along the nozzle
        real, intent(in) :: A_star          ! Input = reference area (has the highest Mach number)
        real, intent(in) :: gam             ! Input = fluid specific weight
        real, intent(in) :: tol             ! Input = the tolerance where the Mach number is close enough
        integer, intent(in) :: max_itter    ! Input = determines the spatial position along the nozzle
        integer :: itter = 1                ! Temp = make sure the function doesn't cause an infinite loop
        real :: c_mid                       ! Temp = determine the new midpoint value
        real :: new_al                      ! Temp = determine the new lower bound value
        real :: new_bh                      ! Temp = determine the new upper bound value
        real :: y                           ! Output = value of the difference between guessed and real Mach

        new_al = al
        new_bh = bh

        do while (itter .le. max_itter)

            c_mid = (new_al+new_bh)/2
            if ( (M_val(c_mid,x,A_star,gam) == 0) .or. (((new_bh-new_al)/2) .lt. tol) ) then
                y = c_mid
                EXIT
            end if

            itter = itter + 1
            if ( ((M_val(c_mid,x,A_star,gam) .lt. 0) .and. (M_val(new_al,x,A_star,gam) .lt. 0)) &
                .or. ((M_val(c_mid,x,A_star,gam) .gt. 0) .and. (M_val(new_al,x,A_star,gam) .gt. 0)) ) then
                new_al = c_mid
            else
                new_bh = c_mid
            end if

        enddo

    end function M_guess

end program plotting
