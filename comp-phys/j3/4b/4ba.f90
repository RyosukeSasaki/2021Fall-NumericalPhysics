program main
    use differential
    use consts
    implicit none
    INTEGER, PARAMETER :: idim=4096
    COMPLEX(KIND(0d0)) :: cp(idim)
    COMPLEX(KIND(0d0)), PARAMETER :: ci = (0d0, 1d0)
    INTEGER :: jmax,nmax,j1,j2,p,q,s
    DOUBLE PRECISION :: sigma,k0,dx,w,t
    DOUBLE PRECISION :: r(idim),x(idim),V(idim)=0d0
    
    do s=2,2
        sigma = sigma_array(s)
        do q=1,1
            k0 = k0_array(q)
            nmax = int(3d0/2d0*L/k0/dt)
            write(*,'(''# sig   = '',e18.8e3)') sigma
            write(*,'(''# k0    = '',e18.8e3)') k0
            do p=1,7
                jmax = 2*int(L)*x_division(p)
                dx = 2.0d0*L/dble(jmax)
                w = 0.5d0/dx**2
                call main_roop()
            end do
        end do
    enddo


    contains
    subroutine main_roop()
        implicit none
        INTEGER :: j,n,n1
        DOUBLE PRECISION, PARAMETER :: eps = 2d0**(-52)
        DOUBLE PRECISION :: sum,rnf
        !========================================================-
        ! generate grid points
        !========================================================-
        do j=1,jmax+1
            x(j)=dx*dble(j-1)-L
        enddo
        !========================================================-
        ! generate potential array
        !========================================================-
        do j=2,jmax
            if (abs(x(j)-x1).le.eps) then
                V(j) = V0*0.5d0
            else if (abs(x(j)-x2).le.eps) then
                V(j) = V0*0.5d0
            else if ((x(j).gt.x1).and.(x(j).lt.x2)) then
                V(j) = V0
            else
                V(j) = 0d0
            end if
        end do
        j = 0
        do while (V(j).eq.0d0)
            j = j + 1
        end do
        j1 = j
        j = jmax
        do while (V(j).eq.0d0)
            j = j - 1
        end do
        j2 = j
        !========================================================-
        ! inital condition
        !========================================================-
        sum = 0d0
        n = 0
        n1 = 0
        do j=2,jmax
            cp(j) = cdexp(ci*k0*x(j)-((x(j)-x0)/(2d0*sigma))**2d0)
            sum = sum+cp(j)*dconjg(cp(j))
        enddo
        cp(1) = (0d0,0d0)
        cp(jmax+1) = (0d0,0d0)
        rnf = 1.0d0/dsqrt(sum*dx)
        do j=2,jmax
            cp(j)=rnf*cp(j)
        enddo
        !call output_xav_xs()
        !========================================================-
        ! time evolution
        !========================================================-
        do n=1,nmax
            cp = runge_kutta(f, idim, cp, t, dt)
            !call output_xav_xs()
        end do
        do j=2,jmax
            r(j) = cp(j)*dconjg(cp(j))
        enddo
        call output_wave()
        !call output_reflect_ratio()
    end subroutine main_roop
    
    subroutine output_reflect_ratio()
        implicit none
        DOUBLE PRECISION :: Pleft, Pcenter, Pright
        INTEGER :: j
        Pleft = 0d0
        Pcenter = 0d0
        Pright = 0d0
        do j=2, j1-1
            Pleft = Pleft + r(j)*dx
        end do
        do j=j1, j2
            Pcenter = Pcenter + r(j)*dx
        end do
        do j=j2+1, jmax
            Pright = Pright + r(j)*dx
        end do
        write(*,'(5e18.8e3)') dx, Pleft, Pcenter, Pright, Pleft+Pcenter+Pright
    end subroutine output_reflect_ratio
    
    subroutine output_wave()
        implicit none
        INTEGER :: j
        do j=1, jmax+1
            write(*,'(4e18.8e3)') x(j), r(j), real(cp(j)), aimag(cp(j))
        end do
    end subroutine output_wave

    subroutine output_xav_xs()
        implicit none
        INTEGER :: j
        DOUBLE PRECISION :: xav,xs
        xav=0.0d0
        do j=2,jmax
            xav = xav+x(j)*r(j)*dx
        enddo
        xs=0.0d0
        do j=2,jmax
            xs = xs+((x(j)-xav)**2)*r(j)*dx
        enddo
        write(*,'(3e18.8e3)') t,xav,xs
    end subroutine output_xav_xs
    
    function f(t, x, n, const)
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION :: t
        COMPLEX(KIND(0d0)) :: f(n)
        COMPLEX(KIND(0d0)), INTENT(IN) :: x(:)
        COMPLEX(KIND(0d0)), OPTIONAL :: const(:)
        INTEGER :: j
        do j=2,jmax
            f(j)=((x(j+1)-2.0d0*x(j)+x(j-1))*w - V(j)*x(j))*ci
        enddo
        f(1) = (0d0,0d0); f(jmax + 1) = (0d0,0d0)
    end function f
end program main