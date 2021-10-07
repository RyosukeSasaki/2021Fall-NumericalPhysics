program main
    integer(16) :: n
    double precision :: OnePlusEpsilon
    !real :: OnePlusEpsilon
    do n = 1, huge(n)
        OnePlusEpsilon = 1d0 + 2d0**(-n)
        if (OnePlusEpsilon == 1d0) exit
    end do
    
    write(*, '(a)', advance='no') "n = "
    write(*, '(i0)', advance='no') (n - 1)
    write(*, '(a)', advance='no') ", epsilon = "
    write(*, '(E8.2)') (2d0**(-(n - 1)))

end