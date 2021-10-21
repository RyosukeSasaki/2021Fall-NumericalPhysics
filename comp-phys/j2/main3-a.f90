program main
    use differential
    implicit none

    contains
    function v(t, x)
        implicit none
        DOUBLE PRECISION :: v
        DOUBLE PRECISION, INTENT(IN) :: t, x
        v = 1
    end function v
end program main