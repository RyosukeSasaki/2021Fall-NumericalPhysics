module differential
    implicit none
     
    contains
    !runge_kuttaの計算を1ステップ進める
    !arg:微分方程式; n:連立する数; init:初期値; t_begin:計算開始; t_end:計算終了; tau:刻み幅;
    !const:定数(optional); boundary:境界条件(optional)
    function runge_kutta(arg, n, init, t_begin, tau, const, boundary)
        implicit none
        interface
            !微分方程式
            function arg(t, x, n, const)
                INTEGER, INTENT(IN) :: n
                DOUBLE PRECISION :: arg(n)
                DOUBLE PRECISION, INTENT(IN) :: x(:), t
                DOUBLE PRECISION, OPTIONAL :: const(:)
            end function arg
            !境界条件の計算関数
            function boundary(x, n)
                implicit none
                INTEGER, INTENT(IN) :: n
                DOUBLE PRECISION, INTENT(IN) :: x(:)
                DOUBLE PRECISION, OPTIONAL :: boundary(n)
            end function boundary
        end interface
        INTEGER :: i
        !段数
        INTEGER, PARAMETER :: order = 4
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION, INTENT(in) :: init(:), tau
        DOUBLE PRECISION, INTENT(INOUT) :: t_begin
        DOUBLE PRECISION :: runge_kutta(n), x(n), t, s(n), delta(n)
        DOUBLE PRECISION, OPTIONAL :: const(:)
        DOUBLE PRECISION :: a(order), b(order)
        a(1)=0d0; a(2)=0.5d0; a(3)=0.5d0; a(4)=1d0
        b(1)=1d0/6d0; b(2)=1d0/3d0; b(3)=1d0/3d0; b(4)=1d0/6d0
        x(:) = init(:)
        t = t_begin

        s = 0; delta = 0
        do i = 1, order
            !optionalの定数が与えられている場合はそれを含む計算を実行
            if (PRESENT(const)) then 
                delta(:) = arg(t+a(i)*tau, x(:) + a(i) * delta, n, const)*tau
            else
                delta(:) = arg(t+a(i)*tau, x(:) + a(i) * delta, n)*tau
            end if
            s(:) = s(:) + b(i) * delta(:)
        end do
        x(:) = x(:) + s(:)
        t = t + tau
        t_begin = t
        !optionalの境界条件が与えられている場合はそれを考慮
        if (PRESENT(boundary)) then
            x(:) = boundary(x, n)
        end if
        runge_kutta(:) = x(:)
    end function runge_kutta
end module differential