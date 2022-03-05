module config
    implicit none
    public

    !Segment boundaries
    real(8), parameter :: st = 0.5d0, ed = 6.d0
    !Threads nubmber
    integer, parameter :: thr = 4
    !Machine zero
    real(8), parameter :: eps = 1d-4
    !real(8), parameter :: pi = 4*atan(1.d0)
end module config

program output
    use config
    use omp_lib
    implicit none

    procedure(ftst), pointer :: fptr

    fptr => ftst

    call bisection(fptr, st, ed, eps)
    call secant(fptr, st, ed, eps)
    call muller(fptr, st, ed, eps)
    call newton(fptr, st, ed, eps)
    call iterations(fptr, st, ed, eps)
    call vegstein(fptr, st, ed, eps)

contains

    pure real(8) function ftst(x)
        implicit none
        real(8), intent(in) :: x
        ftst = (x-3)**3-(x-3)
        !ftst = cos(x)-x+1
    end function

    subroutine bisection(fptr, a, b, prec)
        !Please, ensure that the function is monotoneous, otherwise
        !the criteria f(a)*f(b)<0 doesn't really make that much sense
        !The a and b aren't initial guesses, those are the borders
        implicit none

        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: a, b, prec
        real(8) :: bt, md, tp
        integer :: j

        if (fptr(a) * fptr(b) .gt. 0) then
            print *, "Bisection error: root criteria isn't met"
            return
        end if

        bt = a
        tp = b
        md = (tp + bt) / 2
        j = 0

        do while ((abs(fptr(md)) .gt. prec).and.((tp-bt) .gt. prec))
        !do while (abs(fptr(md)) .gt. prec)
            md = (tp + bt) / 2
            if (fptr(bt) * fptr(md) .gt. 0) then
                bt = md
            else if (fptr(bt) * fptr(md) .lt. 0) then
                tp = md
            end if
            j = j + 1
        end do

        print *, "Bisection gives x =", md, "as a root, in", j , "steps"

    end subroutine bisection

    subroutine secant(fptr, a, b, prec)
        implicit none

        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: a, b, prec
        real(8) :: bt, md, tp
        integer :: j

        bt = a
        tp = b
        md = tp - fptr(tp)*(tp-bt)/(fptr(tp)-fptr(bt))
        j = 0

        do while (abs(fptr(md)) .gt. prec)
            md = tp - fptr(tp)*(tp-bt)/(fptr(tp)-fptr(bt))
            bt = tp
            tp = md
            j= j + 1
        end do

        print *, "Secant gives    x =", md, "as a root, in", j , "steps"
    end subroutine secant

    subroutine muller(fptr, a, b, prec)
        implicit none

        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: a, b, prec
        real(8) :: bt, md, tp, w, st
        integer :: j

        if (fptr(a) * fptr(b) .gt. 0) then
            print *, "Muller error: root criteria isn't met"
            return
        end if

        !The algorithm requires 3 initial guesses
        bt = a
        tp = b
        md = (a + b) / 2

        j = 0
        do while (abs(fptr(md)) .gt. prec)
            st = bt
            bt = md
            md = tp
            w = div(fptr, md, bt) + div(fptr, md, st) - div(fptr, bt, st)
            if (abs(w + cnd(fptr, w, md, bt, st)) > abs(w - cnd(fptr, w, md, bt, st))) then
                tp = md - 2*fptr(md)/(w + cnd(fptr, w, md, bt, st))
            else
                tp = md - 2*fptr(md)/(w - cnd(fptr, w, md, bt, st))
            end if
            j = j +1
        end do

        print *, "Muller gives    x =", md, "as a root, in", j , "steps"
    end subroutine muller

    real(8) function div(fptr, x1, x2)
        !First divided difference
        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: x1, x2
        div = (fptr(x2)-fptr(x1))/(x2-x1)
    end function div

    real(8) function ddiv(fptr, x1, x2, x3)
        !Second divided difference
        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: x1, x2, x3
        ddiv = (div(fptr, x2, x3)-div(fptr, x1, x2))/(x3-x1)
    end function ddiv

    real(8) function cnd(fptr, w, md, bt, st)
        !To simplify the if-statement in muller()
        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: w, md, bt, st
        cnd = sqrt(w**2-4*fptr(md)*ddiv(fptr, md, bt, st))
    end function cnd

    subroutine newton(fptr, a, b, prec)
        !Simple iterations with lam = 1/f'(cr)
        implicit none

        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: a, b, prec
        real(8) :: cr
        integer :: j

        !Algorithm requires 1 initial guess
        cr = (a + b) / 2

        j = 0
        do while (abs(fptr(cr)) .gt. prec)
            cr = cr - fptr(cr)/diff(fptr, cr, prec)
            j = j + 1
        end do

        print *, "Newton gives    x =", cr, "as a root, in", j , "steps"
    end subroutine newton

    real(8) function diff(fptr, arg, step)
        !First derivative
        procedure(ftst), pointer :: fptr
        real(8), intent(in):: arg, step
        diff = (fptr(arg+step) - fptr(arg-step))/(2*step)
    end function diff

    subroutine iterations(fptr, a, b, prec)
        implicit none

        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: a, b, prec
        real(8) :: cr
        integer :: j

        !Algorithm requires 1 initial guess
        cr = (a + b) / 2

        if (abs(diff(fptr, cr, prec)) .ge. 1) then
            print *, "Iterations error: derivative criteria isn't met"
            return
        end if

        j = 0
        do while (abs(fptr(cr)) .gt. prec)
            cr = cr - lam(cr)*fptr(cr)
            j = j + 1
        end do

        print *, "Iterations give x =", cr, "as a root, in", j , "steps"
    end subroutine iterations

    real(8) function lam(arg)
        real(8), intent(in):: arg
        lam = -1 !Gives the simplest iterations method
    end function lam

    subroutine vegstein(fptr, a, b, prec)
        implicit none

        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: a, b, prec
        real(8) :: tp, bt, st
        integer :: j

        bt = a
        tp = b

        j = 0
        do while (abs(fptr(tp)) .gt. prec)
            st = tp
            tp = tp - (fptr(tp)*(tp-bt))/(fptr(tp)-fptr(bt))
            bt = st
            j = j + 1
        end do

        print *, "Vegstein gives  x =", tp, "as a root, in", j , "steps"
    end subroutine vegstein

end program output
