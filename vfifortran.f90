program vfi
  implicit none
  
  real :: start, finish
  integer, parameter :: A = 10
  real, parameter :: alpha = 0.5
  real, parameter :: beta = 0.9
  real, parameter :: tol = 1e-6
  integer, parameter :: maxiter = 1000

  real :: kmin, kmax
  integer :: prec
  real, allocatable :: kgrid(:)
  real, allocatable :: gk(:)

  real, allocatable :: Vk0(:) , Vk(:) , Vkprim(:), squared_valuef(:)
  real, allocatable :: value_array(:, :)
  integer :: n_iter = 0
  integer :: i, iprim, col
  real :: k, kprim, c, kstar
  real :: norm

  character(len=*), parameter :: OUT_FILE = 'datavfi.txt'
  character(len=*), parameter :: PLT_FILE = 'plotvfi.plt'
  integer::fu, m

  call cpu_time(start)


  kmin=1
  kmax=25
  prec=1000
  
  allocate(kgrid(prec))
  allocate(gk(prec))
  allocate(Vk0(prec))
  allocate(Vk(prec))
  allocate(Vkprim(prec))
  allocate(squared_valuef(prec))
  allocate(value_array(prec, prec))

  norm = 1e5
  Vk0(:) = 1.0

  call linspace(kmin, kmax, prec, kgrid)
  call linspace(kmin, kmax, prec, gk)


  Vk = Vk0
  
  do while (n_iter < maxiter .and. norm > tol)
    do iprim = 1, prec
      kprim = kgrid(iprim)
      do i = 1, prec
        k = kgrid(i)
        c = A*(k**alpha) - kprim
        if (c > 0) then
          value_array(i, iprim) = log(c) + beta*Vk(iprim)
        else
          value_array(i, iprim) = -1e6
        end if
      end do
    end do
    

    do col=1, prec
      gk(col) = kgrid(maxloc(value_array(col, :), dim=1))
      Vkprim(col) = maxval(value_array(col, :))
    end do

    call pow_vector(Vkprim - Vk, prec, 2.0, squared_valuef)
    norm = sum(squared_valuef)
    Vk = Vkprim
    n_iter = n_iter + 1
    print*, "Iteration: ", n_iter, " norm: ", norm
  end do

  kstar = kgrid(minloc(abs(gk - kgrid), dim=1)+1)
  print*, "The Steady state value of capital is: ", kstar


  open (action='write', file=OUT_FILE, newunit=fu, status='replace')
  do m=1, prec 
    write(fu, *) kgrid(m), gk(m)
  end do
  close(fu)

  call execute_command_line('gnuplot -c plotvfi.plt')
  
  call cpu_time(finish)
  

  print '("VFI in FORTRAN executed in ",f6.3," seconds.")', finish-start
  print*, "Press enter to leave"
  read*
end program vfi


subroutine linspace(a, b, s, linvec)  
  implicit none
  real, intent(in) :: a, b
  integer, intent(in) :: s
  real, dimension(s)::vec
  real, dimension(s), intent(out)::linvec 
  integer :: i
  
  do i = 0, s-1
    vec(i+1) = a + i*(b-a)/(s-1)
  end do       
  
  linvec = vec
  return
end


subroutine pow_vector(a, n, p, ap)
  implicit none
  real, intent(in), dimension(n) :: a
  integer, intent(in) :: n
  real, intent(in) :: p
  integer::i 
  real, intent(out), dimension(n) :: ap

  do i = 1, n 
    ap(i) = a(i)**p 
  end do 
  return
end 

