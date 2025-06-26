module octree_module
  implicit none
  real, parameter :: G = 6.67430e-11
  real, parameter :: softening = 1.0e-5
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nq = 1000  ! Number of samples between q = 0 and q = 2
  real(dp), allocatable :: w_table(:), dw_table(:)
  real(dp), parameter :: dq = 2.0_dp / nq, smoothing = 10.0_dp


  !types here act like dictionaries with "sets" of arrays
  type :: particle
    real(dp) :: mass, density, internal_energy, pressure
    real(dp) :: internal_energy_rate
    real(dp), dimension(3) :: position
    real(dp), dimension(3) :: velocity
    real(dp), dimension(3) :: acceleration
  end type particle

  type :: branch
    real(dp), dimension(3) :: center
    real(dp) :: size
    integer :: n_particles
    type(Particle), allocatable :: particles(:)
    real(dp) :: mass_total
    real(dp), dimension(3) :: mass_center
    type(branch), allocatable :: children(:)
  end type branch

  contains
  subroutine init_kernel_table()
    integer :: i
    real(dp) :: q

    allocate(w_table(0:nq))
    allocate(dw_table(0:nq))

    do i = 0, nq
      q = i * dq
      if (q >= 0.0_dp .and. q <= 1.0_dp) then
        w_table(i) = 1.0_dp - 1.5_dp*q**2 + 0.75_dp*q**3
        dw_table(i) = -3.0_dp*q + 2.25_dp*q**2
      else if (q > 1.0_dp .and. q <= 2.0_dp) then
        w_table(i) = 0.25_dp * (2.0_dp - q)**3
        dw_table(i) = -0.75_dp * (2.0_dp - q)**2
      else
        w_table(i) = 0.0_dp
        dw_table(i) = 0.0_dp
      end if
    end do
  end subroutine init_kernel_table

  subroutine lookup_kernel(r, hi, Wi, dWi)
    real(dp), intent(in) :: r, hi
    real(dp), intent(out) :: Wi, dWi
    integer :: i
    real(dp) :: alpha, qi

    qi = r / hi
    if (qi >= 0.0_dp .and. qi <= 2.0_dp) then
      i = int(qi / dq)
      alpha = (qi - i * dq) / dq
      Wi = (1.0_dp - alpha) * w_table(i) + alpha * w_table(i+1)
      dWi = (1.0_dp - alpha) * dw_table(i) + alpha * dw_table(i+1)
    else
      Wi = 0.0_dp
      dWi = 0.0_dp
    end if
  end subroutine lookup_kernel

  !All this does is check an volume for particles, then subdivide into 8
  !and allocate particles to the right octant, nothing more 
  recursive subroutine build_tree(node, depth, max_particles)
    implicit none
    type(branch), intent(inout) :: node
    integer, intent(in) :: depth, max_particles

    integer :: i, j, n
    real(dp), dimension(3) :: offset
    type(Particle), allocatable :: child_particles(:)
    integer :: child_index
    integer :: p

    !compute total mass and center of mass
    node%mass_total = 0.0_dp
    node%mass_center = 0.0_dp
    do i = 1, size(node%particles)
      node%mass_total = node%mass_total + node%particles(i)%mass
      node%mass_center = node%mass_center + node%particles(i)%mass * node%particles(i)%position
    end do
    if (node%mass_total > 0.0_dp) then
      node%mass_center = node%mass_center / node%mass_total
    else
      node%mass_center = node%center
    end if

    !base case: do not subdivide further
    if (size(node%particles) <= max_particles .or. depth == 0) return

    !allocate 8 children
    allocate(node%children(8))
    do i = 1, 8
      node%children(i)%size = node%size / 2.0_dp
      node%children(i)%n_particles = 0
      allocate(node%children(i)%particles(0))
      node%children(i)%center = node%center
      !octant shifts: ±x ±y ±z combinations (might be faster to write out the full 8 combs)
      do j = 1, 3
        if (btest(i - 1, j - 1)) then
          offset(j) = 0.25_dp * node%size
        else
          offset(j) = -0.25_dp * node%size
        end if
      end do
      node%children(i)%center = node%center + offset
    end do

    !distribute particles to children
    do p = 1, size(node%particles)
      child_index = 1
      do j = 1, 3
        if (node%particles(p)%position(j) > node%center(j)) then
          child_index = ibset(child_index - 1, j - 1) + 1
        end if
      end do

      n = size(node%children(child_index)%particles)
      call move_alloc(node%children(child_index)%particles, child_particles)
      allocate(node%children(child_index)%particles(n + 1))
      if (n > 0) node%children(child_index)%particles(1:n) = child_particles
      node%children(child_index)%particles(n + 1) = node%particles(p)
      node%children(child_index)%n_particles = node%children(child_index)%n_particles + 1
    end do


    !recursively build children
    do i = 1, 8
      if (node%children(i)%n_particles > 0) then
        call build_tree(node%children(i), depth - 1, max_particles)
      end if
    end do
  end subroutine build_tree

  ! this is just a standard force calc for an octree, just checks if a node satisfies a theta condition to be treated as a single point.
  recursive subroutine navigate_tree(node, body, theta)
  implicit none
  type(branch), intent(in) :: node
  type(particle), intent(inout) :: body(:)
  real(dp), intent(in) :: theta

  real(dp) :: dist, d2
  real(dp), dimension(3) :: direction
  integer :: i, k

  do i = 1, size(body)
    direction = body(i)%position - node%mass_center
    d2 = sum(direction**2) + softening**2
    dist = sqrt(d2)

    if ((node%size / dist) < theta .or. .not. allocated(node%children)) then
      if (node%mass_total > 0.0_dp .and. dist > 0.0_dp) then
        body(i)%acceleration = body(i)%acceleration - (G * node%mass_total * direction / (dist**3))
      end if

    else if (allocated(node%children)) then
      do k = 1, size(node%children)
        if (node%children(k)%n_particles > 0) then
          call navigate_tree(node%children(k), body(i:i), theta)
        end if
      end do
    end if

  end do
  end subroutine navigate_tree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PAST THIS POINT THE STUFF WORKS FINE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_SPH(root, body)
  implicit none
  type(branch), intent(in) :: root
  type(particle), intent(inout) :: body(:)
  integer :: i

  do i = 1, size(body)
    call SPH_tree_search(root, body(i))
  end do
end subroutine get_SPH


recursive subroutine SPH_tree_search(node, body)
  implicit none
  type(branch), intent(in) :: node
  type(particle), intent(inout) :: body
  real(dp), dimension(3) :: over_dr, nr, dWj, vij
  real(dp) :: Wj, dr, mj, dWj_mag, vdotW
  integer :: j
  logical :: has_children

  has_children = allocated(node%children)
  !overlaping distance
  over_dr = (body%position - node%center)  

  if (node%n_particles > 1 .and. all(over_dr < (2*smoothing + node%size/2)) .and. has_children) then
    do j = 1, size(node%children)
      call SPH_tree_search(node%children(j), body)
    end do

  else if (node%n_particles == 1 .and. all(over_dr < (2*smoothing + node%size/2))) then
    !THIS IS JUST ALL THE VECTORS THAT MATTER
    nr = (body%position - node%particles(1)%position)
    vij = (body%velocity - node%particles(1)%velocity)
    dr = sqrt(sum(nr**2)) !might need to switch i and j
    if (dr == 0.0_dp) return
    nr = nr / dr

    call lookup_kernel(dr, smoothing,  Wj, dWj_mag)
    Wj = Wj / (3.14159265359_dp * smoothing**3) 
    dWj = nr * dWj_mag / (3.14159265359_dp * smoothing**4)
    mj = node%particles(1)%mass

    
    !now do forces and shit
    !pressure force
    body%acceleration = body%acceleration - mj * ((body%pressure/(body%density * body%density)) + (node%particles(1)%pressure / (node%particles(1)%density * node%particles(1)%density))) * dWj
    !drho/dt
    vdotW = sum(vij * dWj)
    !rate change in internal energy (might be faster to take out p/rho term and do as a bulk calc after)
    body%internal_energy_rate = body%internal_energy_rate + (body%pressure/body%density)*mj*(vdotW)
    return
  end if
end subroutine SPH_tree_search

!call to find density of all bodies
subroutine get_density(root, body)
  implicit none
  type(branch), intent(inout) :: root
  type(particle), intent(inout) :: body(:)
  integer :: i
  logical :: has_children

  has_children = allocated(root%children)
  do i = 1, size(body)
    body(i)%density = 0.0_dp
    call density_tree_search(root, body(i))
  end do
end subroutine get_density


recursive subroutine density_tree_search(node, body)
  implicit none
  type(branch), intent(inout) :: node
  type(particle), intent(inout) :: body
  real(dp), dimension(3) :: over_dr, nr
  real(dp) :: Wj, dr, mj, dWj_mag
  integer :: j
  logical :: has_children

  has_children = allocated(node%children)
  !overlaping distance
  over_dr = (body%position - node%center)  

  if (node%n_particles > 1 .and. all(over_dr < (2*smoothing + node%size/2)) .and. has_children) then
    do j = 1, size(node%children)
      call density_tree_search(node%children(j), body)
    end do

  else if (node%n_particles == 1 .and. all(over_dr < (2*smoothing + node%size/2))) then
    nr = (body%position - node%particles(1)%position)
    dr = sqrt(sum(nr**2)) !might need to switch i and j
    call lookup_kernel(dr, smoothing,  Wj, dWj_mag)
    Wj = Wj / (3.14159265359_dp * smoothing**3) 
    mj = node%particles(1)%mass
      
    body%density = body%density + mj * Wj
    node%particles(1)%density = node%particles(1)%density + mj * Wj
    return
  end if
end subroutine density_tree_search
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simulate(bodies)
  implicit none
  type(particle), intent(inout) :: bodies(:)
  type(branch), allocatable :: root
  real(dp) :: t, dt, end_time
  integer :: i

  end_time = 1
  dt = 0.1
  do while (t < end_time)
    print *, bodies(1)%position(1),bodies(1)%position(2),bodies(1)%position(3)
    allocate(root)
    ! Initialize root node
    root%center = [150.0_dp, 150.0_dp, 150.0_dp]
    root%size = 10000.0_dp
    root%n_particles = size(bodies)
    allocate(root%particles(size(bodies)))
    root%particles = bodies

    call build_tree(root, 1000, 1)
    call get_density(root, bodies)
    do i = 1, root%n_particles
      bodies(i)%pressure = (0.66666666666666667) * bodies(i)%internal_energy * bodies(i)%density
    end do

    call navigate_tree(root, bodies, 0.5_dp)
    call get_SPH(root, bodies)

    do i = 1, root%n_particles
      bodies(i)%velocity = bodies(i)%velocity + (bodies(i)%acceleration)*dt/2
      bodies(i)%internal_energy = bodies(i)%internal_energy + (bodies(i)%internal_energy_rate)*dt/2
    end do

    do i = 1, root%n_particles
      bodies(i)%position = bodies(i)%position + (bodies(i)%velocity)*dt/2
    end do
    deallocate(root)

    allocate(root)
    ! Initialize root node
    root%center = [150.0_dp, 150.0_dp, 150.0_dp]
    root%size = 10000.0_dp
    root%n_particles = size(bodies)
    allocate(root%particles(size(bodies)))
    root%particles = bodies

    call build_tree(root, 1000, 1)
    call get_density(root, bodies)

    do i = 1, root%n_particles
      bodies(i)%pressure = (0.66666666666666667) * bodies(i)%internal_energy * bodies(i)%density
    end do

    call navigate_tree(root, bodies, 0.5_dp)
    call get_SPH(root, bodies)

    do i = 1, root%n_particles
      bodies(i)%velocity = bodies(i)%velocity + (bodies(i)%acceleration)*dt/2
      bodies(i)%internal_energy = bodies(i)%internal_energy + (bodies(i)%internal_energy_rate)*dt/2

      bodies(i)%position = bodies(i)%position + (bodies(i)%velocity)*dt/2
    end do
    deallocate(root)
    t = t + dt
  end do
end subroutine
end module octree_module

program barnes_hut
  use octree_module
  implicit none
  integer :: n, total_particles,i
  type(particle), allocatable :: bodies(:)
  type(branch), allocatable :: root

  call init_kernel_table()

  !make some random points
  total_particles = 2000
  allocate(bodies(total_particles))

  do n = 1, total_particles
    call random_number(bodies(n)%position) ! Generate random positions
    bodies(n)%position = 1000.0_dp * bodies(n)%position ! Scale to a range
    bodies(n)%mass = 0.1_dp
    call random_number(bodies(n)%velocity)
    bodies(n)%velocity = 0.0_dp * bodies(n)%velocity
    bodies(n)%pressure = 1.0_dp
    bodies(n)%internal_energy = 1.0_dp
    bodies(n)%density = 0.0_dp
  end do


  allocate(root)
  ! Initialize root node
  root%center = [150.0_dp, 150.0_dp, 150.0_dp]
  root%size = 10000.0_dp
  root%n_particles = size(bodies)
  allocate(root%particles(size(bodies)))
  root%particles = bodies


  call build_tree(root, 1000, 1)
  call get_density(root, bodies)
  do i = 1, root%n_particles
      bodies(i)%pressure = (0.66666666666666667) * bodies(i)%internal_energy * bodies(i)%density
  end do
  call navigate_tree(root, bodies, 0.5_dp)
  call get_SPH(root, bodies)
  deallocate(root)


  call simulate(bodies)
  print *, bodies(1)%position(1),bodies(1)%position(2),bodies(1)%position(3)
end program barnes_hut




