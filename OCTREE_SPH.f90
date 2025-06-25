module octree_module
  implicit none
  real, parameter :: G = 6.67430e-11
  real, parameter :: softening = 1.0e-5
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nq = 1000  ! Number of samples between q = 0 and q = 2
  real(dp), allocatable :: w_table(:), dw_table(:)
  real(dp), parameter :: dq = 2.0_dp / nq


  !types here act like dictionaries with "sets" of arrays
  type :: particle
    real(dp) :: mass, density, internal_energy, pressure, smoothing
    real(dp) :: density_rate, internal_energy_rate, pressure_rate
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

    !clean up the parent node - THIS BREAKS EVERYTHING IN NEIGHBOUR SEARCH
    !deallocate(node%particles)

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

! Navigates the tree to get to the nodes within 2h of each particle
recursive subroutine kernel_thing(node, body)
  implicit none
  type(branch), intent(in) :: node
  type(particle), intent(inout) :: body(:)
  real(dp) :: com_dist, dr, Wj, dWj_mag, mj, vdotW
  real(dp), dimension(3) :: nr, dWj, vij
  integer :: i, j, k
  logical :: has_children

  has_children = allocated(node%children)

  do i = 1, size(body)
    com_dist = sqrt(sum((body(i)%position - node%mass_center)**2))
    if (node%n_particles > 0 .and. com_dist < 0.70710678118_dp * node%size .and. com_dist < 2 * body(i)%smoothing) then
      do j = 1, node%n_particles
        !SPH forces/ kernel calc
        nr = (body(i)%position - node%particles(j)%position)
        vij = (body(i)%velocity - node%particles(j)%velocity)
        dr = sqrt(sum(nr**2)) !might need to switch i and j
        if (dr == 0.0_dp) cycle
        nr = nr / dr
        call lookup_kernel(dr, body(i)%smoothing,  Wj, dWj_mag)
        Wj = Wj / (3.14159265359_dp * body(i)%smoothing**3) 
        dWj = nr * dWj_mag / (3.14159265359_dp * body(i)%smoothing**4)

        mj = node%particles(j)%mass
        !now do forces and shit
        !pressure force
        body(i)%acceleration = body(i)%acceleration - mj * ((body(i)%pressure/(body(i)%density*body(i)%density)) + (node%particles(j)%pressure / (node%particles(j)%density*node%particles(j)%density))) * dWj

        !drho/dt
        vdotW = sum(vij * dWj)
        body(i)%density_rate = body(i)%density_rate + mj * vdotW

        !rate change in internal energy (might be faster to take out p/rho term and do as a bulk calc after)
        body(i)%internal_energy_rate = body(i)%internal_energy_rate + (body(i)%pressure/body(i)%density)*mj*(vdotW)
      end do

    else if (has_children .and. com_dist < 2 * body(i)%smoothing) then
      do k = 1, size(node%children)
        call kernel_thing(node%children(k), body(i:i))
      end do
    end if
  end do

end subroutine kernel_thing

subroutine simulate()
  implicit none

end subroutine simulate  
end module octree_module

program barnes_hut
  use octree_module
  implicit none
  integer :: n, total_particles, a
  real(dp) :: time_i, time_f, Wi, dWi
  type(particle), allocatable :: bodies(:)
  type(branch), allocatable :: root

  call init_kernel_table()

  !make some random points
  total_particles = 1000
  allocate(bodies(total_particles))

  do n = 1, total_particles
    call random_number(bodies(n)%position) ! Generate random positions
    bodies(n)%position = 300.0_dp * bodies(n)%position ! Scale to a range
    bodies(n)%mass = 10.0_dp
    call random_number(bodies(n)%velocity)
    bodies(n)%velocity = 30.0_dp * bodies(n)%velocity
    bodies(n)%smoothing = 10.0_dp
    bodies(n)%pressure = 100.0_dp
    bodies(n)%internal_energy = 1.0_dp
    bodies(n)%density = 0.0_dp
  end do
  do n = 1, total_particles
    do a = 1, total_particles
      call lookup_kernel(norm2(bodies(n)%position - bodies(a)%position), bodies(n)%smoothing, Wi, dWi)
      bodies(n)%density = bodies(n)%density + bodies(a)%mass * Wi / (3.14159265359_dp * bodies(n)%smoothing**3)
    end do
  end do

  !testing the time to perform multiple loops

  call cpu_time(time_i)
  do a = 1, 10
    allocate(root)
    ! Initialize root node
    root%center = [150.0_dp, 150.0_dp, 150.0_dp]
    root%size = 300.0_dp
    root%n_particles = size(bodies)
    allocate(root%particles(size(bodies)))
    root%particles = bodies


    call build_tree(root, 1000, 1)

    
    call navigate_tree(root, bodies, 0.5_dp)
    call kernel_thing(root, bodies)
    deallocate(root)
  end do
  call cpu_time(time_f)
  print *, 'program finished', time_f - time_i
  print *, '1day real to', (86400 * 10 / (time_f - time_i)) / 365.25 , ' years in sim'

  print *, bodies%acceleration(1)
end program barnes_hut





