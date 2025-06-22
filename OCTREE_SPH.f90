module octree_module
  implicit none
  real, parameter :: G = 6.67430e-11
  real, parameter :: softening = 1.0e-5
  integer, parameter :: dp = kind(1.0d0)

  !types here act like dictionaries with "sets" of arrays
  type :: particle
    real(dp) :: mass, density, internal_energy, pressure, smoothing
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

  type :: kkgrad !information on the kernel and grad(kernel)
    real(dp), allocatable :: k(:)
    real(dp), allocatable, dimension(3) :: kgrad(:,:)
  end type kkgrad

  type :: settings
    real(dp) :: dt, scale_dt, size
  end type settings

  contains

  ! All this does is check an volume for particles, then subdivide into 8
  ! and allocate particles to the right octant, nothing more 
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

  ! Used by the kernel_thing subroutine when actually smoothing all particles
subroutine cubic_spline(body, node, w)
  implicit none
  type(particle), intent(in) :: body
  type(branch), intent(in) :: node
  type(kkgrad), intent(inout) :: w

  integer :: i, n
  real(dp) :: r, h, q, norm
  real(dp), dimension(3) :: dr

  n = node%n_particles
  if (n == 0) return

  ! Allocate memory for kernel weights and gradients
  allocate(w%k(n))
  allocate(w%kgrad(3, n))

  h = body%smoothing
  norm = 1.0_dp / (3.14159265359_dp * h**3)

  do i = 1, n
    dr = body%position - node%particles(i)%position
    r = sqrt(sum(dr**2))
    if (r == 0.0_dp) then
      q = 0.0_dp
    else
      q = r / h
    end if

    !smoothing kernel and its gradient
    if (q >= 0.0_dp .and. q <= 1.0_dp) then
      w%k(i) = norm * (1.0_dp - 1.5_dp*q**2 + 0.75_dp*q**3)
      w%kgrad(:,i) = (dr/r) * norm * (-3.0_dp*q/h + 2.25_dp*q**2/h)
    else if (q > 1.0_dp .and. q <= 2.0_dp) then
      w%k(i) = norm * 0.25_dp * (2.0_dp - q)**3
      w%kgrad(:,i) = (dr/r) * norm * (-0.75_dp * (2.0_dp - q)**2 / h)
    else
      w%k(i) = 0.0_dp
      w%kgrad(:,i) = 0.0_dp
    end if
  end do

end subroutine cubic_spline

! Navigates the tree to get to the nodes within 2h of each particle
recursive subroutine kernel_thing(node, body)
  implicit none
  type(branch), intent(in) :: node
  type(particle), intent(inout) :: body(:)
  real(dp) :: com_dist
  integer :: i, k
  type(kkgrad) :: w
  logical :: has_children

  has_children = allocated(node%children)

  do i = 1, size(body)
    com_dist = sqrt(sum((body(i)%position - node%mass_center)**2))

    if (node%n_particles > 0 .and. com_dist < 0.70710678118_dp * node%size .and. com_dist < 4 * body(i)%smoothing) then
      ! Evaluate SPH forces
      call cubic_spline(body(i), node, w)
      !print *, w%k(:)
      ! Deallocate after use
      if (allocated(w%k)) deallocate(w%k)
      if (allocated(w%kgrad)) deallocate(w%kgrad)

      return
    else if (has_children .and. com_dist < 2 * body(i)%smoothing) then
      ! Recurse through children
      do k = 1, size(node%children)
        call kernel_thing(node%children(k), body(i:i))
      end do
    end if
  end do

end subroutine kernel_thing


end module octree_module

program barnes_hut
  use octree_module
  implicit none
  integer :: n, total_particles, a
  real :: time_i, time_f
  type(particle), allocatable :: bodies(:)
  type(branch), allocatable :: root


  !make some random points
  total_particles = 5000
  allocate(bodies(total_particles))

  do n = 1, total_particles
    call random_number(bodies(n)%position) ! Generate random positions
    bodies(n)%position = 300.0_dp * bodies(n)%position ! Scale to a range
    bodies(n)%mass = 100.0_dp
    bodies(n)%velocity = [0.0_dp, 0.0_dp, 0.0_dp]
    bodies(n)%smoothing = 140.0_dp
    bodies(n)%pressure = 1.0_dp
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


    call build_tree(root, 100, 1)

    
    call navigate_tree(root, bodies, 0.5_dp)
    call kernel_thing(root, bodies)
    !print *, root%n_particles
    deallocate(root)
  end do
  call cpu_time(time_f)
  print *, 'program finished', time_f - time_i
  print *, '30day real to', (86400 * 10 * 30 / (time_f - time_i)) / 365.25 , ' years in sim'

  
end program barnes_hut


