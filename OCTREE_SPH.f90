module octree_module
  implicit none
  ! Global Parameters:
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: G = 39.47841760435743 !AU^3/(Msun*yr^2) Gravitational constant (in SI units)
  integer, parameter :: nq = 1000, max_depth = 10000          ! Number of samples for the SPH kernel lookup tables.
                                           ! This determines the resolution of the pre-computed kernel values.
  real(dp), allocatable :: w_table(:), dw_table(:), grav_table(:) ! Allocatable arrays to store pre-computed SPH kernel (W)
                                                   ! and its derivative (dW/dr) values.
  real(dp), parameter :: dq = 2.0_dp / nq 
  real(dp), parameter :: smoothing = 9.5e-5_dp, bounding_size = 100000.0_dp

  ! Represents a single SPH particle with its physical properties.
  type :: particle
    integer :: number !this is just an identifier for when we sync the tree with the particles
    real(dp) :: mass
    real(dp) :: density
    real(dp) :: internal_energy
    real(dp) :: pressure
    real(dp) :: sound_speed
    real(dp) :: internal_energy_rate  
    real(dp), dimension(3) :: position 
    real(dp), dimension(3) :: velocity 
    real(dp), dimension(3) :: acceleration 
  end type particle

  !sink particle type
  type :: sink
    real(dp) :: mass
    real(dp) :: radius
    real(dp), dimension(3) :: spin
    real(dp), dimension(3) :: position 
    real(dp), dimension(3) :: velocity 
    real(dp), dimension(3) :: acceleration 
  end type sink  

  ! Represents a node (branch) in the Barnes-Hut octree.
  type :: branch
    real(dp), dimension(3) :: center
    real(dp) :: size
    integer :: n_particles  
    type(Particle), allocatable :: particles(:) 
    type(sink), allocatable :: sinks(:)
    real(dp) :: mass_total               
    real(dp), dimension(3) :: mass_center 
    type(branch), allocatable :: children(:) ! Array of 8 child branches

  end type branch

  contains

!---------------------------Kernel Tables and Lookup routines-----------------------------
  ! Initializes the global lookup tables (w_table and dw_table) for the SPH kernel
  ! and its derivative.
  subroutine init_kernel_table()
    integer :: i
    real(dp) :: q      ! Normalized distance (r/h), ranging from 0 to 2

    ! Allocate memory for the kernel and derivative tables
    allocate(w_table(0:nq))
    allocate(dw_table(0:nq))

    do i = 0, nq
      q = i * dq ! Calculate current normalized distance
      if (q >= 0.0_dp .and. q <= 1.0_dp) then
        w_table(i) = 1.0_dp - 1.5_dp*q**2 + 0.75_dp*q**3
        dw_table(i) = -3.0_dp*q + 2.25_dp*q**2

      else if (q > 1.0_dp .and. q <= 2.0_dp) then
        w_table(i) = 0.25_dp * (2.0_dp - q)**3
        dw_table(i) = -0.75_dp * (2.0_dp - q)**2

      else
        !Outside the compact support (r/h > 2), kernel and derivative are zero
        w_table(i) = 0.0_dp
        dw_table(i) = 0.0_dp
      end if
    end do
  end subroutine init_kernel_table

  subroutine init_grav_kernel_table()
    integer :: i
    real(dp) :: q      ! Normalized distance (r/h), ranging from 0 to 2

    ! Allocate memory for the kernel and derivative tables
    allocate(grav_table(0:nq))

    do i = 0, nq
      q = i * dq ! Calculate current normalized distance
      if (q >= 0.0_dp .and. q <= 1.0_dp) then
        grav_table(i) = ((40.0_dp*q**3) - (36.0_dp * q**5) + (15.0_dp * q**6))/30.0_dp

      else if (q > 1.0_dp .and. q <= 2.0_dp) then
        grav_table(i) = ((80.0_dp * q**3) - (90.0_dp * q**4) + (36.0_dp * q**5) - (5*q**6) - 2)/30.0_dp

      else
        !Outside the compact support (r/h > 2), kernel and derivative are zero
        grav_table(i) = 1.0_dp
      end if
    end do
  end subroutine init_grav_kernel_table

  ! Retrieves interpolated SPH kernel values (W) and their derivatives (dW/dr)
  ! from the pre-computed tables for a given distance 'r' and smoothing length 'hi'.
  subroutine lookup_kernel(r, hi, Wi, dWi)
    real(dp), intent(in) :: r, hi 
    real(dp), intent(out) :: Wi, dWi 
    integer :: i
    real(dp) :: alpha, qi

    qi = r / hi
    ! Check if the normalized distance is within the kernel's compact support [0, 2]
    if (qi >= 0.0_dp .and. qi <= 2.0_dp) then
      i = int(qi / dq) 
      alpha = (qi - i * dq) / dq 
      ! Perform linear interpolation to get the kernel and derivative values
      Wi = (1.0_dp - alpha) * w_table(i) + alpha * w_table(i+1)
      dWi = (1.0_dp - alpha) * dw_table(i) + alpha * dw_table(i+1)
    else
      ! If outside the compact support, kernel values are zero
      Wi = 0.0_dp
      dWi = 0.0_dp
    end if
  end subroutine lookup_kernel

  subroutine lookup_grav_kernel(r, hi, Wi)
    real(dp), intent(in) :: r, hi 
    real(dp), intent(out) :: Wi
    integer :: i
    real(dp) :: alpha, qi

    qi = r / hi
    ! Check if the normalized distance is within the kernel's compact support [0, 2]
    if (qi >= 0.0_dp .and. qi <= 2.0_dp) then
      i = int(qi / dq) 
      alpha = (qi - i * dq) / dq 
      ! Perform linear interpolation to get the kernel and derivative values
      Wi = (1.0_dp - alpha) * grav_table(i) + alpha * grav_table(i+1)
    else
      ! If outside the compact support, kernel values are 1
      Wi = 1.0_dp
    end if
  end subroutine lookup_grav_kernel

!--------------------------Tree building and grav interaction------------------------------
  recursive subroutine build_tree(node, depth, max_particles)
    implicit none
    type(branch), intent(inout) :: node
    integer, intent(in) :: depth, max_particles ! Input: current recursion depth, max particles allowed in a leaf node.
    integer :: i, j, n                        ! Loop indices and temporary size variable
    real(dp), dimension(3) :: offset         ! Spatial offset for calculating child node centers
    type(Particle), allocatable :: child_particles(:)
    integer :: child_index
    integer :: p

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

    ! Base case for recursion:
    ! If the number of particles in the node is less than or equal to `max_particles` (making it a leaf node),
    ! OR if the maximum recursion `depth` has been reached, stop subdividing.
    if (size(node%particles) <= max_particles .or. depth == 0) return

    ! If not a base case, allocate 8 child nodes for this branch.
    allocate(node%children(8))
    do i = 1, 8
      node%children(i)%size = node%size / 2.0_dp
      node%children(i)%n_particles = 0           !Initialize particle count for child to zero
      allocate(node%children(i)%particles(0))    !Allocate an empty array for particles in the child node
      node%children(i)%center = node%center      !Start with parent's center for calculating child's center

      do j = 1, 3
        if (btest(i - 1, j - 1)) then ! If j-th bit is set (for x, y, or z)
          offset(j) = 0.25_dp * node%size ! Positive offset for this dimension
        else
          offset(j) = -0.25_dp * node%size ! Negative offset for this dimension
        end if
      end do
      node%children(i)%center = node%center + offset ! Set the child's actual center
    end do

    ! Distribute particles from the current node to its newly created children.
    do p = 1, size(node%particles)
      child_index = 1 ! Initialize child index (corresponds to octant 1)
      ! Determine which child octant the current particle `node%particles(p)` belongs to.
      do j = 1, 3
        if (node%particles(p)%position(j) > node%center(j)) then
          child_index = ibset(child_index - 1, j - 1) + 1
        end if
      end do

      ! Add the particle to the selected child node's particle list.
      n = size(node%children(child_index)%particles) ! Current number of particles in the target child
      ! `move_alloc` is used for efficient reallocation of arrays, preventing full data copy.
      call move_alloc(node%children(child_index)%particles, child_particles)
      allocate(node%children(child_index)%particles(n + 1)) ! Allocate space for one more particle
      if (n > 0) node%children(child_index)%particles(1:n) = child_particles ! Copy existing particles back
      node%children(child_index)%particles(n + 1) = node%particles(p) ! Add the new particle
      node%children(child_index)%n_particles = node%children(child_index)%n_particles + 1 ! Increment child's particle count
    end do

    ! Recursively call build_tree for each child node that contains particles.
    do i = 1, 8
      if (node%children(i)%n_particles > 0) then
        call build_tree(node%children(i), depth - 1, max_particles)
      end if
    end do
  end subroutine build_tree

  ! Recursive Subroutine: navigate_tree
  ! Traverses the octree to calculate the gravitational acceleration on a given 'body' particle
  ! using the Barnes-Hut approximation.
  recursive subroutine particle_gravforces(node, body, theta)
  implicit none
  type(branch), intent(in) :: node
  type(particle), intent(inout) :: body(:)
  real(dp), intent(in) :: theta
  real(dp) :: dist, d2, W
  real(dp), dimension(3) :: direction
  integer :: i, k

  ! Loop over the (single) particle provided in the 'body' array slice.
  do i = 1, size(body)
    direction = body(i)%position - node%mass_center ! Vector pointing from particle to the node's center of mass
    d2 = sum(direction**2) + smoothing**2 
    dist = sqrt(d2)  

    ! Barnes-Hut criterion:
    if ((node%size / dist) < theta .or. .not. allocated(node%children)) then
      ! Calculate gravitational acceleration if the node has mass and distance is positive.
      if (node%mass_total > 0.0_dp .and. dist > 0.0_dp) then
        call lookup_grav_kernel(dist, smoothing, W)
        body(i)%acceleration = body(i)%acceleration - (G * node%mass_total * W *direction / (dist**3))
      end if
    else if (allocated(node%children)) then
      ! If the criterion is not met and the node has children, recurse into each child node.
      do k = 1, size(node%children)
        if (node%children(k)%n_particles > 0) then
          call particle_gravforces(node%children(k), body(i:i), theta) ! Recurse, passing the same particle
        end if
      end do
    end if
  end do
  end subroutine particle_gravforces

!--------------------------SPH and Density finding subroutines---------------------------------

  ! Orchestrates the calculation of SPH forces and internal energy rates for all particles.
  subroutine get_SPH(root, body)
    implicit none
    type(branch), intent(in) :: root           ! The root node of the octree.
    type(particle), intent(inout) :: body(:)  ! Array of all particles in the simulation.

    integer :: i

    ! Loop through each particle in the main `body` array.
    do i = 1, size(body)
      call SPH_tree_search(root, body(i))
    end do
  end subroutine get_SPH

  ! Searches the octree for neighbors of a given 'body' particle within its smoothing length,
  ! and calculates the SPH interaction terms (pressure force and internal energy rate).
  recursive subroutine SPH_tree_search(node, body)
    implicit none
    type(branch), intent(in) :: node          ! Current octree node being examined.
    type(particle), intent(inout) :: body     ! The single particle for which SPH interactions are calculated.
    real(dp), dimension(3) :: over_dr, nr, dWj, vij ! over_dr: the 'overlap distance'
    real(dp) :: Wj, dr, mj, dWj_mag, vdotW, vdotr, vis_nu,viscous_cont, avg_sound_speed
    integer :: j
    logical :: has_children

    has_children = allocated(node%children)
    ! Calculate vector from the 'body' particle's position to the current node's center.
    over_dr = (body%position - node%center)

    ! Recursion condition:
    ! If the node contains multiple particles AND the 'body' particle's smoothing sphere
    ! potentially overlaps this node's bounding box (indicated by `abs(over_dr) < (2*smoothing + node%size/2)`),
    ! AND the node has children, then recurse into its children.
    if (node%n_particles > 1 .and. all(abs(over_dr) < (2.0_dp*smoothing + node%size/2.0_dp)) .and. has_children) then
      do j = 1, size(node%children)
        ! Only recurse if the child node actually contains particles.
        if (node%children(j)%n_particles > 0) then
          call SPH_tree_search(node%children(j), body)
        end if
      end do

    ! Base case:
    ! If the node is a leaf node (contains exactly one particle) AND the 'body' particle's smoothing sphere
    ! potentially overlaps this node's bounding box, then this single particle in the leaf node
    ! is a candidate neighbor for 'body'.
    else if (node%n_particles == 1 .and. all(abs(over_dr) < (2.0_dp*smoothing + node%size/2.0_dp))) then
      ! This is where the interaction with a single neighbor particle (node%particles(1)) occurs.

      if (node%particles(1)%number == body%number) return
      
      nr = (body%position - node%particles(1)%position) 
      dr = sqrt(sum(nr**2))
      if (dr <= 1.0e-10) return
      vij = (body%velocity - node%particles(1)%velocity)
      vdotr = sum(vij * nr)

      if (vdotr >= 0) vdotr = 0.0_dp

      nr = nr / dr ! Normalize the separation vector

      ! Lookup kernel (W) and derivative (dW/dr) values for the calculated distance 'dr'.
      call lookup_kernel(dr, smoothing, Wj, dWj_mag)
      ! Normalize kernel and its gradient for 3D cubic spline:
      Wj = Wj / (3.14159265359_dp * smoothing**3)
      dWj = nr * dWj_mag / (3.14159265359_dp * smoothing**4)
      mj = node%particles(1)%mass ! Mass of the neighbor particle

      !get viscous contributions
      vis_nu = (smoothing * vdotr)/(dr*dr + 0.1*smoothing*smoothing)
      avg_sound_speed = - 0.25 * (body%sound_speed + node%particles(1)%sound_speed) !average sound speed multiplied by - 1/2
      viscous_cont = (avg_sound_speed * vis_nu + vis_nu*vis_nu) / (body%density + node%particles(1)%density)

      ! Accumulate SPH pressure force:
      body%acceleration = body%acceleration - mj * ((body%pressure/(body%density * body%density)) + &
                                                  (node%particles(1)%pressure / (node%particles(1)%density * node%particles(1)%density)) + viscous_cont) * dWj

      vdotW = sum(vij * dWj)
      ! Accumulate rate of change in internal energy:
      body%internal_energy_rate = body%internal_energy_rate + ((body%pressure/body%density)*mj*(vdotW)) + (0.5 * viscous_cont * mj* vdotW)
      return 
    end if
    ! If neither recursion nor leaf node interaction conditions are met, simply return.
  end subroutine SPH_tree_search

  ! Subroutine: get_density
  subroutine get_density(root, body)
    implicit none
    type(branch), intent(inout) :: root
    type(particle), intent(inout) :: body(:)    ! Array of all particles.

    integer :: i                                  ! Loop index
    logical :: has_children                       ! Local flag for node having children

    has_children = allocated(root%children)
    ! Loop through each particle to calculate its density.
    do i = 1, size(body)
      body(i)%density = 0.0_dp ! Initialize density to zero before accumulating contributions for particle `i`.
      ! Search the tree to find all neighbors and accumulate density contributions for `body(i)`.
      call density_tree_search(root, body(i))
    end do

    !now distribute the density to the tree
    call sync_density_to_tree(root, body)

  end subroutine get_density

  ! Recursively traverses the octree to find neighbors for a given 'body' particle
  recursive subroutine density_tree_search(node, body)
    implicit none
    type(branch), intent(inout) :: node 
    type(particle), intent(inout) :: body 
    real(dp), dimension(3) :: over_dr, nr
    real(dp) :: Wj, dr, mj, dWj_mag
    integer :: j
    logical :: has_children

    has_children = allocated(node%children)
    ! Calculate vector from the 'body' particle's position to the current node's center.
    over_dr = (body%position - node%center)

    ! Recursion condition:
    ! potentially overlaps this node's bounding box, AND the node has children, then recurse.
    if (node%n_particles > 1 .and. all(abs(over_dr) < (2.0_dp*smoothing + node%size/2.0_dp)) .and. has_children) then
      do j = 1, size(node%children)
        ! Only recurse if the child node actually contains particles.
        if (node%children(j)%n_particles > 0) then
          call density_tree_search(node%children(j), body)
        end if
      end do

    ! Base case:
    ! If the node is a leaf node (contains exactly one particle) AND the 'body' particle's smoothing sphere
    ! potentially overlaps this node's bounding box, then this single particle in the leaf node
    ! is a candidate neighbor for 'body'.
    else if (node%n_particles == 1 .and. all(abs(over_dr) < (2.0_dp*smoothing + node%size/2.0_dp))) then
      ! Get the neighbor particle from the leaf node (node%particles(1)).
      nr = (body%position - node%particles(1)%position) 
      dr = sqrt(sum(nr**2))

      ! Lookup kernel (W) value for the calculated distance 'dr'.
      call lookup_kernel(dr, smoothing, Wj, dWj_mag) ! dWj_mag is not used for density, but still returned.
      ! Normalize kernel for 3D cubic spline.
      Wj = Wj / (3.14159265359_dp * smoothing**3)
      mj = node%particles(1)%mass ! Mass of the neighbor particle

      ! Accumulate density for the 'body' particle: rho_i = sum_j m_j * W_ij
      body%density = body%density + mj * Wj
      return
    end if
  end subroutine density_tree_search

  recursive subroutine sync_density_to_tree(node, bodies) !the tree needs to be synced up to the calculated densities in bodies
    implicit none
    type(branch), intent(inout) :: node
    type(particle), intent(inout)  :: bodies(:)
    integer :: i, num

    ! For each particle in this node, find the corresponding particle in bodies by position and set density
    do i = 1, size(node%particles)
      num = node%particles(i)%number
      
      node%particles(i)%density = bodies(num)%density
      if (node%particles(i)%density == 0.0_dp) stop
      
      node%particles(i)%pressure = (0.66666666666666667_dp) * bodies(num)%internal_energy * bodies(num)%density
      bodies(num)%pressure = node%particles(i)%pressure

      node%particles(i)%sound_speed = sqrt(1.66666666666666667_dp*bodies(num)%pressure/bodies(num)%density)
      bodies(num)%sound_speed = node%particles(i)%sound_speed
    end do

    ! Recurse into children if present
    if (allocated(node%children)) then
      do i = 1, size(node%children)
        if (node%children(i)%n_particles > 0) then
          call sync_density_to_tree(node%children(i), bodies)
        end if
      end do
    end if
  end subroutine sync_density_to_tree
  
!------------------Sink and particle deletion subroutines-----------------
  subroutine check_bounds(bodies)
    implicit none
    type(particle), allocatable, intent(inout) :: bodies(:)
    logical :: keep_mask(size(bodies))
    integer :: i

    ! Create logical mask: .TRUE. for bodies inside bounding box
    keep_mask = [(all(abs(bodies(i)%position) <= bounding_size), i = 1, size(bodies))]

    ! Filter and assign
    bodies = pack(bodies, keep_mask)
    ! reassign numbers
    do i = 1, size(bodies)
        bodies(i)%number = i
    end do
  end subroutine check_bounds

  !subroutine check4sinkcreate(bodies)
  !  implicit none 
  !  type(particle), intent(inout) :: bodies(:)
  !  type(sink), allocatable :: sinks
  !
  !
  !end subroutine check4sinkcreate

  subroutine initiate_sink_accretion(sinks, bodies, root)
    implicit none
    type(sink), intent(inout) :: sinks(:)
    type(particle), intent(inout), allocatable :: bodies(:)
    type(branch), intent(in) :: root
    logical :: keep_mask(size(sinks)+1, size(bodies))
    integer :: i, j

    !###########IDENTIFY ACCRETABLE PARTICLES##############
    do i = 1, size(sinks)
      keep_mask(i,:) = [((i == i), j = 1, size(bodies))] !sets all to .true. by default as to NOT accrete
      call sink2gasdists(sinks(i), root, keep_mask(i,:))

      sinks(i)%mass = sinks(i)%mass + sum(pack(bodies%mass, .not. keep_mask(i,:)))
      sinks(i)%position = (sinks(i)%mass * sinks(i)%position + [ &
        sum(pack(bodies%mass * bodies%position(1), .not. (keep_mask(i,:)))), &
        sum(pack(bodies%mass * bodies%position(2), .not. (keep_mask(i,:)))), &
        sum(pack(bodies%mass * bodies%position(3), .not. (keep_mask(i,:)))) ]) / (sinks(i)%mass)
      
      !also need somthing to track the angular momentum

    end do

    !###########PACK AND UPDATE SINK####################
    ! Create logical mask: .TRUE. for bodies inside bounding box
    keep_mask(size(sinks)+1,:) = [(all(abs(bodies(j)%position) <= bounding_size), j = 1, size(bodies))] ! .true. if inside bounds, .false. if outside

    call pack_sinks(bodies, keep_mask)
  end subroutine initiate_sink_accretion

  recursive subroutine sink2gasdists(sink_i, node, mask)
    implicit none
    type(sink), intent(inout) :: sink_i
    type(branch), intent(in) :: node
    logical, intent(inout) :: mask(:)
    real(dp) :: dr
    real(dp), dimension(3) :: over_dr
    integer :: j

    over_dr = node%center - sink_i%position

    !recurse if not deep enough
    if (node%n_particles > 1 .and. all(abs(over_dr) < (sink_i%radius + node%size/2.0_dp)) .and. allocated(node%children)) then
      do j = 1, size(node%children)
        if (node%n_particles > 0) call sink2gasdists(sink_i, node%children(j), mask)
      end do
      return
    
    !base case of not needing to recurse
    else if (node%n_particles == 1 .and. all(abs(over_dr) < (2*sink_i%radius + node%size/2.0_dp))) then
      dr = sum(sqrt(node%center*node%center - sink_i%position * sink_i%position))
      if (dr < sink_i%radius) then
        mask(node%particles%number) = .false.
      end if
      return
    end if
    return
  end subroutine sink2gasdists

  subroutine pack_sinks(bodies, mask)
    implicit none
    type(particle), intent(inout), allocatable :: bodies(:)
    logical, intent(in) :: mask(:,:)
    logical :: d1mask(size(mask(1,:)))
    integer :: i

    !pack bodies with vertically or'd mask
    d1mask = .not.any(.not.mask, dim=1)
    bodies = pack(bodies, d1mask)
    !print *, d1mask

    do i = 1, size(bodies)
      bodies%number = i
    end do
  end subroutine pack_sinks

  !direct sink forces calculation
  subroutine sink_gravforces(bodies, sinks)
    implicit none
    type(particle), intent(inout) :: bodies(:)
    type(sink), intent(inout) :: sinks(:)
    real(dp) :: dr
    real(dp), dimension(3) :: vect_dr, dist_weighting
    integer :: i, j

    do i = 1, size(sinks)
      do j = 1, size(bodies)
        vect_dr = bodies(j)%position - sinks(i)%position
        dr = sqrt(sum((vect_dr**2)))

        dist_weighting = G * vect_dr / (dr*dr*dr)
        !print *, 'weit', vect_dr
        sinks(i)%acceleration = [0.0, 0.0,0.0 ]!sinks(i)%acceleration + (bodies(j)%mass * dist_weighting)
        bodies(j)%acceleration = bodies(j)%acceleration - (sinks(i)%mass * dist_weighting)
      end do
    end do
  end subroutine sink_gravforces

!------------------File Reading and writing---------------------------------
  subroutine read_data_from_file(filename, bodies, max_lines)
    implicit none

    ! Declare input argument
    character(len=*), intent(in) :: filename

    ! Define maximum number of lines and parameters
    integer :: max_lines! Adjust this based on your file size
    integer :: status, num_lines, i
    character(len=256) :: header_line
    real(kind=8), allocatable, dimension(:) :: x, y, z, vx, vy, vz, energy, mass
    type(particle), allocatable, intent(inout) :: bodies(:)

    ! Allocate arrays
    ALLOCATE(x(MAX_LINES), y(MAX_LINES), z(MAX_LINES), &
             vx(MAX_LINES), vy(MAX_LINES), vz(MAX_LINES), &
             energy(MAX_LINES), mass(MAX_LINES), STAT=status)
    IF (status /= 0) THEN
        WRITE(*,*) 'Error allocating memory for arrays.'
        RETURN ! Exit the subroutine if allocation fails
    END IF

    ! Open the file
    OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)

    IF (status /= 0) THEN
        WRITE(*,*) 'Error opening file: ', TRIM(filename)
        WRITE(*,*) 'Check if the file exists and has correct permissions.'
        DEALLOCATE(x, y, z, vx, vy, vz, energy, mass) ! Deallocate on error
        RETURN ! Exit the subroutine on file open error
    END IF

    ! Read the header line (and discard it)
    READ(10, '(A)') header_line

    ! Loop to read data
    num_lines = 0
    DO
        READ(10, *, IOSTAT=status) x(num_lines+1), y(num_lines+1), z(num_lines+1), &
                                   vx(num_lines+1), vy(num_lines+1), vz(num_lines+1), &
                                   energy(num_lines+1), mass(num_lines+1)

        IF (status /= 0) EXIT ! Exit loop on end-of-file or error

        num_lines = num_lines + 1

        IF (num_lines > MAX_LINES) THEN
            WRITE(*,*) 'Warning: Maximum number of lines (', MAX_LINES, ') reached for file ', TRIM(filename), '.'
            WRITE(*,*) 'Some data might not be read. Consider increasing MAX_LINES.'
            EXIT
        END IF
    END DO

    ! Close the file
    CLOSE(UNIT=10)

    ! Print a confirmation message and some data (optional)
    write(*,*) 'Successfully read ', num_lines, ' lines of data from ', TRIM(filename), '.'
    if (num_lines > 0) then
    else
        write(*,*) 'No data points were read from ', TRIM(filename), '.'
    end if

    allocate(bodies(size(x)))

    bodies%position(1) = x
    bodies%position(2) = y
    bodies%position(3) = z
    bodies%velocity(1) = vx
    bodies%velocity(2) = vy
    bodies%velocity(3) = vz
    bodies%internal_energy = energy
    bodies%mass = mass
    DEALLOCATE(x, y, z, vx, vy, vz, energy, mass)

    do i =1, size(bodies)
      bodies(i)%number = i
    end do 

  END SUBROUTINE read_data_from_file

  subroutine make_save(bodies, number)
    implicit none
    integer :: io, i, number
    type(particle) :: bodies(:)
    character(len = 256) :: savename

    write(savename, '(A,I0,A)') 'save', number, '.txt'

    open(newunit=io, file=savename, status="new", action = "write")
    write(io, *) 'x  ','y  ','z  ','vx  ','vy ','vz ','energy ','density  ','mass  '
    do i = 1, size(bodies)
      write(io, *)  bodies(i)%position(1),bodies(i)%position(2),bodies(i)%position(3),bodies(i)%velocity(1),bodies(i)%velocity(2),bodies(i)%velocity(3), bodies(i)%internal_energy, bodies(i)%density, bodies(i)%mass
    end do  
    close(io)

  end subroutine

!----------------------- kick drift kick------------------------
  !velocity update
  subroutine kick(bodies, sinks, dt)
    implicit none
    type(particle), intent(inout) :: bodies(:)
    type(sink), intent(inout) :: sinks(:)
    real(dp), intent(in) :: dt

    !v(i+0.5)
    bodies%velocity(1) = bodies%velocity(1) + 0.5_dp*bodies%acceleration(1)*dt
    bodies%velocity(2) = bodies%velocity(2) + 0.5_dp*bodies%acceleration(2)*dt
    bodies%velocity(3) = bodies%velocity(3) + 0.5_dp*bodies%acceleration(3)*dt

    bodies%internal_energy = bodies%internal_energy + 0.5_dp*bodies%internal_energy_rate*dt
  end subroutine

  !position update
  subroutine drift(bodies, sinks, dt)
    implicit none
    type(particle), intent(inout) :: bodies(:)
    type(sink), intent(inout) :: sinks(:)
    real(dp), intent(in) :: dt

    !x(i+1)
    bodies%position(1) = bodies%position(1) + 0.5_dp*bodies%velocity(1)*dt
    bodies%position(2) = bodies%position(2) + 0.5_dp*bodies%velocity(2)*dt
    bodies%position(3) = bodies%position(3) + 0.5_dp*bodies%velocity(3)*dt

    bodies%internal_energy = bodies%internal_energy + 0.5_dp*bodies%internal_energy_rate*dt
  end subroutine
!-----------------------Simulation Loop Subroutine---------------------------

  subroutine simulate(bodies, sinks)
    implicit none
    type(particle), intent(inout), allocatable :: bodies(:) ! Array of all particles in the simulation.
    type(branch), allocatable :: root           ! The root node of the octree. Allocated and deallocated within the loop.
    type(sink), intent(inout) :: sinks(:)
    real(dp) :: t, dt, dt_candidate, end_time
    real(dp), allocatable :: vel_squared(:)
    integer :: i, number_bodies , t_test                  ! Loop index.
    t_test = 0
    t = 0.0_dp         ! Initialize simulation time.
    end_time = 0.002_dp ! Set simulation end time.
    dt = 0.0000001_dp          ! Set time step size.
    number_bodies = size(bodies) !total number of pariticles

    ! Main simulation loop: continue as long as current time is less than end time.
    do while (t < end_time)
      ! === First Half-Step of Integration ===
      number_bodies = size(bodies)
      print *,"SPH Particles:", number_bodies, "dt :", dt, "time : ", t
      !print *, 'sink mass:', sinks(1)%mass, 'sink pos',sinks(1)%position
      ! 1. Allocate and initialize the root node for tree building.
      allocate(root)
      ! Initialize root node's bounding box based on the min/max positions of all particles.
      ! This ensures the tree spans the entire particle distribution.
      root%center = [(maxval(bodies%position(1)) + minval(bodies%position(1)))/2.0_dp, &
                     (maxval(bodies%position(2)) + minval(bodies%position(2)))/2.0_dp, &
                     (maxval(bodies%position(3)) + minval(bodies%position(3)))/2.0_dp]
      root%size = maxval([(maxval(bodies%position(1)) - minval(bodies%position(1))), &
                          (maxval(bodies%position(2)) - minval(bodies%position(2))), &
                          (maxval(bodies%position(3)) - minval(bodies%position(3)))])
      
      root%n_particles = number_bodies ! Set the number of particles in the root node.
      allocate(root%particles(number_bodies)) ! Allocate space for copies of all particles in the root.
      root%particles = bodies ! Copy current particle states into the root node for tree construction.

      ! 2. Build the octree from the current particle positions.
      call build_tree(root, max_depth, 1) 

      ! 3. Calculate densities for all particles using the newly built tree.
      call get_density(root, bodies)

      ! 4. Calculate gravitational acceleration for all particles using Barnes-Hut.
      ! Initial acceleration is reset before accumulation.
      do i = 1, root%n_particles
        bodies(i)%acceleration = [0.0_dp, 0.0_dp, 0.0_dp]
      end do

      !call particle_gravforces(root, bodies, 0.5_dp) ! theta criterion 0.5.
      call sink_gravforces(bodies, sinks)

      !print *, 'maxacc', sqrt(sum(bodies(1)%acceleration**2))

      ! 5. Calculate SPH accelerations (pressure forces) and internal energy rates.
      ! Initial internal energy rate is reset within `get_SPH` before accumulation.
      call get_SPH(root, bodies)
      ! 6. Update velocities (first half-kick) and internal energies (half-step).
      call kick(bodies, sinks, dt)

      ! 7. Update positions (first half-drift), and reset accelerations/internal_energy_rates for next force calculation.
      call drift(bodies, sinks, dt)

      deallocate(root) ! Deallocate the current tree before rebuilding for the second half.
      ! === Second Half-Step of Integration ===
      ! Recalculate forces based on the new (half-drifted) positions for the second kick.
      ! 8. Re-allocate and re-initialize the root node with the updated positions.
      allocate(root)
      root%center = [(maxval(bodies%position(1)) + minval(bodies%position(1)))/2.0_dp, &
                     (maxval(bodies%position(2)) + minval(bodies%position(2)))/2.0_dp, &
                     (maxval(bodies%position(3)) + minval(bodies%position(3)))/2.0_dp]
      root%size = maxval([(maxval(bodies%position(1)) - minval(bodies%position(1))), &
                          (maxval(bodies%position(2)) - minval(bodies%position(2))), &
                          (maxval(bodies%position(3)) - minval(bodies%position(3)))])

      root%n_particles = number_bodies
      allocate(root%particles(number_bodies))
      root%particles = bodies

      ! 9. Rebuild the octree.
      call build_tree(root, max_depth, 1)

      do i = 1, root%n_particles
        bodies(i)%acceleration = [0.0_dp, 0.0_dp, 0.0_dp]
      end do

      ! 10. Recalculate densities.
      call get_density(root, bodies)

      ! 11. Recalculate gravitational acceleration.
      !call particle_gravforces(root, bodies, 0.5_dp)
      call sink_gravforces(bodies, sinks)
      !print *, 'maxacc', sqrt(sum(bodies(1)%acceleration**2))

      ! 12. Recalculate SPH accelerations and internal energy rates.
      call get_SPH(root, bodies)

      ! 13. Final update of velocities (second half-kick) and internal energies (second half-step).
      ! The positions are also updated here with the second half-drift to complete the full step.
      call kick(bodies, sinks, dt)

      t = t + dt 

      !find next timestep candidate
      allocate(vel_squared(number_bodies))
      do i = 1, number_bodies
        vel_squared(i) = sqrt(sum(bodies(i)%velocity * bodies(i)%velocity)/sum(bodies(i)%acceleration * bodies(i)%acceleration))
      end do
      dt_candidate = minval(vel_squared) * 0.01

      deallocate(vel_squared)

      if (dt_candidate > 2*dt .and. 1.5 * dt < 1) then
        dt = 1.5 * dt
      else if (dt_candidate < 0.5 * dt .and. dt * 0.5 >0.0000001) then
        dt = 0.5 * dt
      end if

      ! Sink accretion and boundary check
      call initiate_sink_accretion(sinks, bodies, root)
      deallocate(root) 

      if (t > 0.3333_dp*end_time.and. t_test==0) then 
        call make_save(bodies, 2)
        t_test = 1
      end if
      if (t > 0.6667_dp*end_time.and. t_test==1) then 
        call make_save(bodies, 3)
        t_test = 2
      end if
    end do

    call make_save(bodies, 4)
  end subroutine simulate
end module octree_module

program barnes_hut
  use octree_module
  implicit none
  character(len=256) :: filename
  type(particle), allocatable :: bodies(:) ! Allocatable array to hold all particles in the simulation
  type(sink), allocatable :: sinks(:)
  ! Initialize the SPH kernel lookup tables (W and dW/dr). This needs to be done once at the start.
  call init_kernel_table()
  call init_grav_kernel_table()

  filename = 'keplerian_ring_big.txt'

  call read_data_from_file(filename,bodies,5000)

  allocate(sinks(1))
  sinks(1)%position = [0.0_dp,0.0_dp,0.0_dp]
  sinks(1)%velocity = [0.0_dp,0.0_dp,0.0_dp]
  sinks(1)%radius = 0.0000000000001_dp
  sinks(1)%mass = 1.0_dp

  ! Start the main simulation loop.
  call simulate(bodies,sinks)

end program barnes_hut
