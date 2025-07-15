module octree_module
  implicit none
  ! Global Parameters:
  real, parameter :: G = 6.67430e-12 !should be -12        ! Gravitational constant (in SI units)
  real, parameter :: softening = 300.0
  integer, parameter :: dp = kind(1.0d0)    ! Defines double precision
  integer, parameter :: nq = 1000, max_depth = 10000          ! Number of samples for the SPH kernel lookup tables.
                                           ! This determines the resolution of the pre-computed kernel values.
  real(dp), allocatable :: w_table(:), dw_table(:) ! Allocatable arrays to store pre-computed SPH kernel (W)
                                                   ! and its derivative (dW/dr) values.
  real(dp), parameter :: dq = 2.0_dp / nq 
  real(dp), parameter :: smoothing = 100.0_dp 

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

  ! Represents a node (branch) in the Barnes-Hut octree.
  type :: branch
    real(dp), dimension(3) :: center      
    real(dp) :: size                     
    integer :: n_particles              
    type(Particle), allocatable :: particles(:) 
    real(dp) :: mass_total               
    real(dp), dimension(3) :: mass_center 
    type(branch), allocatable :: children(:) ! Array of 8 child branches

  end type branch

  contains

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
  recursive subroutine navigate_tree(node, body, theta)
  implicit none
  type(branch), intent(in) :: node
  type(particle), intent(inout) :: body(:)
  real(dp), intent(in) :: theta
  real(dp) :: dist, d2 
  real(dp), dimension(3) :: direction
  integer :: i, k

  ! Loop over the (single) particle provided in the 'body' array slice.
  do i = 1, size(body)
    direction = body(i)%position - node%mass_center ! Vector pointing from particle to the node's center of mass
    d2 = sum(direction**2) + softening**2 
    dist = sqrt(d2)  

    ! Barnes-Hut criterion:
    if ((node%size / dist) < theta .or. .not. allocated(node%children)) then
      ! Calculate gravitational acceleration if the node has mass and distance is positive.
      if (node%mass_total > 0.0_dp .and. dist > 0.0_dp) then
        body(i)%acceleration = body(i)%acceleration - (G * node%mass_total * direction / (dist**3))
      end if
    else if (allocated(node%children)) then
      ! If the criterion is not met and the node has children, recurse into each child node.
      do k = 1, size(node%children)
        if (node%children(k)%n_particles > 0) then
          call navigate_tree(node%children(k), body(i:i), theta) ! Recurse, passing the same particle
        end if
      end do
    end if
  end do
  end subroutine navigate_tree

  ! Subroutine: get_SPH
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
      nr = (body%position - node%particles(1)%position) ! Vector from neighbor to the 'body' particle
      dr = sqrt(sum(nr**2))                             ! Distance between them
      !if (dr == 0.0_dp) return

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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulate(bodies)
    implicit none
    type(particle), intent(inout) :: bodies(:) ! Array of all particles in the simulation.
    type(branch), allocatable :: root           ! The root node of the octree. Allocated and deallocated within the loop.
    real(dp) :: t, dt, dt_candidate, end_time
    real(dp), allocatable :: vel_squared(:)
    integer :: i, io, number_bodies                           ! Loop index.

    t = 0.0_dp         ! Initialize simulation time.
    end_time = 10000.0_dp ! Set simulation end time.
    dt = 1.0_dp          ! Set time step size.
    number_bodies = size(bodies) !total number of pariticles
    allocate(vel_squared(number_bodies))

    ! Main simulation loop: continue as long as current time is less than end time.
    do while (t < end_time)
      ! === First Half-Step of Integration ===

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

      ! 5. Calculate gravitational acceleration for all particles using Barnes-Hut.
      ! Initial acceleration is reset before accumulation.
      do i = 1, root%n_particles
        bodies(i)%acceleration = [0.0_dp, 0.0_dp, 0.0_dp]
      end do

      call navigate_tree(root, bodies, 0.5_dp) ! theta criterion 0.5.

      ! 6. Calculate SPH accelerations (pressure forces) and internal energy rates.
      ! Initial internal energy rate is reset within `get_SPH` before accumulation.
      call get_SPH(root, bodies)

      ! 7. Update velocities (first half-kick) and internal energies (half-step).
      do i = 1, root%n_particles
        bodies(i)%velocity = bodies(i)%velocity + (bodies(i)%acceleration)*dt/2.0_dp
        bodies(i)%internal_energy = bodies(i)%internal_energy + (bodies(i)%internal_energy_rate)*dt/2.0_dp
      end do

      ! 8. Update positions (first half-drift), and reset accelerations/internal_energy_rates for next force calculation.
      do i = 1, root%n_particles
        bodies(i)%position = bodies(i)%position + (bodies(i)%velocity)*dt/2.0_dp
        bodies(i)%acceleration = [0.0_dp, 0.0_dp, 0.0_dp]
        bodies(i)%internal_energy_rate = 0.0_dp
      end do

      deallocate(root) ! Deallocate the current tree before rebuilding for the second half.

      ! === Second Half-Step of Integration ===
      ! Recalculate forces based on the new (half-drifted) positions for the second kick.

      ! 9. Re-allocate and re-initialize the root node with the updated positions.
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

      ! 10. Rebuild the octree.
      call build_tree(root, max_depth, 1)

      ! 11. Recalculate densities.
      call get_density(root, bodies)

      ! 13. Recalculate gravitational acceleration.
      call navigate_tree(root, bodies, 0.5_dp)

      ! 14. Recalculate SPH accelerations and internal energy rates.
      call get_SPH(root, bodies)

      ! 15. Final update of velocities (second half-kick) and internal energies (second half-step).
      ! The positions are also updated here with the second half-drift to complete the full step.
      do i = 1, root%n_particles
        bodies(i)%velocity = bodies(i)%velocity + (bodies(i)%acceleration)*dt/2.0_dp
        bodies(i)%internal_energy = bodies(i)%internal_energy + (bodies(i)%internal_energy_rate)*dt/2.0_dp
        bodies(i)%position = bodies(i)%position + (bodies(i)%velocity)*dt/2.0_dp ! This completes the full position update (first half was above).
      end do

      deallocate(root) 
      t = t + dt 

      !find next timestep candidate
      do i = 1, number_bodies
        vel_squared(i) = sqrt(sum(bodies(i)%velocity * bodies(i)%velocity)/sum(bodies(i)%acceleration * bodies(i)%acceleration))
      end do
      dt_candidate = minval(vel_squared) * 0.5

      if (dt_candidate > 2*dt .and. 1.5 * dt < 150) then
        dt = 1.5 * dt
      else if (dt_candidate < 0.5 * dt) then
        dt = 0.5 * dt
      end if
      print *,"dt :", dt, "time : ", t

    end do

    deallocate(vel_squared)

    open(newunit=io, file="log-150725.txt", status="new", action = "write")
    write(io, *) 'x  ','y  ','z  ', 'vx  ','vy ','vz ','energy ','density  ', 'mass '
    do i = 1, number_bodies
      write(io, *)  bodies(i)%position(1),bodies(i)%position(2),bodies(i)%position(3),bodies(i)%velocity(1),bodies(i)%velocity(2),bodies(i)%velocity(3), bodies(i)%internal_energy, bodies(i)%density, bodies(i)%mass
    end do  
    close(io)
    print *, maxval([(maxval(bodies%position(1)) - minval(bodies%position(1))), &
                          (maxval(bodies%position(2)) - minval(bodies%position(2))), &
                          (maxval(bodies%position(3)) - minval(bodies%position(3)))])
  end subroutine simulate
end module octree_module

program barnes_hut
  use octree_module  
  implicit none      

  integer :: n, total_particles   ! Loop counter and total number of particles
  type(particle), allocatable :: bodies(:) ! Allocatable array to hold all particles in the simulation
  logical :: exists

  ! Initialize the SPH kernel lookup tables (W and dW/dr). This needs to be done once at the start.
  call init_kernel_table()

  ! Define the total number of particles for the simulation.
  total_particles = 2000
  ! Allocate memory for the array of particles.
  allocate(bodies(total_particles))

  ! Initialize the properties of each particle.
  do n = 1, total_particles
    call random_number(bodies(n)%position) ! Generate random numbers (0 to 1) for initial positions.
    ! Scale and shift positions to be within a cube (e.g., from 0 to 12).
    bodies(n)%position = 20000.0_dp * (bodies(n)%position - [0.5,0.5,0.5])
    bodies(n)%number = n
    bodies(n)%mass = 10000.0_dp           ! Assign a mass to each particle.
    call random_number(bodies(n)%velocity) ! Generate random numbers for initial velocities.
    ! Scale velocities to be very small, effectively starting from rest or very slow movement.
    bodies(n)%velocity = 0.100_dp * (bodies(n)%velocity- [0.5,0.5,0.5])

    bodies(n)%pressure = 1.0_dp
    bodies(n)%internal_energy = 500.0_dp 
    bodies(n)%density = 0.0_dp         
  end do

  !----------------------------------
  !------------take file------------- 
  if (exists) then
    
  end if
  ! Start the main simulation loop.
  call simulate(bodies)

end program barnes_hut
