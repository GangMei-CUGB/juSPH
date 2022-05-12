function fluid_initialisation(particles, phase, flag, rho_0, dx, d, vel)

    int=collect(size(particles,1)+1:size(particles,1)+size(phase,1))
    particles=zeros(Float64,size(particles,1)+length(int),8)

    particles[int, 1] .= phase[int,1]              # x coordinates of phase particles
    particles[int, 2] .= phase[int,2]              # y coordinates of phase particles
    particles[int, 3] .= flag                      # integer flag, 1 for boundary
    particles[int, 5] .= rho_0[flag]               # density
    particles[int, 4] .= particles[int,5] * dx^d   # mass of particles is fixed
    particles[int, 6] .= vel[1]                    # initial velocity x-component
    particles[int, 7] .= vel[2]                    # initial velocity y-component
    particles[int, 8] .= 0.0                       # pressure (will be calculated from density)
      
  return particles,int
end