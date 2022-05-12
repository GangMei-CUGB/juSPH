function boundary_initialisation(particles, phase, flag,v_wall,vel)

    int=collect(size(particles,1)+1:size(particles,1)+size(phase,1))
    particle=zeros(Float64,length(int),8)
   
    for n = 1 : length(int)
        particle[n, 1] = phase[n,1]              # x coordinates of phase particles
        particle[n, 2] = phase[n,2]              # y coordinates of phase particles
        particle[n, 3] = flag                    # integer flag, 1 for boundary
        particle[n, 4] = v_wall[1]               # the prescribed wall velocity x-component
        particle[n, 5] = v_wall[2]               # the prescribed wall velocity y-component
        particle[n, 6] = vel[1]                  # initial wall velocity x-component
        particle[n, 7] = vel[2]                  # initial wall velocity x-component
        particle[n, 8] = 0.0                     # pressure (will be calculated from density)
    end
    particles = [particles;particle]
return particles,int
end