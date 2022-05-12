function step_size(c_0, particles, h, alpha, d, g)

    dt = zeros(Float64,3,1)

    # CFL condition 
    c_max = maximum(c_0)
    vel = zeros(Float64,length(particles),1)
    vel = sqrt.(particles[:,6].^2 + particles[:,7].^2) # velocity magnitude of each particle
    v_max = maximum(vel)
    dt[1] = 0.25 * h / (c_max + (v_max))
    
    # viscous condition 
    my = 0.5/(d+2) * maximum(alpha) * h * c_max  # artificial
    dt[2] = 0.125 * h^2 / my
    
    # body force condition 
    function norm(x,p=2)  
        sum(x.^p)^(1/p)  
    end
    dt[3] = 0.25 * sqrt(h/norm(g))

    # global timestep
    dt = minimum(dt)
return dt
end