include("link_list.jl")
include("kernel_derivative.jl")

function velocity_derivative(particles, int_fluid, h, d, r_c,alpha, c_0, g, rho_0, p_0, Xi,gamma,xMin,yMin,Nx,Ny,dx_x,dy_y)
    
    # initialisation
    dvelx_by_dt = zeros(Float64,length(int_fluid),1)
    dvely_by_dt = zeros(Float64,length(int_fluid),1)
    dx = h        # h is smoothing length
    epsilon = 0.01 # parameter to avoid zero denominator 
      
    # Nearest neighbor particle search(link list search method)     
    # Search for particles adjacent to fluid particles in the range of all particles 
        # establish fluid particles' neighbors
        # Estimate the number of adjacent particles per fluid particle (estimated by maximum number of adjacent particles)  
        Max_number = floor(Int,(r_c*2/h+1)^2)
        fluidParticles_neighbors = zeros(Int,Max_number,length(int_fluid))
        # link list search
        numBins,adjacentBins,Gridnum,ParticlesList = link_list(particles,xMin,yMin,Nx,Ny,dx_x,dy_y)
        # Build the list of neighbors for each boundary particle
        for z = 1:numBins    
            # Get all boundary particles in bin z 
            if z == 1
                local AList = filter(b->b<=length(int_fluid),ParticlesList[1:Gridnum[z]])
            else
                local AList = filter(b->b<=length(int_fluid),ParticlesList[Gridnum[z-1]+1:Gridnum[z]]) 
            end
            # Get vector of bin neighorIDs adjacent to bin z
            local AdjacentBins = filter(!(isequal(0)), adjacentBins[:,z])
            local w = copy(AdjacentBins)
            push!(w,z)
            local consider = Int[] 
            for a = 1 : length(w)
                # Get all fluid particles in bin w
                if w[a] == 1     
                    append!(consider,ParticlesList[1:Gridnum[w[a]]])
                else
                    append!(consider,ParticlesList[Gridnum[w[a]-1]+1:Gridnum[w[a]]])
                end
            end

            for k = 1:length(AList) # for each particle in bin z
                local neighbors = Int[]
                for j = 1:length(consider) # for each possible neighboring particle
                    xDiff = particles[AList[k],1] - particles[consider[j],1]
                    yDiff = particles[AList[k],2] - particles[consider[j],2]
                    dist = sqrt(xDiff^2 + yDiff^2)
                    if dist <= r_c && AList[k] != consider[j]
                        # Particle consider[j] is a neighbor of AList[k]                        
                        push!(neighbors,consider[j])               
                    end
                end
                for i = 1:length(neighbors)
                    fluidParticles_neighbors[i,AList[k]] = neighbors[i]
                end 
            end
        end
    
    # calculate the acceleration of each fluid particle 
    for a = 1 : length(int_fluid)   # loop over all fluid particles
        local fluidParticles_Neighors = filter(!(isequal(0)), fluidParticles_neighbors[:,a])
        # calculate interaction of one fluid particle with all particles within range r_c
        for m = 1:length(fluidParticles_Neighors) 
            b = fluidParticles_Neighors[m]

            # distance between particles
            drx = particles[a,1] - particles[b,1]
            dry = particles[a,2] - particles[b,2]
            rad = sqrt(drx^2 + dry^2)
               
            # nondimensional distance
            q = rad / h
            # derivative of Kernel
            der_W = kernel_derivative(d, h, q) / h
                  
            # velocity difference between particles
            dvx = particles[a,6] - particles[b,6]
            dvy = particles[a,7] - particles[b,7]
            
            # pressure of particles
            p_a = particles[a,8]
            p_b = particles[b,8]
            
            # particle type
            flag_a = particles[a,3]
            flag_a = convert(Int,flag_a)
            flag_b = particles[b,3]
            flag_b = convert(Int,flag_b)
            
            # densities, mass
            rho_a = particles[a,5]
            m_a   = particles[a,4]
                if flag_b == 1 # if b is a boundary particle
                    # calculate density
                    rho_b = rho_0[flag_a] * ((particles[a,8] - Xi[flag_a]) / p_0[flag_a] + 1) .^(1/gamma[flag_a])
                    m_b = rho_0[flag_a] * dx * dx     
                else # if b is a fluid particle
                    rho_b = particles[b,5]
                    m_b   = particles[b,4]
                end
            rho_ab = 0.5 * (rho_a + rho_b)   
            
            # momentum equation Pressure Gradient part
            p_ab = (rho_b * p_a + rho_a * p_b) / (rho_a + rho_b)   
            pressure_fact = - 1/m_a * ((m_a/rho_a)^2 + (m_b/rho_b)^2) * (p_ab) * der_W
             
            # Monaghan formulation
            #pressure_fact = - m_b * (p_a / rho_a^2 + p_b / rho_b^2) * der_W  # momentum equation Pressure Gradient part
            
            # acceleration due to pressure gradient (eq7) 
            dvelx_by_dt[a] = dvelx_by_dt[a] + pressure_fact * drx / rad
            dvely_by_dt[a] = dvely_by_dt[a] + pressure_fact * dry / rad
                                   
            # if free slip condition aplies only consider particles which are not fluid particles                       
            # artificial viscosity 
            if flag_b != 1 # so if particle b is a fluid
                alpha_ab = 0.5 * (alpha[flag_a] + alpha[flag_b])
                c_ab = 0.5 * (c_0[flag_a] + c_0[flag_b])
            else
                alpha_ab = alpha[flag_a]
                c_ab = c_0[flag_a]
            end
                    
            if (drx * dvx + dry * dvy) < 0
                visc_art_fact = m_b * alpha_ab * h * c_ab * (dvx * drx + dvy * dry)/(rho_ab * (rad^2 + epsilon * h^2)) * der_W
            else
                visc_art_fact = 0               
            end                   
                    
            dvelx_by_dt[a] = dvelx_by_dt[a] + visc_art_fact * drx / rad
            dvely_by_dt[a] = dvely_by_dt[a] + visc_art_fact * dry / rad
        end  
            # Gravity
            dvelx_by_dt[a] = dvelx_by_dt[a] + g[1]
            dvely_by_dt[a] = dvely_by_dt[a] + g[2]
    end

   return dvelx_by_dt, dvely_by_dt

end