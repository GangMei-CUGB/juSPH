include("kernel_derivative.jl")

function density_derivative(particles, int_fluid, h, d, r_c,rho_0, p_0, Xi, gamma, dx,xMin,yMin,Nx,Ny,dx_x,dy_y)

    # initialisation
    drho_by_dt  = zeros(Float64,length(int_fluid),1)

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
        
    for a = 1 : length(int_fluid)   # loop over all fluid particles
        local fluidParticles_Neighors = filter(!(isequal(0)), fluidParticles_neighbors[:,a])
        
        # calculate interaction of one fluid particle with all other particles
        for m = 1:length(fluidParticles_Neighors) 
            b = fluidParticles_Neighors[m]
    
            flag_b = particles[b,3]
            flag_b = convert(Int,flag_b)
            flag_a = particles[a,3]
            flag_a = convert(Int,flag_a)
            
            # distance between particles
            drx = particles[a,1] - particles[b,1]
            dry = particles[a,2] - particles[b,2]
            rad = sqrt(drx^2 + dry^2)
            
            # nondimensional distance
            q = rad / h
            # kernel and derivative values
            DWab = kernel_derivative(d, h, q) / h
            
            # kernel derivative with respect to x_a
            Fab = [(drx/rad);(dry/rad)]*DWab
            
            # densities, mass
            rho_a = particles[a,5]
            if particles[b,3] == 1 # if b boundary particle
                # calculate density of boundary particles, depending on interacting fluid particle
                rho_b = rho_0[flag_a] * ((particles[a,8] - Xi[flag_a]) / p_0[flag_a] + 1) .^(1/gamma[flag_a])   # calculate density from pressure for wall particles 
                m_b = rho_0[flag_a] *dx*dx
                # velocity difference between particles(vel_wall from boundary_update)
                dvx = particles[a,6] - particles[b,6]
                dvy = particles[a,7] - particles[b,7]

            else # straightforward if fluid particle
                rho_b = particles[b,5]
                m_b = particles[b,4]
                # velocity difference between particles
                dvx = particles[a,6] - particles[b,6]
                dvy = particles[a,7] - particles[b,7]
            end
                       
            # continuity equation       
            drho_by_dt[a] = drho_by_dt[a] + rho_a * m_b / rho_b * [dvx;dvy]'*Fab
               
        end
    end
  return drho_by_dt
end
    