include("kernel_function.jl")
include("link_list.jl")

function boundary_update(particles, int_fluid, int_boundary, d,r_c, h, g,xMin,yMin,Nx,Ny,dx_x,dy_y)
    
    # prescribed wall velocity(no slip boundary condition)
    v_wall = zeros(Float64,size(particles,1),2) 
    v_wall[int_boundary,1:2] = particles[int_boundary,4:5] 

    # velocity components
    vel_x = zeros(Float64,size(particles,1),1) 
    vel_y = zeros(Float64,size(particles,1),1) 
    vel_x = particles[:,6]
    vel_y = particles[:,7]
    vel_x = SharedArray(vel_x)
    vel_y = SharedArray(vel_y)
   
    # initialisation pressure
    pressure = SharedArray{Float64}(length(int_boundary),1)
      
    # Nearest neighbor particle search(link list search method) 
    # Search for particles adjacent to the boundary particle in the fluid particle range 
        # Establish boundary particles' neighbors
        # Estimate the number of adjacent particles per boundary particle (estimated by maximum number of adjacent particles)  
        Max_number = floor(Int,(r_c*2/h+1)^2)
        boundaryParticles_neighbors = SharedArray{Int}(Max_number,length(int_boundary))
        # link list search
        numBins,adjacentBins,Gridnum,ParticlesList = link_list(particles,xMin,yMin,Nx,Ny,dx_x,dy_y)
        # Build the list of neighbors for each boundary particle
        @sync @distributed for z = 1:numBins    
            # Get all boundary particles in bin z 
            if z == 1
                local AList = filter(b->b>length(int_fluid),ParticlesList[1:Gridnum[z]])
            else
                local AList = filter(b->b>length(int_fluid),ParticlesList[Gridnum[z-1]+1:Gridnum[z]]) 
            end
            # Get vector of bin neighorIDs adjacent to bin z
            local AdjacentBins = filter(!(isequal(0)), adjacentBins[:,z])
            local w = copy(AdjacentBins)
            push!(w,z)
            local consider = Int[]  
            @inbounds for a = 1 : length(w)
                # Get all fluid particles in bin w
                if w[a] == 1     
                    append!(consider,filter(b->b<=length(int_fluid),ParticlesList[1:Gridnum[w[a]]]))
                else
                    append!(consider,filter(b->b<=length(int_fluid),ParticlesList[Gridnum[w[a]-1]+1:Gridnum[w[a]]])) 
                end
            end

            @inbounds for k = 1:length(AList) # for each particle in bin z
                local neighbors = Int[]
                @inbounds for j = 1:length(consider) # for each possible neighboring particle
                    xDiff = particles[AList[k],1] - particles[consider[j],1]
                    yDiff = particles[AList[k],2] - particles[consider[j],2]
                    dist = sqrt(xDiff^2 + yDiff^2)
                    if dist <= r_c
                        # Particle consider[j] is a neighbor of AList[k]
                        push!(neighbors,consider[j])                   
                    end
                end

                for i = 1:length(neighbors)
                    boundaryParticles_neighbors[i,findfirst(in(AList[k]),int_boundary)] = neighbors[i]
                end 
            end
        end

    # calculate wall pressure and velocity
    @inbounds @sync @distributed for n = 1:length(int_boundary) # over all boundary particles
        w = int_boundary[n]
        # initialisation wall pressure
        sum_pW = 0.0
        sum_rhorW = 0.0
        sum_W = 0.0
        # initialisation wall velocity
        sum_vWx = 0.0
        sum_vWy = 0.0 
        
        local boundaryParticles_Neighors = filter(!(isequal(0)), boundaryParticles_neighbors[:,n])
        
        @inbounds for m = 1:length(boundaryParticles_Neighors) # interaction with all fluid particles in viscinity
            f = boundaryParticles_Neighors[m]

            # distance between particles
            drx = particles[w,1] - particles[f,1]
            dry = particles[w,2] - particles[f,2]
            rad = sqrt(drx^2 + dry^2)
            
            # nondimensional distance
            q = rad / h
            # Kernel function
            W = kernel_function(d, h, q)

            # build up pressure of wall particle
            sum_pW = sum_pW + particles[f,8] * W
            sum_rhorW = sum_rhorW + particles[f,5] * W *(drx * g[1] + dry * g[2])
            sum_W = sum_W + W
            
            
            # building up the SPH average for wall velocity(no slip boundary condition)
            sum_vWx = sum_vWx + vel_x[f] * W
            sum_vWy = sum_vWy + vel_y[f] * W
            
        end

        if sum_W == 0
            pressure[n] = 0
        else
            # combining terms to get pressure of wall particle 
            pressure[n] = ( sum_pW + sum_rhorW ) / sum_W 

            # calculate wall velocity(no slip boundary condition)
            vel_x[w] = 2 * v_wall[w,1] - sum_vWx / sum_W
            vel_y[w] = 2 * v_wall[w,2] - sum_vWy / sum_W
            
        end        
    end

  return pressure,vel_x,vel_y

end