# Calculate the time for linked list search to search and store neighborhood particles

using Distributed
using SharedArrays
addprocs()

include("dam_break.jl")
include("link_list.jl")
@everywhere include("get_adjacentBins.jl")  

# Initialisation
# setup variables and domain etc...
    # dimensions
    const d = 2
    # distance between particles
    const dx = 0.006
    # smoothing length
    const h = dx
    # cutoff radius
    const r_c = 3 * h
    # viscosity coefficient value
    const alph = 0.1  
    # dambreak testcase
    particles, rho_0,gamma,c_0,p_0,Xi,my,alpha,a_wall, int_fluid, int_boundary, G, Height = dam_break(dx, d, alph)    
    # acceleration 
    dvelx_by_dt = zeros(Float64,length(int_fluid),1)
    dvely_by_dt = zeros(Float64,length(int_fluid),1)
    
    # reference values  
    H = maximum(particles[int_fluid,2]) - minimum(particles[int_fluid,2]) + dx # height of fluid
    if abs(H - Height) >= 1.0e-6
        println("wrong height specified,please input dx again")
    end    
    const v_ref = sqrt(2 *abs(G[2]) * H )
    const t_ref = H / v_ref 
    const p_ref = rho_0 * abs(G[2]) * H
    
    # Set grid size
    const xMax = 1.61+3*dx
    const yMax = 0.8+3*dx
    const xMin = 0.0-3*dx
    const yMin = 0.0-3*dx
    const Nx = cld((xMax-xMin),r_c)
    const Ny = cld((yMax-yMin),r_c)
    const dx_x = (xMax-xMin)/Nx
    const dy_y = (yMax-yMin)/Ny
    
    # time initialisation
    t = 0      # time
    n_dt = 0   # number of current timestep  
    const time_nondim = 5    # nondimensional running time (t/t_ref)
    const t_damp_fact = 0    # damping coefficiant
    const t_damp = t_damp_fact * t_ref  # time over which gravitiy is introduced
    const t_end = time_nondim * t_ref
        
    # saving initialisation
    const dt_save = 0.01   # saving time increment
    n_save = Int64(ceil(time_nondim/dt_save))
    if n_save >= 1001
        println("Too many saving steps specified!")
    end
    save_pos_t = zeros(Float64,n_save,2)
    n_save = 1
    
    # initial particles values
    particles_ini = particles
    particles = SharedArray(particles) 
       
    # Search for particles adjacent to fluid particles in the range of all particles 
    function search_neighbors(r_c,h,int_fluid,particles,xMin,yMin,Nx,Ny,dx_x,dy_y)
        # Establish fluid particles' neighbors
        # Estimate the number of adjacent particles per fluid particle (estimated by maximum number of adjacent particles)  
        Max_number = floor(Int,(r_c*2/h+1)^2)
        fluidParticles_neighbors = SharedArray{Int}(Max_number,length(int_fluid))
        # link list search
        numBins,adjacentBins,Gridnum,ParticlesList = link_list(particles,xMin,yMin,Nx,Ny,dx_x,dy_y)
        # Build the list of neighbors for each boundary particle
        @sync @distributed for z = 1:numBins    
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
            @inbounds for a = 1 : length(w)
                # Get all fluid particles in bin w
                if w[a] == 1     
                    append!(consider,ParticlesList[1:Gridnum[w[a]]])
                else
                    append!(consider,ParticlesList[Gridnum[w[a]-1]+1:Gridnum[w[a]]]) 
                end
            end

            @inbounds for k = 1:length(AList) # for each particle in bin z
                local neighbors = Int[]
                @inbounds for j = 1:length(consider) # for each possible neighboring particle
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
    end

# Nearest neighbor particle search(link list search method)
@time fluidParticles_neighbors = search_neighbors(r_c,h,int_fluid,particles,xMin,yMin,Nx,Ny,dx_x,dy_y)    

