module SPH

    using DelimitedFiles
    
    include("step_size.jl")
    include("boundary_update.jl")
    include("velocity_derivative.jl")
    include("density_derivative.jl")
    include("save_vtu.jl")

    export main_SPH

    # Simulate Particle Movement
    function main_SPH(d,r_c, dx, h, dt_save,particles, rho_0,gamma,c_0,p_0,Xi,alpha, int_fluid, int_boundary,g,t_end,t,n_dt,dvelx_by_dt,dvely_by_dt,xMin,yMin,Nx,Ny,dx_x,dy_y,n_save,t_ref,save_pos_t)
    
        while t <= t_end # loop over timestep
  
            # display progress
            prog = ceil(100*t/t_end)    
            println("progress = $prog%,time = $t,timestep = $n_dt")
     
            # determine time step size 
            dt = step_size(c_0, particles, h, alpha, d, g)
            t = t + dt
            n_dt = n_dt + 1
             
            # time stepping schemes    
            # velocity-Verlet scheme
                # Calculate acceleration
                if n_dt == 1 
                    particles[int_boundary,8],particles[:,6],particles[:,7] = boundary_update(particles, int_fluid, int_boundary, d,r_c, h, g,xMin,yMin,Nx,Ny,dx_x,dy_y)
                    dvelx_by_dt, dvely_by_dt = velocity_derivative(particles, int_fluid, h, d, r_c,alpha, c_0, g, rho_0, p_0, Xi,gamma,xMin,yMin,Nx,Ny,dx_x,dy_y)
                end
         
                # Velocity Update
                particles[int_fluid,6] = particles[int_fluid,6] + dt/2 * dvelx_by_dt
                particles[int_fluid,7] = particles[int_fluid,7] + dt/2 * dvely_by_dt

                # Position Update
                particles[int_fluid,1] = particles[int_fluid,1] + dt/2 * particles[int_fluid,6]
                particles[int_fluid,2] = particles[int_fluid,2] + dt/2 * particles[int_fluid,7]

                # Density and Pressure update    
                particles[int_boundary,8],particles[:,6],particles[:,7] = boundary_update(particles, int_fluid, int_boundary, d,r_c, h, g,xMin,yMin,Nx,Ny,dx_x,dy_y)
                drho_by_dt = density_derivative(particles, int_fluid, h, d, r_c,rho_0, p_0, Xi, gamma, dx,xMin,yMin,Nx,Ny,dx_x,dy_y)
                particles[int_fluid,5] = particles[int_fluid,5] + dt * drho_by_dt 
                particles[int_fluid,8] = p_0[2] * ((particles[int_fluid,5]/rho_0[2]).^gamma[2] .- 1 ) .+ Xi[2]

                # Position Update
                particles[int_fluid,1] = particles[int_fluid,1] + dt/2 * particles[int_fluid,6]
                particles[int_fluid,2] = particles[int_fluid,2] + dt/2 * particles[int_fluid,7]

                # Velocity update 
                particles[int_boundary,8],particles[:,6],particles[:,7] = boundary_update(particles, int_fluid, int_boundary, d,r_c, h, g,xMin,yMin,Nx,Ny,dx_x,dy_y)
                dvelx_by_dt, dvely_by_dt = velocity_derivative(particles, int_fluid, h, d, r_c,alpha, c_0, g, rho_0, p_0, Xi,gamma,xMin,yMin,Nx,Ny,dx_x,dy_y)
                particles[int_fluid,6] = particles[int_fluid,6] + dt/2 * dvelx_by_dt
                particles[int_fluid,7] = particles[int_fluid,7] + dt/2 * dvely_by_dt
                
            # saving particle data in time saving increments 'dt_save'
             n_save = save_data(t,n_save,dt_save,t_ref,save_pos_t,particles,int_fluid)        
        end
       return particles
    end

    function save_data(t,n_save,dt_save,t_ref,save_pos_t,particles,int_fluid)
        
        if t >= n_save*dt_save*t_ref 
        
            # save time
            save_pos_t[n_save,1] = t
            # save x front
            save_pos_t[n_save,2] = maximum(particles[int_fluid,1])-0.6
                   
            # save fluid particle matrix as .vtu-file
            save_vtu(particles[int_fluid,1:8],n_save,dirname)
            n_save = n_save + 1
            
        end
            pos_csv = open("save_pos_t.csv","w")
            data = save_pos_t
            writedlm(pos_csv,data,',')
            close(pos_csv) 
      return n_save
    end
end    