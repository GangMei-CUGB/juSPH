include("SPH.jl")
include("dam_break.jl")
include("save_vtu.jl")

using .SPH
using DelimitedFiles

# Initialisation
# setup variables and domain etc...
    # dimensions
    const d = 2
    # distance between particles
    const dx = 0.006
    # smoothing length
    const h = dx
    # cutoff radius
    const r_c = 3 * h     # Quintic Spline
    # viscosity coefficient value
    const alph = 0.02 
    # dambreak testcase
    particles, rho_0,gamma,c_0,p_0,Xi,my,alpha,a_wall, int_fluid, int_boundary, g, Height = dam_break(dx, d, alph)    
    # acceleration 
    dvelx_by_dt = zeros(Float64,length(int_fluid),1)
    dvely_by_dt = zeros(Float64,length(int_fluid),1)

    # reference values  
    H = maximum(particles[int_fluid,2]) - minimum(particles[int_fluid,2]) + dx # height of fluid
    if abs(H - Height) >= 1.0e-6
        println("wrong height specified,please input dx again")
    end    
    const v_ref = sqrt(2 *abs(g[2]) * H )
    const t_ref = H / v_ref 
    const p_ref = rho_0 * abs(g[2]) * H

    # Set grid size
    const xMax = 1.61+3*dx
    const yMax = 0.8+dx
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

# Simulate Particle Movement
@time begin
    particles = main_SPH(d, r_c, dx, h, dt_save,particles, rho_0,gamma,c_0,p_0,Xi,alpha, int_fluid, int_boundary,g,t_end,t,n_dt,dvelx_by_dt,dvely_by_dt,xMin,yMin,Nx,Ny,dx_x,dy_y,n_save,t_ref,save_pos_t)
end

# final savings
# save boundary particle matrix as .vtu-file
    save_vtu(particles[int_boundary,1:8],0,dirname)

    println("progress = 100%, Simulation completed \n")
