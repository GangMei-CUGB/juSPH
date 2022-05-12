include("create_fluid.jl")
include("fluid_initialisation.jl")
include("boundary_initialisation.jl")
include("create_boundary.jl")

function dam_break(dx, d, alph)
# Parameter set up 
    # dambreak testcase        
    g = 1 * [ 0, -9.81]
    Height = 0.3
    v_max = sqrt(2*abs(g[2]) * Height) 

# create particle matrix and fluid properties
    # rectangular domain
    # first fluid phase
    # specify coordinates of edges for fluid
    f_lowleft =  [-0.0  0.0]
    f_lowright = [ 0.6  0.0]
    f_upleft =   [-0.0  0.3]
    f_upright =  [ 0.6  0.3]
    
    # define parameter
    rho_0 = zeros(Float64,2,1)
    gamma = zeros(Int64,2,1)
    c_0   = zeros(Float64,2,1)
    p_0   = zeros(Float64,2,1)
    Xi    = zeros(Float64,2,1)
    my    = zeros(Float64,2,1)
    alpha = zeros(Float64,2,1)

    # properties of fluid (will be stored)
    flag = 2
    rho_0[flag] = 1000.0                                # density    
    gamma[flag] = 7                                     # pressure exponent
               # artificial speed of sound, c_0 = 10 * v_max = 10 * sqrt(2*g*H)
    c_0[flag] = 10*v_max
    p_0[flag] = rho_0[flag] * c_0[flag]^2 / gamma[flag] # reference pressure 
    Xi[flag] = 0.0 * p_0[flag]                          # background pressure
    my[flag] = 0.01                                     # viscosity
    alpha[flag] = alph                                  # artificial visc factor
    # initial velocities
    vel = [0.0,0.0]
    # create particle matrix
    fluid = Float64[]
    fluid = create_fluid(dx, f_lowleft, f_lowright, f_upleft, f_upright)   
    particles = Float64[]
    particles,int_fluid = fluid_initialisation(particles, fluid, flag, rho_0, dx, d, vel)

    # specify coordinates of edges for boundary
    width = dx*ceil(1.61/dx)
    b_lowleft =  [-0.0  0.0]
    b_lowright = [width  0.0]
    b_upleft =   [-0.0  0.8]
    b_upright =  [width  0.8]
    # properties
    v_wall = [0.0,0.0]  # prescribed wall velocity
    a_wall = [0.0,0.0]  # wall acceleration
    flag = 1            # for boundary
    # set artificial viscosity of boundary to 0
    alpha[flag] = 0.0
    # create boundary matrix
    boundary = Float64[]
    boundary = create_boundary(dx, b_lowleft, b_lowright, b_upleft, b_upright)
    particles, int_boundary = boundary_initialisation(particles, boundary, flag,v_wall,vel)

# initial pressure, hydrostatic
    particles[int_fluid,8] .= rho_0[2] .* abs(g[2]) .* (Height .- particles[int_fluid,2])
    particles[int_fluid,5] = rho_0[2] .* ( (particles[int_fluid,8] .- Xi[2]) / p_0[2] .+ 1) .^(1/gamma[2])
  
  return particles, rho_0,gamma,c_0,p_0,Xi,my,alpha,a_wall, int_fluid, int_boundary, g,Height
  
end