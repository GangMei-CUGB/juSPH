function create_fluid(dx, f_lowleft, f_lowright, f_upleft, f_upright)

    # create fluid particle vector
    no_rows = ceil(abs(f_upleft[1,2] - f_lowleft[1,2]) / dx)
    no_cols = ceil(abs(f_upleft[1,1] - f_upright[1,1]) / dx)
    no_rows = convert(Int,no_rows)
    no_cols = convert(Int,no_cols)
    fluid = zeros(Float64,no_rows*no_cols,2)

    for n = 1 : no_rows
        for m = 1 : no_cols
            # x-coordinate
            fluid[m+(n-1)*no_cols,1] = f_lowleft[1,1]+(m-1/2)*dx
            # y-coordinate
            fluid[m+(n-1)*no_cols,2] = f_lowleft[1,2]+(n-1/2)*dx
        end
    end
  return fluid
end