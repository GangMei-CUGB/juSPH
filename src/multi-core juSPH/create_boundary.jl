function create_boundary(dx, b_lowleft, b_lowright, b_upleft, b_upright)

    # number of boundary particles for kernel support
    k = 3    # quintic spline
    
    # create boundary particles vector
    no_rows = ceil(abs(b_upleft[1,2] - b_lowleft[1,2]) / dx)
    no_cols = ceil(abs(b_upleft[1,1] - b_upright[1,1]) / dx + 2*k)
    no_rows = convert(Int,no_rows)
    no_cols = convert(Int,no_cols)
    boundary = zeros(Float64,2*no_rows*k+k*no_cols,2)

    for n = 1 : no_rows
        for m = 1 : 2*k
        # x-coordinate
        if m < k+1
            boundary[m+(n-1)*2*k ,1] = b_upleft[1,1]-(k*dx-dx/2)+(m-1)*dx
        else
            boundary[m+(n-1)*2*k ,1] = b_upright[1,1]+dx/2+(m-k-1)*dx
        end
        # y-coordinate
        boundary[m+(n-1)*2*k ,2] = b_lowleft[1,2] + (n-1/2) * dx
        end
    end

    for n = 1 : k
        for m = 1 : no_cols
        # x-coordinate
        boundary[2*no_rows*k+m+(n-1)*no_cols ,1] = b_lowleft[1,1] -(k*dx-dx/2) + (m-1) * dx
        # y-coordinate
        boundary[2*no_rows*k+m+(n-1)*no_cols ,2] = b_lowleft[1,2] - (n-1/2) * dx
        end
    end
    
return boundary
end