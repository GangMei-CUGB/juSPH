function kernel_derivative(d, h, q)
# Quintic Spline
    # cutoff radius r_c = 3 * h

    # normalisation parameter
    if d == 2
        #alpha_d = -5 * 7 / (478 * pi * h^2)  # for 2D
        alpha_d = -0.0233072092393989 / (h * h)
    elseif d == 3
        alpha_d = -5 / (120 * pi * h^3) # for 3D
    elseif d == 1
        alpha_d = -5 / (120 * h)
    end

    # derivative of Weighting function
    if q < 3 && q >= 2
        der_W = alpha_d * (3-q)^4
    elseif q < 2 && q >= 1
        der_W = alpha_d * ((3-q)^4 - 6 * (2-q)^4)
    elseif q < 1
        der_W = alpha_d * ((3-q)^4 - 6 * (2-q)^4 + 15 * (1-q)^4)
    elseif q >= 3
        der_W = 0
    end
   
   return der_W
end