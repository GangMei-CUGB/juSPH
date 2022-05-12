function kernel_function(d, h, q)
# Quintic Spline
    # cutoff radius r_c = 3 * h;

    # normalisation parameter alpha_d
    if d == 2
        #alpha_d = 7 / (478 * pi * h^2); # for 2D
        alpha_d = 0.004661441847880 / (h * h)
    elseif d == 3
        alpha_d = 1 / (120 * pi * h^3) # for 3D
    elseif d == 1
        alpha_d = 1 / (120 * h)
    end

    # Weighting function
    if q < 3 && q >= 2
        W = alpha_d * (3-q)^5
    elseif q < 2 && q >= 1
        W = alpha_d * ((3-q)^5 - 6 * (2-q)^5)
    elseif q < 1
        W = alpha_d * ((3-q)^5 - 6 * (2-q)^5 + 15 * (1-q)^5)
    elseif q >= 3
        W = 0
    end

return W

end