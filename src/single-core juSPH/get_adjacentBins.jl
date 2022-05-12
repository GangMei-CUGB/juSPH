function get_adjacentBins(binNum,binNumMax,M,N)
    # return the getAdjacentBins for a certain bin
    # M*N number of bins
        adjacentBins = Int[]
        # Left bins
        if binNum <= N
            push!(adjacentBins , binNum + N)
            # Bottom left corner
            if binNum == 1
            push!(adjacentBins , binNum + 1)
            push!(adjacentBins , binNum + N + 1)
            # Top left corner
            elseif binNum == N
            push!(adjacentBins , binNum - 1)
            push!(adjacentBins , binNum + N - 1)
            # Middle of left bins
            else
            push!(adjacentBins , binNum - 1)
            push!(adjacentBins , binNum + 1)
            push!(adjacentBins , binNum + N - 1)
            push!(adjacentBins , binNum + N + 1)
            end
        # Right bins
        elseif binNum > binNumMax - N
            push!(adjacentBins , binNum - N)
            # Right Bottom corner
            if binNum == binNumMax - N + 1
            push!(adjacentBins , binNum + 1)
            push!(adjacentBins , binNum - N + 1)
             # Right top corner
            elseif binNum == binNumMax
            push!(adjacentBins , binNum - 1)
            push!(adjacentBins , binNum - N - 1)
            # Middle of right bins
            else
            push!(adjacentBins , binNum - 1)
            push!(adjacentBins , binNum + 1)
            push!(adjacentBins , binNum - N - 1)
            push!(adjacentBins , binNum - N + 1)
            end
        # Top bins
        elseif mod(binNum,N) == 0
            push!(adjacentBins , binNum - 1)
            # Top left corner
            if binNum == N
            push!(adjacentBins , binNum + N - 1)
            push!(adjacentBins , binNum + N)
            # Top right corner
            elseif binNum == binNumMax
            push!(adjacentBins , binNum - N)
            push!(adjacentBins , binNum - N - 1)
            # Middle of Top bins
            else
            push!(adjacentBins , binNum - N - 1)
            push!(adjacentBins , binNum - N)
            push!(adjacentBins , binNum + N - 1)
            push!(adjacentBins , binNum + N)
            end
        # Bottom bins
        elseif mod(binNum,N) == 1
             push!(adjacentBins , binNum + 1)
            # Bottom left corner
            if binNum == 1
            push!(adjacentBins , binNum + N)
            push!(adjacentBins , binNum + N + 1)
            # Bottom right corner
            elseif binNum == (M-1)*N + 1
            push!(adjacentBins , binNum - N)
            push!(adjacentBins , binNum - N + 1)
            # Middle of Bottom bins
            else
            push!(adjacentBins , binNum - N + 1)
            push!(adjacentBins , binNum - N)
            push!(adjacentBins , binNum + N + 1)
            push!(adjacentBins , binNum + N)
            end
        # Middle bins
        else
        push!(adjacentBins , binNum - N - 1)
        push!(adjacentBins , binNum - N)
        push!(adjacentBins , binNum - N + 1)
        push!(adjacentBins , binNum - 1)
        push!(adjacentBins , binNum + 1)
        push!(adjacentBins , binNum + N - 1)
        push!(adjacentBins , binNum + N)
        push!(adjacentBins , binNum + N + 1)
        end
        return adjacentBins
    end
    