include("get_adjacentBins.jl")
function link_list(particles,xMin,yMin,Nx,Ny,dx_x,dy_y)

# Nearest neighbor particle search(link list search method) 
# Calculate the binNum of each particle
    
    # establish Nx*Ny bins
        # Calculate the bin arrangement
        numBins = Nx*Ny
        numBins = convert(Int,numBins)
  
        # establish adjacentBins of every bin
        adjacentBins = zeros(Int,8,numBins)
        for i = 1:numBins
            local AdjacentBins = get_adjacentBins(i,numBins,Nx,Ny)
            for j = 1:length(AdjacentBins)
                adjacentBins[j,i] = AdjacentBins[j]   
            end
        end

        # Calculate the binNum of each particle
        function number(particles,k)
            binNum = (cld((particles[k,1]-xMin),dx_x) - 1) * Ny + cld((particles[k,2]-yMin),dy_y)
            binNum = convert(Int,binNum)
          return binNum
        end
       
        Num = zeros(Int,size(particles,1))      
        for k = 1 : size(particles,1)
           num = number(particles,k)
           Num[k]=num
        end

        # Store the end index position of the number of particles under each bin number
        Gridnum = zeros(Int,numBins) 
        SortNum = sort(Num)
        for i = 1:numBins
            Gridnum[i] = searchsortedlast(SortNum,i)
        end

        # Store the correct number of particles under each bin number
        ParticlesList = Int64[] 
        for i = 1:numBins
            append!(ParticlesList,findall(in(i),Num))
        end
              
   return numBins,adjacentBins,Gridnum,ParticlesList

end