include("get_adjacentBins.jl")
function link_list(particles,xMin,yMin,Nx,Ny,dx_x,dy_y)

# Nearest neighbor particle search(link list search method) 
# Calculate the binNum of each particle
    
    # establish Nx*Ny bins
        # Calculate the bin arrangement
        numBins = Nx*Ny
        numBins = convert(Int,numBins)
  
        # establish adjacentBins of every bin
        adjacentBins = SharedArray{Int}(8,numBins)
        @sync @distributed for i = 1:numBins
            local AdjacentBins = get_adjacentBins(i,numBins,Nx,Ny)
            @inbounds for j = 1:length(AdjacentBins)
                adjacentBins[j,i] = AdjacentBins[j]   
            end
        end

        # Calculate the binNum of each particle
        function number(particles,k)
            binNum = (cld((particles[k,1]-xMin),dx_x) - 1) * Ny + cld((particles[k,2]-yMin),dy_y)
            binNum = convert(Int,binNum)
          return binNum
        end
       
        Num = zeros(Int64,size(particles,1))
        Num = SharedArray(Num)
        @inbounds @sync @distributed for k = 1 : size(particles,1)
           num = number(particles,k)
           Num[k]=num
        end

        # Store the end index position of the number of particles under each bin number
        Gridnum = SharedArray{Int64}(numBins,1)
        SortNum = sort(Num)
        @sync @distributed for i = 1:numBins
            Gridnum[i] = searchsortedlast(SortNum,i)
        end

        # Store the correct number of particles under each bin number
        ParticlesList = SharedArray{Int64}(size(particles,1))
        @sync @distributed for i = 1:numBins
            local List = findall(in(i),Num)
            local j = 0
            if i==1
                @inbounds for k = 1:Gridnum[i]
                    j = j+1
                    if Gridnum[i]!=0
                        ParticlesList[k] = List[j]
                    end
                end
            else
                @inbounds for k = Gridnum[i-1]+1:Gridnum[i]
                    j = j+1
                    if Gridnum[i-1] != Gridnum[i]
                        ParticlesList[k] = List[j] 
                    end
                end
            end            
        end
                    
   return numBins,adjacentBins,Gridnum,ParticlesList

end