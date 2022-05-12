# Script to plot data from SPH particles moving simulation
using DelimitedFiles
function save_vtu(particles,n_save,dirname)
# don't have return value
# specify n_save = 0 for boundary particles

## saving directory
   # check for existence of paraviewfiles/vtu directory. this is the directory where
   # the .vtu files will be stored. if it does not exist create it
    dirname="VTU_Results"
    dirstatus=string(isdir(dirname))
    if dirstatus=="false"
       mkdir(dirname)  #创建一个新的文件夹
    end
    
	cd("/juSPH/src/single-core juSPH/$dirname")
    ## test vtu output (ascii)
    i = n_save
	# specify file name   
    if n_save != 0
        name = string("slosh","$i",".dat")
        namestr = string("PART","$i",".vtu")
        vtu_name = open("$namestr","w")
    else # for boundary particles
        name = string("slosh","$i",".dat")
        namestr = string("BOUND","$i",".vtu")
        vtu_name = open("$namestr","w")
        i = 1
    end
        
    # specify data to store / output
    np = size(particles,1)
    xp=particles[:,1]   # position
    zp=particles[:,2]
    up=particles[:,6]   # velocity
    wp=particles[:,7]
    fp=particles[:,3]   # flag
    rhop=particles[:,5] # density
    P=particles[:,8]    # pressure
    
    # output to file in vtu format
    vtu_name = open("$namestr","a")
    write(vtu_name,"<?xml version=\"1.0\"?>\r\n")
    write(vtu_name,"<VTKFile type= \"UnstructuredGrid\"  version= \"0.1\"  byte_order= \"BigEndian\">\r\n")
    write(vtu_name," <UnstructuredGrid>\r\n")
    write(vtu_name,"  <Piece NumberOfPoints=\"$np\" NumberOfCells=\"$np\">\r\n")
    
    # write in pressure data
    write(vtu_name,"   <PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\r\n")
    write(vtu_name,"   <DataArray type=\"Float64\" Name=\"Pressures\" format=\"ascii\">\r\n")
    for ii=1:np
        writedlm(vtu_name,P[ii],"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")

    # write density data
    write(vtu_name,"    <DataArray type=\"Float64\" Name=\"Density\" format=\"ascii\">\r\n")
    for ii=1:np
        writedlm(vtu_name,rhop[ii],"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")

    # this section is used to color different particles based the input idiv specified above.
    write(vtu_name,"    <DataArray type=\"Float64\" Name=\"Scalarplot\" format=\"ascii\">\r\n")
    for ii=1:np
        writedlm(vtu_name,fp[ii],"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")

    # write velocity data
    write(vtu_name,"    <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\r\n")
    for ii=1:np
        vel=[up[ii] 0 wp[ii]]
        writedlm(vtu_name,vel,"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")
    write(vtu_name,"   </PointData>\r\n")

    # write particle position data
    write(vtu_name,"   <Points>\r\n")
    write(vtu_name,"    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\r\n")
    for ii=1:np
        pos=[xp[ii] 0 zp[ii]]
        writedlm(vtu_name,pos,"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")
    write(vtu_name,"   </Points>\r\n")

    # write cell data. cell is of type vertex.
    write(vtu_name,"   <Cells>\r\n")
    write(vtu_name,"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\r\n")
    for ii=1:np
        writedlm(vtu_name,ii-1,"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")
    write(vtu_name,"\r\n")
    write(vtu_name,"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\r\n")
    for ii=1:np
        writedlm(vtu_name,ii,"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")
    write(vtu_name,"\r\n")
    write(vtu_name,"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\r\n")
    for ii=1:np
        writedlm(vtu_name,1,"\t")
    end
    write(vtu_name,"\r\n")
    write(vtu_name,"    </DataArray>\r\n")
    write(vtu_name,"   </Cells>\r\n")
    write(vtu_name,"  </Piece>\r\n")
    write(vtu_name," </UnstructuredGrid>\r\n")
    write(vtu_name,"</VTKFile>")
    close(vtu_name)
    cd("/juSPH/src/single-core juSPH")
end