using HDF5

function unlink(filename)
    assert(length(filename)>0)
    r=ccall((:unlink, "libc"), Int32, (Ptr{Uint8},), filename)
    if r!=0
        n=errno()
        error("errno: ", n, " ", strerror(n))
    end
end


dataname="z.h5"
if isfile(dataname)
    unlink(dataname)
end

f=HDF5.h5open(dataname,"w")
ds=HDF5.dataspace([3])
traj_type=HDF5.t_create(HDF5.H5T_COMPOUND, 3)

s="s"
ps=convert(Ptr{Uint8}, s)
HDF5.h5t_insert(traj_type, ps, 0, HDF5.H5T_STD_I64LE)
HDF5.close(f)
