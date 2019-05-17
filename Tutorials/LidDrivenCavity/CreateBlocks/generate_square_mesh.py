import numpy as np

Nblocks = 4
nx, ny, nz  = (129, 129, 2)

x = np.linspace(0,  1,  nx)
y = np.linspace(0,  1,  ny)
z = np.linspace(0,0.1,  nz)

gridX, gridY = np.meshgrid(x,y)

print gridX
print gridY
gridZ1 = np.full((nx*ny), 0.0)
gridZ2 = np.full((nx*ny), 0.1)
print gridZ1, gridZ2
gridZ = np.stack((gridZ1, gridZ2)).flatten()
print gridZ.size
gridX = gridX.flatten()
gridY = gridY.flatten()
gridX = np.stack((gridX, gridX)).flatten()
gridY = np.stack((gridY, gridY)).flatten()
print gridX
print gridY
print gridZ.flatten()


## write grid file
def write_dat_file(filename):
   print 'writing ', filename
   header="  ".join([str(nx), str(ny), str(nz)])
   np.savetxt(filename, np.c_[gridX, gridY, gridZ], fmt='%26.16e %26.16e %26.16e', delimiter="  ",header=header, comments='')

write_dat_file("grid_00.txt")
