import numpy as np
import matplotlib.pyplot as plt


# variables that can be changed  
imax = 97
jmax = 49
kmax = 2
blocks = 2

# fixed variables
x_min = -1.5
x_max = 1.5
y_min = 0.0
y_max = 0.8
z_min = 0.0
z_max = 0.1

# number of points in one block in x direction 
nd = (imax-1)//2 + 1

# output array
points = np.empty( [blocks,nd,jmax,kmax,3] )

def bump(x):
  return 0.06250*np.exp(1)**(-25*(x**2))

def gen_grid():

  x_lower = np.zeros((blocks,nd))
  xmin = x_min
  xmax = x_max
  delx = (x_max - x_min)/blocks
  for each in range(blocks):
    xmin = x_min + delx*each
    xmax = xmin + delx
    x_lower[each,:] = np.linspace(xmin, xmax, nd, endpoint=True)
  y_lower = bump(x_lower)
  y_upper = y_max
  #print x_lower
  #print y_lower

  k = 0
  for each in range(blocks):
    for i in range(nd):
      for j in range(jmax):
        points[each,i,j,0,0] = x_lower[each,i]
        points[each,i,j,0,1] = y_lower[each,i] + (j)*(y_upper - y_lower[each,i])/(jmax-1.)
    for k in range(kmax):
      points[each,:,:,k,2] = z_min + k*(z_max-z_min)/(kmax-1.)
      points[each,:,:,k,0] = points[each,:,:,0,0]
      points[each,:,:,k,1] = points[each,:,:,0,1]
  return

def write_grid(format, block):

#  f = open('grid_bump_'+str(imax)+'x'+str(jmax)+'x'+str(kmax)+'_'+str(block)+'.dat', 'w')
  f = open('grid/grid_'+str(block).zfill(2)+'.txt', 'w')

  if (format=='tecplot'):
    head = 'variables= x y z'+'\n'+\
    'zone T="Bump_grid", i='+str(nd)+', j='+str(jmax)+ ', k=' + str(kmax)+'\n'
    f.write(head)
#    for line in points:
#      for item in line:
#        np.savetxt(f, item, fmt='%.16e', delimiter='\t')
    for k in range(kmax):
      for j in range(jmax):
        for i in range(nd):
          s = "{:.16f}".format(points[block,i,j,k,0])+'\t'+"{:.16f}".format(points[block,i,j,k,1])+'\t'+"{:.16f}".format(points[block,i,j,k,2])+'\n'
          f.write(s)

  if (format=='solver'):
    head = str(nd) + ' ' + str(jmax) + ' ' + str(kmax) + '\n'
    f.write(head)
    for k in range(kmax):
      for j in range(jmax):
        for i in range(nd):
          s = "{:.18e}".format(points[block,i,j,k,0])+'\t'+"{:.18e}".format(points[block,i,j,k,1])+'\t'+"{:.18e}".format(points[block,i,j,k,2])+'\n'
          f.write(s)
  
  f.close()

def main():
  gen_grid()
  for each in range(blocks):
    write_grid('solver', each)

main()

