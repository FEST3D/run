# raetaf.x
The grid for RAE2822 test case in Plot3D (.x file) format is
located in current folder. The purpose of the blocking_point.f90
code is to decompose grid into mulitple blocks and write separate
grid file for each block in particular format required by the FEST-3D
solver.

# Blocking_point.f90
In the blocking_point.f90, you can change the number of block you
want based on the requirement. It a better to keep the number
a multiple of 2, which reduces the chance of this code to fail.


# Create block grid.
In order to create block grid just run <br>
```$make```
and the block grid will be placed in the grid/ folder
