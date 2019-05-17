all: generate_square_mesh.py blocking_point.f90
	rm -rf grid/
	python generate_square_mesh.py
	gfortran -fdefault-real-8 blocking_point.f90
	./a.out
	rm a.out
