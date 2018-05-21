input: clean
	gmsh -3 ./src/pile_cuboid.geo
	cp src/pile_cuboid.msh .

clean: 
	rm -f src/*.msh *.msh *.halo *.vtu *.pvtu *.s  *.stat *.log-0 *.err-0 matrixdump*