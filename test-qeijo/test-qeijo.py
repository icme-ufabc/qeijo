# Launching a calculation with pw.x of H2 molecule, reading the input from 
# an existing input file

from qeijo import pw

cmd="mpirun -np 2 pw.x"
h2_calc=pw.calc()
h2_calc.read_input("h2.inp")
inpstr=h2_calc.build_input()
h2_out=h2_calc.run(command_line=cmd,input_string=inpstr,saveout=True,
                   outfile="h2.out",savecoords=True,coordfile="h2.xyz")

# Print the total energy in eV
print("Total energy: %f" % h2_out.energy[-1])

# Print the relaxed coordinates of the the two H atoms
print("H1: %f %f %f" % (h2_out.x[0],h2_out.y[0],h2_out.z[0]))
print("H2: %f %f %f" % (h2_out.x[1],h2_out.y[1],h2_out.z[1]))
