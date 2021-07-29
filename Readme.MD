First of all, *qeijo* is not an acronym. Okay, we have the Q and E for Quantum Espresso, for sure, but that's it. 
The point was just to use a funny name for the package. 

*qeijo* is an incorrect way - on purpose, of course - of spelling "queijo", which is simply the Portuguese word 
for "cheese".

This package contains a lightweight library that allows the user to launch electronic structure calculations with 
Quantum Espresso and then collect the results, all done in a Pythonic way. Just to be clear, I don't 
claim that this package can compete in terms of features with much larger frameworks for atomistic 
simulations, like ASE. All you can expect from *qeijo* is that it is light and easy to use.

In version 0.1, there is only one module, *pw*, which must be used to build an input to the pw.x program, 
launch a total energy calculation and get from it important quantities, like the total energy of the 
system, the relaxed atomic coordinates, magnetization, etc. In future releases, I plan to add interfaces 
to other programs from the Quantum Espresso package, as well as work a bit more on user documentation.

### Basic Usage Example

---

    # Launching a calculation with pw.x of H2 molecule, reading the input from an existing input file
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

