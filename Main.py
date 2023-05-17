from Initialize import Initialize

print ("\n"*80)

Lattice = Initialize()
print("Simulation Started...:")
Lattice.Simulate()
Lattice.Output(3)
print("Simulation Finished:")