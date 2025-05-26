###### custom functions
# Define.Var = A,0.1
# Define.Func = ABC@{[(A*PHI<1>)]}@
# default field variable: "PHI", "dPHI_dt", "lap_PHI", "PHI_X", "dPHI_X_dt", "X", "dX_dt", "T", "dT_dt", "lap_T", "P", "dP_dt", "lap_P", "PHI_P", "dPHI_P_dt", "lap_PHI_P"
# default functions		: "pow", "sqrt", "abs", "exp", "ln", "log", "sin", "cos", "tan", "asin", "acos", "atan"

###### Init Solver
Solver.Mesh.Nx = 64
Solver.Mesh.Ny = 64
Solver.Mesh.Nz = 64

# Solver.SimulationSystem.PhaseNames = (name_0, name_1, ... ) 
Solver.SimulationSystem.PhaseNames = (Grain0)

# Solver.ControlEquation.PhaseField = 0 - Const, 1 - AllenCahn_Pairwise, ... 
Solver.ControlEquation.PhaseField = 0

Solver.ControlEquation.PhaseField.PhaseNumber = 2
# Solver.ControlEquation.PhaseField.Property.PhiName = ( Phi_index_1, Phi_index_2, ... ) 
Solver.ControlEquation.PhaseField.Property.Grain0 = (0,1)
Solver.ControlEquation.PhaseField.VTS.index = true

# .matrix_phi = ( phi_index, phi_value ) 
Preprocess.Microstructure.matrix_phi = (0, 1.0)

# Preprocess.Microstructure.porous = (first_phi_index, second_phi_index, porosity, noise_level) 
Preprocess.Microstructure.porous = (0,1,0.5,0.1)

Solver.Output.VTS.frequence = 1000
Solver.Output.VTS.show_with_boundary = false
