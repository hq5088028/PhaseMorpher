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
Solver.SimulationSystem.PhaseNames = (Grain0, Grain1)

# Solver.ControlEquation.PhaseField = 0 - Const, 1 - AllenCahn_Pairwise, ... 
Solver.ControlEquation.PhaseField = 0

Solver.ControlEquation.PhaseField.PhaseNumber = 2
# Solver.ControlEquation.PhaseField.Property.PhiName = ( Phi_index_1, Phi_index_2, ... ) 
Solver.ControlEquation.PhaseField.Property.Grain0 = (0)
Solver.ControlEquation.PhaseField.Property.Grain1 = (1)
Solver.ControlEquation.PhaseField.VTS.index = true

# .matrix_phi = ( phi_index, phi_value ) 
Preprocess.Microstructure.matrix_phi = (0, 1.0)

Preprocess.Microstructure.geometry_layer_number = 1
# .property = (phi_index, geometry_type, rotation_gauge, reverse_region) 
#              geometry_type  : 0 - None, 1 - Ellipsoid, 2 - Polyhedron, 3 - Geo_SegmentedCylinder, 4 - RectangularCuboid 
#              rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ 
#                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX 
#  The origin of the simulation mesh is ( x = 1 , y = 1 , z = 1 ). 
Preprocess.Microstructure.geometry_layer_0.property = (1,3,1,false)
# .segmented_cylinder = [(radius),(central_axis_point_x,central_axis_point_y,central_axis_point_z), ... at least two point ... ,(rotation_angle_1,rotation_angle_2,rotation_angle_3)] 
Preprocess.Microstructure.geometry_layer_0.segmented_cylinder = \\
                    [(5),\\
                    (32,32,5),(32,5,32),(32,32,58),(32,58,32),(32,32,5),\\
                    (20,0,0)]
Preprocess.Microstructure.geometry_layer_0.phi = 1.000000
Preprocess.Microstructure.geometry_layer_0.is_normalized = true

Solver.Output.VTS.frequence = 1000
Solver.Output.VTS.show_with_boundary = false
