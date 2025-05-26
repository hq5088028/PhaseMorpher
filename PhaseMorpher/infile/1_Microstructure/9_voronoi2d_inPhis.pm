###### custom functions
# Define.Var = A,0.1
# Define.Func = ABC@{[(A*PHI<1>)]}@
# default field variable: "PHI", "dPHI_dt", "lap_PHI", "PHI_X", "dPHI_X_dt", "X", "dX_dt", "T", "dT_dt", "lap_T", "P", "dP_dt", "lap_P", "PHI_P", "dPHI_P_dt", "lap_PHI_P"
# default functions		: "pow", "sqrt", "abs", "exp", "ln", "log", "sin", "cos", "tan", "asin", "acos", "atan"

###### Init Solver
Solver.Mesh.Nx = 100
Solver.Mesh.Ny = 100
Solver.Mesh.Nz = 1

# Solver.SimulationSystem.PhaseNames = (name_0, name_1, ... ) 
Solver.SimulationSystem.PhaseNames = (Grain0)

# Solver.ControlEquation.PhaseField = 0 - Const, 1 - AllenCahn_Pairwise, ... 
Solver.ControlEquation.PhaseField = 0

Solver.ControlEquation.PhaseField.PhaseNumber = 20
# Solver.ControlEquation.PhaseField.Property.PhiName = ( Phi_index_1, Phi_index_2, ... ) 
Solver.ControlEquation.PhaseField.Property.Grain0 = ($TUBE[0,19,1]$)
Solver.ControlEquation.PhaseField.VTS.index = true

# .matrix_phi = ( phi_index, phi_value ) 
Preprocess.Microstructure.matrix_phi = (0, 1.0)

Preprocess.Microstructure.geometry_layer_number = 1
# .property = (phi_index, geometry_type, rotation_gauge, reverse_region) 
#              geometry_type  : 0 - None, 1 - Ellipsoid, 2 - Polyhedron, 3 - Geo_SegmentedCylinder, 4 - RectangularCuboid 
#              rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ 
#                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX 
#  The origin of the simulation mesh is ( x = 1 , y = 1 , z = 1 ). 
Preprocess.Microstructure.geometry_layer_0.property = (1,1,1,false)
# .ellipsoid = [(core_x,core_y,core_z),(radius_x,radius_y,radius_z),(rotation_angle_1,rotation_angle_2,rotation_angle_3)] 
Preprocess.Microstructure.geometry_layer_0.ellipsoid = [(51,51,1),(45,30,0),(0,10,0)]
Preprocess.Microstructure.geometry_layer_0.phi = 1.000000
Preprocess.Microstructure.geometry_layer_0.is_normalized = true

# .voronoi = [(phi_index_begin, phi_index_end), ( box_origin_point ),( box_end_point )] 
Preprocess.Microstructure.voronoi = [(1,19),(1,1,1),(100,100,1)]
Preprocess.Microstructure.Voronoi.rand_seed = 10
Preprocess.Microstructure.Voronoi.const_distance = 13
Preprocess.Microstructure.Voronoi.VTS.debug = true
# .in_phi_indexs = (in phi_index1, ... ) 
Preprocess.Microstructure.Voronoi.in_phi_indexs = (1)
Preprocess.Microstructure.Voronoi.is_box_periodic = false

Solver.Output.VTS.frequence = 1000
Solver.Output.VTS.show_with_boundary = false
