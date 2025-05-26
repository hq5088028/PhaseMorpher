###### whether to test this input file
InputFile.debug = FALSE

###### Init Solver
Solver.Loop.begin_step = 0
Solver.Loop.end_step = 1000
Solver.Mesh.Nx = 100
Solver.Mesh.Ny = 1
Solver.Mesh.Nz = 1
Solver.Loop.dt = 0.1
# Solver.Mesh.BoundaryCondition : 0 - FIXED , 1 - PERIODIC , 2 - ADIABATIC
Solver.Mesh.BoundaryCondition.x_up = 2
Solver.Mesh.BoundaryCondition.x_down = 2
Solver.Mesh.BoundaryCondition.y_up = 1
Solver.Mesh.BoundaryCondition.y_down = 1
Solver.Mesh.BoundaryCondition.z_up = 1
Solver.Mesh.BoundaryCondition.z_dowm = 1

# Solver.SimulationSystem.PhaseNames = (name_0, name_1, ... ) 
Solver.SimulationSystem.PhaseNames = (Grain0)

Solver.ControlEquation.PhaseField.PhaseNumber = 2
# Solver.ControlEquation.PhaseField.Property.PhiName = ( Phi_index_1, Phi_index_2, ... ) 
Solver.ControlEquation.PhaseField.Property.Grain0 = (0,1)

# Solver.ControlEquation.PhaseField = 0 - Const, 1 - AllenCahn_Pairwise, ... 
Solver.ControlEquation.PhaseField = 1
Solver.ControlEquation.Pairwise.Acceleration.ContainerSize = 2
# Solver.ControlEquation.Pairwise.Lij.const  = Lij_value 
#                                .matrix = [(phi_i, phi_j, Lij_value), ... ] 
#                                .block = [(phi_begin, phi_end, Lij_value), ... ] 
Solver.ControlEquation.Pairwise.Lij.const = 1.0
# Solver.ControlEquation.Pairwise.InterfaceEnergy.int_gradient : 0 - Steinbach_1996 , 1 - Steinbach_1999 , 2 - Steinbach_G2009 , 3 - DanielCrack_G2016
Solver.ControlEquation.Pairwise.InterfaceEnergy.int_gradient = 2
# Solver.ControlEquation.Pairwise.InterfaceEnergy.int_potential : 0 - Nestler_Well , 1 - Nestler_Obstacle , 2 - Steinbach_P2009 , 3 - DanielCrack_P2016
Solver.ControlEquation.Pairwise.InterfaceEnergy.int_potential = 2
Solver.ControlEquation.Pairwise.interface_width = 10.0
# Solver.ControlEquation.Pairwise.xi_ab.const  = xi_ab 
#                                      .matrix = [(phi_a_name, phi_b_name, xi_ab_value), ...] 
Solver.ControlEquation.Pairwise.xi_ab.const = 1.0
# Solver.ControlEquation.Pairwise.xi_abc.const  = xi_ab 
#                                       .matrix = [(phi_a_name, phi_b_name, phi_c_name, xi_abc_value), ...] 
Solver.ControlEquation.Pairwise.xi_abc.const = 0.0

# .matrix_phi = ( phi_index, phi_value ) 
Preprocess.Microstructure.matrix_phi = (0, 1.0)

Preprocess.Microstructure.geometry_layer_number = 1
# .property = (phi_index, geometry_type, rotation_gauge, reverse_region) 
#              geometry_type  : 0 - None, 1 - Ellipsoid, 2 - Polyhedron, 3 - Geo_SegmentedCylinder, 4 - RectangularCuboid 
#              rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ 
#                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX 
#  The origin of the simulation mesh is ( x = 1 , y = 1 , z = 1 ). 
Preprocess.Microstructure.geometry_layer_0.property = (1,2,1,false)
Preprocess.Microstructure.geometry_layer_0.polyhedron =  {[(75,0,0)],[(50.5,0,0),(50.5,1,0),(50.5,0,1)],[(0,0,0)]}
Preprocess.Microstructure.geometry_layer_0.phi = 1.0
Preprocess.Microstructure.geometry_layer_0.is_normalized = true

Solver.Output.LOG.loop_info_step = 10
Solver.Output.VTS.frequence = 10
Solver.Output.VTS.show_with_boundary = false
Solver.ControlEquation.PhaseField.VTS.all_phi = true





