#pragma once
#include "../postprocess_modules/WriteMeshData.h"
#include "../../tools/mathTools.h"
#include "../input_modules/ioFiles_Params.h"
#include "../input_modules/inputfiles/selectFile.h"
#include "../model_modules/Model_Params.h"
#include "../model_modules/PhaseField/PhaseField_Params.h"
#include "../model_modules/ConcentrationField/ConcentrationField_Params.h"
#include "../model_modules/TemperatureField/TemperatureField_Params.h"
#include "../Module.h"
namespace pf {
	enum Geometry { Geo_None, Geo_Ellipsoid, Geo_Polyhedron, Geo_SegmentedCylinder, Geo_RectangularCuboid };
	enum NucleationPosition { NP_Anywhere, NP_Bulk, NP_Interface };
	enum NucleationConcentration { NC_Default, NC_Defined };
	struct Point {
		REAL x;
		REAL y;
		REAL z;
		Point() {
			x = 0;
			y = 0;
			z = 0;
		}
		Point(REAL _x, REAL _y, REAL _z) {
			x = _x;
			y = _y;
			z = _z;
		}
		Point(int _x, int _y, int _z) {
			x = REAL(_x);
			y = REAL(_y);
			z = REAL(_z);
		}
		void set(REAL _x, REAL _y, REAL _z) {
			x = _x;
			y = _y;
			z = _z;
		}
		REAL to_length(REAL _x, REAL _y, REAL _z) {
			return std::sqrt((_x - x) * (_x - x) + (_y - y) * (_y - y) + (_z - z) * (_z - z));
		}
		REAL to_length(Point p) {
			return std::sqrt((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) + (p.z - z) * (p.z - z));
		}
		REAL to_cosine(Point direction1, Point direction2) {
			Vector3 vec1(direction1.x - x, direction1.y - y, direction1.z - z), vec2(direction2.x - x, direction2.y - y, direction2.z - z);
			return vec1 * vec2 / vec1.length() / vec2.length();
		}
		void do_boundary(BoundaryCondition x_up_bc, BoundaryCondition y_up_bc, BoundaryCondition z_up_bc, BoundaryCondition x_down_bc, BoundaryCondition y_down_bc, BoundaryCondition z_down_bc, REAL x_limit, REAL y_limit, REAL z_limit) {
			if (x_down_bc == BoundaryCondition::ADIABATIC && x < 0)
				x = 0;
			else if (x_up_bc == BoundaryCondition::ADIABATIC && x >= x_limit)
				x = x_limit - 1;
			else if (x_down_bc == BoundaryCondition::FIXED && x < 0)
				x = 0;
			else if (x_up_bc == BoundaryCondition::FIXED && x >= x_limit)
				x = x_limit - 1;
			else if (x_down_bc == BoundaryCondition::PERIODIC && x < 0)
				x += x_limit;
			else if (x_up_bc == BoundaryCondition::PERIODIC && x >= x_limit)
				x -= x_limit;
			else if (y_down_bc == BoundaryCondition::ADIABATIC && y < 0)
				y = 0;
			else if (y_up_bc == BoundaryCondition::ADIABATIC && y >= y_limit)
				y = y_limit - 1;
			else if (y_down_bc == BoundaryCondition::FIXED && y < 0)
				y = 0;
			else if (y_up_bc == BoundaryCondition::FIXED && y >= y_limit)
				y = y_limit - 1;
			else if (y_down_bc == BoundaryCondition::PERIODIC && y < 0)
				y += y_limit;
			else if (y_up_bc == BoundaryCondition::PERIODIC && y >= y_limit)
				y -= y_limit;
			else if (z_down_bc == BoundaryCondition::ADIABATIC && z < 0)
				z = 0;
			else if (z_up_bc == BoundaryCondition::ADIABATIC && z >= z_limit)
				z = z_limit - 1;
			else if (z_down_bc == BoundaryCondition::FIXED && z < 0)
				z = 0;
			else if (z_up_bc == BoundaryCondition::FIXED && z >= z_limit)
				z = z_limit - 1;
			else if (z_down_bc == BoundaryCondition::PERIODIC && z < 0)
				z += z_limit;
			else if (z_up_bc == BoundaryCondition::PERIODIC && z >= z_limit)
				z -= z_limit;
			else
				return;
			return do_boundary(x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc, x_limit, y_limit, z_limit);
		}
		Point& operator=(const Point& n) {
			this->x = n.x;
			this->y = n.y;
			this->z = n.z;
			return *this;
		}
		Point operator+(const Point& n) {
			Point re;
			re.x = this->x + n.x;
			re.y = this->y + n.y;
			re.z = this->z + n.z;
			return re;
		}
		Point operator-(const Point& n) {
			Point re;
			re.x = this->x - n.x;
			re.y = this->y - n.y;
			re.z = this->z - n.z;
			return re;
		}
		Point operator*(const REAL n) {
			Point re;
			re.x = this->x * n;
			re.y = this->y * n;
			re.z = this->z * n;
			return re;
		}
		Point operator/(const REAL n) {
			Point re;
			re.x = this->x / n;
			re.y = this->y / n;
			re.z = this->z / n;
			return re;
		}
		Point& operator+=(const Point& n) {
			this->x += n.x;
			this->y += n.y;
			this->z += n.z;
			return *this;
		}
		Point& operator-=(const Point& n) {
			this->x -= n.x;
			this->y -= n.y;
			this->z -= n.z;
			return *this;
		}
		Point& operator*=(const REAL n) {
			this->x *= n;
			this->y *= n;
			this->z *= n;
			return *this;
		}
		Point& operator/=(const REAL n) {
			this->x /= n;
			this->y /= n;
			this->z /= n;
			return *this;
		}
		Point operator+(const Vector3& n) {
			Point re;
			re.x = this->x + n[0];
			re.y = this->y + n[1];
			re.z = this->z + n[2];
			return re;
		}
		Point operator-(const Vector3& n) {
			Point re;
			re.x = this->x - n[0];
			re.y = this->y - n[1];
			re.z = this->z - n[2];
			return re;
		}
		Point& operator+=(const Vector3& n) {
			this->x += n[0];
			this->y += n[1];
			this->z += n[2];
			return *this;
		}
		Point& operator-=(const Vector3& n) {
			this->x -= n[0];
			this->y -= n[1];
			this->z -= n[2];
			return *this;
		}
	};
	inline Vector3 get_vector(Point Tail, Point Head) {
		return Vector3(Head.x - Tail.x, Head.y - Tail.y, Head.z - Tail.z);
	}
	struct surf_func_3D {
		// temp: a * x + b * y + c * z + d = 0
		surf_func_3D() {
			a = 0;
			b = 0;
			c = 0;
			d = 0;
		}
		surf_func_3D(Point p1, Point p2, Point p3) {
			init(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
		}
		surf_func_3D(Vector3 normal, Point point) {
			a = normal[0];
			b = normal[1];
			c = normal[2];
			d = -(a * point.x + b * point.y + c * point.z);
		}
		void init(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2, REAL x3, REAL y3, REAL z3) {
			Vector3 t1(x2 - x1, y2 - y1, z2 - z1);
			Vector3 t2(x3 - x1, y3 - y1, z3 - z1);
			t1 = t1.cross(t2);
			a = t1[0];
			b = t1[1];
			c = t1[2];
			d = -(a * x1 + b * y1 + c * z1);
		}
		void move(REAL x, REAL y, REAL z) {
			d += -(a * x + b * y + c * z);
		}
		REAL distance(REAL _x, REAL _y, REAL _z) {
			return std::abs(a * _x + b * _y + c * _z + d) / std::sqrt(a * a + b * b + c * c);
		}
		bool is_point_on_surf(REAL x, REAL y, REAL z) {
			REAL val = a * x + b * y + c * z + d;
			if (isTwoNumEquality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_on_surf(Point p) {
			REAL val = a * p.x + b * p.y + c * p.z + d;
			if (isTwoNumEquality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_above_surf(REAL x, REAL y, REAL z) {
			if ((a * x + b * y + c * z + d) > 0.0)
				return true;
			else
				return false;
		}
		bool is_point_above_surf(Point p) {
			if ((a * p.x + b * p.y + c * p.z + d) > 0.0)
				return true;
			else
				return false;
		}
		surf_func_3D& operator=(const surf_func_3D& n) {
			a = n.a;
			b = n.b;
			c = n.c;
			d = n.d;
			return *this;
		}
		REAL a, b, c, d;
	};
	struct Ellipsoid {
		//> (x - core.x) * (x - core.x) / radius_x / radius_x + (y - core.y) * (y - core.y) / radius_y / radius_y + (z - core.z) * (z - core.z) / radius_z / radius_z = 1.0
		Ellipsoid() {
			radius_x = 0;
			radius_y = 0;
			radius_z = 0;
			radian_x = 0;
			radian_y = 0;
			radian_z = 0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		void set_core(REAL x, REAL y, REAL z) {
			core.x = x;
			core.y = y;
			core.z = z;
		}
		void set_radius(REAL x_radius, REAL y_radius, REAL z_radius) {
			radius_x = abs(x_radius) + SYS_EPSILON;
			radius_y = abs(y_radius) + SYS_EPSILON;
			radius_z = abs(z_radius) + SYS_EPSILON;
		}
		void move(REAL x, REAL y, REAL z) {
			core.x += x;
			core.y += y;
			core.z += z;
		}
		void add_rotation_radian(REAL radian1, REAL radian2, REAL radian3) {
			radian_x += -radian1;
			radian_y += -radian2;
			radian_z += -radian3;
		}
		void set_rotation_radian_and_rotation_gauge(REAL _radian[3], RotationGauge _rotationGauge) {
			radian_x = -_radian[0];
			radian_y = -_radian[1];
			radian_z = -_radian[2];
			rotationGauge = _rotationGauge;
		}
		bool check_point_inside_ellipsoid(REAL x, REAL y, REAL z) {
			REAL test = (x - core.x) * (x - core.x) / radius_x / radius_x + (y - core.y) * (y - core.y) / radius_y / radius_y + (z - core.z) * (z - core.z) / radius_z / radius_z;
			if (test <= 1.0)
				return true;
			else
				return false;
		}
		bool check_point_inside_ellipsoid(Point p) {
			REAL test = (p.x - core.x) * (p.x - core.x) / radius_x / radius_x + (p.y - core.y) * (p.y - core.y) / radius_y / radius_y + (p.z - core.z) * (p.z - core.z) / radius_z / radius_z;
			if (test <= 1.0)
				return true;
			else
				return false;
		}
		Ellipsoid& operator=(const Ellipsoid& n) {
			this->core = n.core;
			this->radius_x = n.radius_x;
			this->radius_y = n.radius_y;
			this->radius_z = n.radius_z;
			this->radian_x = n.radian_x;
			this->radian_y = n.radian_y;
			this->radian_z = n.radian_z;
			rotationGauge = n.rotationGauge;
			return *this;
		}
		Point core;
		REAL radius_x;
		REAL radius_y;
		REAL radius_z;
		REAL radian_x;
		REAL radian_y;
		REAL radian_z;
		RotationGauge rotationGauge;
	};
	struct SegmentedCylinder {
		SegmentedCylinder() {
			radius = 0;
			geometric_center = Point(0, 0, 0);
			radian_x = 0;
			radian_y = 0;
			radian_z = 0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		void add_point(Point p) {
			central_axis_line.push_back(p);
			geometric_center = Point(0, 0, 0);
			for (int index = 0; index < central_axis_line.size(); index++)
				geometric_center += central_axis_line[index];
			geometric_center /= REAL(central_axis_line.size());
		}
		void set_radius(REAL _radius) {
			radius = abs(_radius) + SYS_EPSILON;
		}
		void move(REAL x, REAL y, REAL z) {
			geometric_center += Point(x, y, z);
			for (int index = 0; index < central_axis_line.size(); index++)
				central_axis_line[index] += Point(x, y, z);
		}
		void add_rotation_radian(REAL radian1, REAL radian2, REAL radian3) {
			radian_x += -radian1;
			radian_y += -radian2;
			radian_z += -radian3;
		}
		void set_rotation_radian_and_rotation_gauge(REAL _radian[3], RotationGauge _rotationGauge) {
			radian_x = -_radian[0];
			radian_y = -_radian[1];
			radian_z = -_radian[2];
			rotationGauge = _rotationGauge;
		}
		bool check_point_inside_segmented_cylinder(REAL x, REAL y, REAL z) {
			Point p(x, y, z);
			for (int index = 1; index < central_axis_line.size(); index++) {
				REAL cosin1 = central_axis_line[index - 1].to_cosine(p, central_axis_line[index]),
					cosin2 = central_axis_line[index].to_cosine(p, central_axis_line[index - 1]);
				if (cosin1 > 0.0 && cosin2 > 0.0) {
					REAL length = central_axis_line[index].to_length(p), d2 = cosin2 * length,
						p2l_length_2 = length * length - d2 * d2;
					if (p2l_length_2 <= radius * radius)
						return true;
				}
				if (index < central_axis_line.size() - 1) {
					REAL length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
						+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
						+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
					if (length2 <= radius * radius)
						return true;
				}
				if (index == central_axis_line.size() - 1) {
					if (isTwoNumEquality(central_axis_line[0].x, central_axis_line[index].x) 
						&& isTwoNumEquality(central_axis_line[0].y, central_axis_line[index].y) 
						&& isTwoNumEquality(central_axis_line[0].z, central_axis_line[index].z)) {
						REAL length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
							+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
							+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
						if (length2 <= radius * radius)
							return true;
					}
				}
			}
			return false;
		}
		bool check_point_inside_segmented_cylinder(Point p) {
			for (int index = 1; index < central_axis_line.size(); index++) {
				REAL cosin1 = central_axis_line[index - 1].to_cosine(p, central_axis_line[index]),
					cosin2 = central_axis_line[index].to_cosine(p, central_axis_line[index - 1]);
				if (cosin1 > 0.0 && cosin2 > 0.0) {
					REAL length = central_axis_line[index].to_length(p), d2 = cosin2 * length,
						p2l_length_2 = length * length - d2 * d2;
					if (p2l_length_2 < radius * radius)
						return true;
				}
				if (index < central_axis_line.size() - 1) {
					REAL length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
						+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
						+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
					if (length2 <= radius * radius)
						return true;
				}
				if (index == central_axis_line.size() - 1) {
					if (isTwoNumEquality(central_axis_line[0].x, central_axis_line[index].x) 
						&& isTwoNumEquality(central_axis_line[0].y, central_axis_line[index].y) 
						&& isTwoNumEquality(central_axis_line[0].z, central_axis_line[index].z)) {
						REAL length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
							+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
							+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
						if (length2 <= radius * radius)
							return true;
					}
				}
			}
			return false;
		}
		SegmentedCylinder& operator=(const SegmentedCylinder& n) {
			this->central_axis_line = n.central_axis_line;
			this->radius = n.radius;
			this->geometric_center = n.geometric_center;
			this->radian_x = n.radian_x;
			this->radian_y = n.radian_y;
			this->radian_z = n.radian_z;
			rotationGauge = n.rotationGauge;
			return *this;
		}
		vector<Point> central_axis_line;
		REAL radius;
		Point geometric_center;
		REAL radian_x;
		REAL radian_y;
		REAL radian_z;
		RotationGauge rotationGauge;
	};
	struct Polyhedron {
		Polyhedron(REAL inside_x = 0, REAL inside_y = 0, REAL inside_z = 0) {
			point_inside_polyhedron.set(inside_x, inside_y, inside_z);
			radian_x = 0;
			radian_y = 0;
			radian_z = 0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		Polyhedron(Point inside_point) {
			point_inside_polyhedron.set(inside_point.x, inside_point.y, inside_point.z);
			radian_x = 0;
			radian_y = 0;
			radian_z = 0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		void move(REAL x, REAL y, REAL z) {
			point_inside_polyhedron += Point(x, y, z);
			for (int index = 0; index < surfaces.size(); index++)
				surfaces[index].move(x, y, z);
		}
		void add_rotation_radian(REAL radian1, REAL radian2, REAL radian3) {
			radian_x += -radian1;
			radian_y += -radian2;
			radian_z += -radian3;
		}
		void set_rotation_radian_and_rotation_gauge(REAL _radian[3], RotationGauge _rotationGauge) {
			radian_x = -_radian[0];
			radian_y = -_radian[1];
			radian_z = -_radian[2];
			rotationGauge = _rotationGauge;
		}
		void add_surf(Point p1, Point p2, Point p3) {
			surf_func_3D surf(p1, p2, p3);
			surfaces.push_back(surf);
		}
		void add_surf(Vector3 norm, Point p) {
			surf_func_3D surf(norm, p);
			surfaces.push_back(surf);
		}
		void set_a_point_inside_polyhedron(REAL x, REAL y, REAL z) {
			point_inside_polyhedron.set(x, y, z);
		}
		void set_a_point_inside_polyhedron(Point p) {
			point_inside_polyhedron.set(p.x, p.y, p.z);
		}
		bool check_point_inside_polyhedron(REAL x, REAL y, REAL z) {
			bool is_inside = true;
			for (auto surf = surfaces.begin(); surf < surfaces.end(); surf++)
				if (surf->is_point_above_surf(x, y, z) != surf->is_point_above_surf(point_inside_polyhedron) && !surf->is_point_on_surf(x, y, z))
					is_inside = false;
			return is_inside;
		}
		bool check_point_inside_polyhedron(Point p) {
			bool is_inside = true;
			for (auto surf = surfaces.begin(); surf < surfaces.end(); surf++)
				if (surf->is_point_above_surf(p) != surf->is_point_above_surf(point_inside_polyhedron) && !surf->is_point_on_surf(p))
					is_inside = false;
			return is_inside;
		}
		Polyhedron& operator=(const Polyhedron& n) {
			this->surfaces = n.surfaces;
			this->point_inside_polyhedron = n.point_inside_polyhedron;
			this->radian_x = n.radian_x;
			this->radian_y = n.radian_y;
			this->radian_z = n.radian_z;
			rotationGauge = n.rotationGauge;
			return *this;
		}
		vector<surf_func_3D> surfaces;
		Point point_inside_polyhedron;
		REAL radian_x;
		REAL radian_y;
		REAL radian_z;
		RotationGauge rotationGauge;
	};
	struct GeometricRegion {
		Geometry geometryProperty;
		Ellipsoid ellipSolid;
		Polyhedron polyhedron;
		SegmentedCylinder cylinder;
		int generate_step;
		int phaseIndex;
		vector<REAL> con;
		REAL temperature;
		REAL phi;
		bool isReverseRegion;
		bool isNormalized;

		GeometricRegion(pf::Geometry _geometry_property = Geo_None, int _generate_step = 0, int _phase_index = 0, bool _isReverseRegion = false) {
			init(_geometry_property, _generate_step, _phase_index, _isReverseRegion);
			temperature = 0;
			phi = 0;
			isNormalized = false;
		}
		void init(pf::Geometry _geometry_property = Geo_None, int _generate_step = 0, int _phase_index = 0, bool _isReverseRegion = false) {
			geometryProperty = _geometry_property;
			generate_step = _generate_step;
			phaseIndex = _phase_index;
			isReverseRegion = _isReverseRegion;
		}
		GeometricRegion& operator=(const GeometricRegion& n) {
			this->generate_step = n.generate_step;
			this->geometryProperty = n.geometryProperty;
			this->temperature = n.temperature;
			this->phaseIndex = n.phaseIndex;
			this->ellipSolid = n.ellipSolid;
			this->polyhedron = n.polyhedron;
			this->cylinder = n.cylinder;
			this->con = n.con;
			this->phi = n.phi;
			this->isReverseRegion = n.isReverseRegion;
			this->isNormalized = n.isNormalized;
			return *this;
		}
	};
	class PointSet {
	public:
		PointSet(int _generate_step = 0, int _phase_index = 0, REAL _temperature = 0) {
			init(_generate_step, _phase_index, _temperature);
			is_normalized = true;
		}
		void init(int _generate_step = 0, int _phase_index = 0, REAL _temperature = 0) {
			generate_step = _generate_step;
			phaseIndex = _phase_index;
			temperature = _temperature;
		}
		~PointSet() {
			points.clear();
		}
		vector<Point> points;
		vector<REAL> points_phi;
		int generate_step;
		int phaseIndex;
		vector<REAL> con;
		REAL temperature;
		bool is_normalized;
		bool is_point_in_set(Point point, int x_limit, BoundaryCondition x_down_bc, BoundaryCondition x_up_bc, int y_limit = 1, 
			BoundaryCondition y_down_bc = BoundaryCondition::ADIABATIC, BoundaryCondition y_up_bc = BoundaryCondition::ADIABATIC,
			int z_limit = 1, BoundaryCondition z_down_bc = BoundaryCondition::ADIABATIC, BoundaryCondition z_up_bc = BoundaryCondition::ADIABATIC) {
			int x = REAL_to_int(point.x), y = REAL_to_int(point.y), z = REAL_to_int(point.z);
			if (x < 0) {
				if (x_down_bc == BoundaryCondition::PERIODIC)
					x += x_limit;
				else
					x = 0;
			}
			else if (x >= x_limit) {
				if (x_up_bc == BoundaryCondition::PERIODIC)
					x -= x_limit;
				else
					x = x_limit - 1;
			}
			if (y < 0) {
				if (y_down_bc == BoundaryCondition::PERIODIC)
					y += y_limit;
				else
					y = 0;
			}
			else if (y >= y_limit) {
				if (y_up_bc == BoundaryCondition::PERIODIC)
					y -= y_limit;
				else
					y = y_limit - 1;
			}
			if (z < 0) {
				if (z_down_bc == BoundaryCondition::PERIODIC)
					z += z_limit;
				else
					z = 0;
			}
			else if (z >= z_limit) {
				if (z_up_bc == BoundaryCondition::PERIODIC)
					z -= z_limit;
				else
					z = z_limit - 1;
			}

			for (auto p = points.begin(); p < points.end(); p++)
				if (p->x == x && p->y == y && p->z == z)
					return true;
			return false;
		}
		void add_point(REAL _x, REAL _y, REAL _z, REAL _phi) {
			points.push_back(Point(_x, _y, _z));
			points_phi.push_back(_phi);
		}
		void add_point(int _x, int _y, int _z, REAL _phi) {
			points.push_back(Point(_x, _y, _z));
			points_phi.push_back(_phi);
		}
		void add_point(Point _point, REAL _phi) {
			points.push_back(_point);
			points_phi.push_back(_phi);
		}
		PointSet& operator=(const PointSet& n) {
			this->generate_step = n.generate_step;
			this->points = n.points;
			this->temperature = n.temperature;
			this->phaseIndex = n.phaseIndex;
			this->con = n.con;
			this->points_phi = n.points_phi;
			this->is_normalized = n.is_normalized;
			return *this;
		}
	};
	struct NucleationBox {
	public:
		std::vector<pf::GeometricRegion> geometry_box;
		std::vector<pf::PointSet> point_set_box;
		NucleationBox() {
			
		};
		NucleationBox& operator=(const NucleationBox& n) {
			geometry_box = n.geometry_box;
			point_set_box = n.point_set_box;
			return *this;
		}
	};
	struct point_in_region_index {
		point_in_region_index(int re, int gi) {
			region = re;
			grain_index = gi;
		}
		int region;
		int grain_index;
	};
	inline void check_phi_index(int phi_index) {
		if (phi_index < 0 || phi_index >= phi_parameters::phi_number) {
			string _report = "> ERROR : Phi index = " + to_string(phi_index) + " is invalid !\n";
			WriteDebugFile(_report);
			WriteLog(_report);
			SYS_PROGRAM_STOP;
		}
	}
	inline void check_con_size(int con_size) {
		if (con_size != con_parameters::con_number) {
			string _report = "> ERROR : Con size = " + to_string(con_size)
				+ " is not equal to " + to_string(con_parameters::con_number) + " !\n";
			WriteDebugFile(_report);
			WriteLog(_report);
			SYS_PROGRAM_STOP;
		}
	}

	namespace microstructure_init {
		// - Nucleation
		inline NucleationBox nucleation_box;
		inline bool is_datafile_init = false;
		inline bool is_read_datafile_by_path = false;
		inline string datafile_path = "DATA.dat";
		inline pf::write_mesh_data::Data_MeshInfo datafile_report;
		// - matrix
		inline int matrix_phi_index = 0;
		inline REAL matrix_phi_value = 0;
		inline vector<REAL> matrix_con;
		inline REAL matrix_temperature = 0.0;
		// - bmp24
		inline bool is_read_bmp24file = false;
		inline string bmp24file_path = "FIG.bmp";
		inline int bmp24_layer = 0;
		inline vector<vector<REAL>> bmp24_threshold;
		inline vector<int> bmp24_phi_index;
		inline vector<REAL> bmp24_phi_value;
		inline vector<bool> bmp24_phi_normalized;
		inline vector<vector<REAL>> bmp24_con;
		inline vector<REAL> bmp24_temperature;
		// - voronoi
		enum VoronoiType { VDT_CONST, VDT_REF_DOT, VDT_REF_SURFACE, VDT_DOTS_MATRIX };
		inline bool is_voronoi = false;
		inline bool is_voronoi_mirror_generation = true;
		inline VoronoiType voronoi_type = VoronoiType::VDT_CONST;
		inline vector<Point> voronoi_points;
		inline REAL voronoi_const_pointsDistance = -1.0;
		inline Point voronoi_reference_dot;
		inline REAL voronoi_reference_dot_distance = -1.0;
		inline REAL voronoi_reference_dot_min_pointsDistance = -1.0;
		inline REAL voronoi_reference_dot_max_pointsDistance = -1.0;
		inline surf_func_3D voronoi_reference_surface;
		inline REAL voronoi_reference_surface_distance = -1.0;
		inline REAL voronoi_reference_surface_min_pointsDistance = -1.0;
		inline REAL voronoi_reference_surface_max_pointsDistance = -1.0;
		inline vector<Point> voronoi_matrix_dots;
		inline vector<REAL> voronoi_matrix_dots_pointsDistance;
		inline bool is_voronoi_rand = true;
		inline int voronoi_rand_seed = 0;
		inline Vector3 voronoi_box_position = Vector3(0, 0, 0);
		inline Vector3 voronoi_box_size = Vector3(0, 0, 0);
		inline vector<int> voronoi_phi_index_range = { 0, 0 };
		inline vector<REAL> voronoi_con;
		inline REAL voronoi_temperature = 0.0;
		inline vector<int> voronoi_phis_indexs;
		// - porous
		inline bool is_porous = false;
		inline bool is_porous_rand = true;
		inline int porous_rand_seed = 0;
		inline REAL porosity = 2;
		inline REAL porous_init_noise = 0;
		inline int porous_first_phi_index = 0;
		inline int porous_second_phi_index = 0;
		inline vector<REAL> porous_first_con;
		inline vector<REAL> porous_second_con;
		inline REAL porous_first_temperature = 0;
		inline REAL porous_second_temperature = 0;
		inline vector<int> porous_phis_indexs;
		inline bool is_porous_normalized = true;
		inline REAL porous_TwoD_d1 = REAL(0.05), porous_TwoD_d5 = REAL(0.0125);
		inline REAL porous_ThreeD_d1 = REAL(0.02), porous_ThreeD_d7 = porous_ThreeD_d1 / 2, porous_ThreeD_d19 = porous_ThreeD_d1 / 8;
	}
}