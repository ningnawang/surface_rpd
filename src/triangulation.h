#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <geogram/delaunay/periodic_delaunay_3d.h>

#include "medial_sphere.h"
#include "params.h"

//////////
// CGAL

// Regular Vertex info (dual to sphere)
class RVI {
 public:
  int tag = -1;  // matching valid_id
  // matching index of all_medial_spheres
  // remains -1 if 8 bbox points
  int all_id = -1;
};

// Regular Tetrahedron info
class RTI {
 public:
  int id = -1;  // assigned while iterating RT, not unique for each run
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kf;
typedef Kf::FT Weight;
typedef Kf::Point_3 Point;
typedef Kf::Weighted_point_3 Weighted_point;

typedef CGAL::Regular_triangulation_vertex_base_3<Kf> Vb0_rt;
typedef CGAL::Triangulation_vertex_base_with_info_3<RVI, Kf, Vb0_rt> Vb_rt;
typedef CGAL::Regular_triangulation_cell_base_3<Kf> Cb0_rt;
typedef CGAL::Triangulation_cell_base_with_info_3<RTI, Kf, Cb0_rt> Cb_rt;
typedef CGAL::Triangulation_data_structure_3<Vb_rt, Cb_rt> Tds_rt;
typedef CGAL::Regular_triangulation_3<Kf, Tds_rt> Rt;
typedef CGAL::Triangulation_3<Kf, Tds_rt> Tr;

typedef Kf::Point_3 Point_rt;
typedef Rt::Vertex_iterator Vertex_iterator_rt;
typedef Rt::Vertex_handle Vertex_handle_rt;
typedef Rt::Cell_iterator Cell_iterator_rt;
typedef Rt::Cell_handle Cell_handle_rt;
typedef Rt::Cell_circulator Cell_circulator_rt;
typedef Rt::Facet_circulator Facet_circulator_rt;

typedef Rt::Finite_cells_iterator Finite_cells_iterator_rt;
typedef Rt::Finite_facets_iterator Finite_facets_iterator_rt;
typedef Rt::Finite_edges_iterator Finite_edges_iterator_rt;
typedef Rt::Finite_vertices_iterator Finite_vertices_iterator_rt;
typedef Rt::Tetrahedron Tetrahedron_rt;

///////////////
class RegularTriangulationNN : public Rt, public GEO::Counted {
 public:
  // RegularTriangulationNN();
  /**
   * \brief RegularTriangulationNN destructor
   */
  ~RegularTriangulationNN(){};

 public:
  ///////////////
  inline void clean() {
    this->clear();
    tag_to_vh.clear();
    nb_vertices = 0;
  }

  // inline void print_info() {
  //   logger().info("------ Regular Triangulation Info ------");
  //   logger().info("nb_vertices: {}", nb_vertices);
  //   logger().info("#v: {}", number_of_vertices());
  //   logger().info("#e: {}", number_of_finite_edges());
  //   logger().info("#f: {}", number_of_finite_facets());
  // }
  // inline void print_dual_info() {
  //   logger().info("------ Regular Triangulation Dual Info ------");
  //   logger().info("#tets dual vs: {}", rt_ts_info.size());
  //   logger().info("#f dual segs: {}", rt_fs_info.size());
  //   logger().info("#e dual polys: {}", rt_es_info.size());
  // }

  // ATTENTION:
  // we store valid_medial_spheres.size() here
  // is because rt tag range is [0, valid_medial_spheres.size()]
  // and RPD would use this tag for flagging
  inline void set_nb_vertices(const int nb_v) { nb_vertices = nb_v; }

  inline int get_nb_vertices() const { return nb_vertices; }

  inline Vector3 to_geo_vec(const Point_rt& wp) const {
    return Vector3(wp.x(), wp.y(), wp.z());
  };

  inline double tet_radius(const Tetrahedron_rt& tet) const {
    return GEO::length(to_geo_vec(tet.vertex(0)) -
                       to_geo_vec(CGAL::circumcenter(tet)));
  };

  inline Vertex_handle_rt get_vh(const GEO::index_t tag) const {
    if (tag_to_vh.empty()) {
      printf("tag_to_vh cannot be empty\n");
      assert(false);
    } else if (tag_to_vh.find(tag) == tag_to_vh.end()) {
      printf("tag %d cannot be found at tag_to_vh \n", tag);
      printf("tag not match tag_to_vh value\n");
      assert(false);
    } else if (tag != tag_to_vh.at(tag)->info().tag) {
      printf("tag %d has tag_to_vh %d \n", tag, tag_to_vh.at(tag)->info().tag);
      printf("tag not match tag_to_vh value\n");
      assert(false);
    } else
      return tag_to_vh.at(tag);
  };

  inline void set_tag_to_vh(int tag, Vertex_handle_rt& vh) {
    tag_to_vh.insert(std::pair<int, Vertex_handle_rt>(tag, vh));
  }

  inline std::vector<double> get_double_vector(const Weighted_point& wp) const {
    std::vector<double> p;
    p.push_back(CGAL::to_double(wp.x()));
    p.push_back(CGAL::to_double(wp.y()));
    p.push_back(CGAL::to_double(wp.z()));
    // this is necessary, since we're calling RPD Vertex::intersect_geom()
    p.push_back(CGAL::to_double(wp.weight()));  // weight,
    return p;
  }

  inline std::vector<double> get_double_vector(const GEO::index_t tag) const {
    const Vertex_handle_rt vh = get_vh(tag);
    return get_double_vector(vh->point());
  };

  /////////////////////////////////////////////////////////
  /////  TODO: try to deprecate following functions
  ///// 		 this type of memory allocation might cause memory leak

  // need to delete[] after calling
  inline double* get_double_data(const Weighted_point& wp) const {
    // new operator creates variable p in the HEAP
    // instead of StACK which willl be destroyed
    // when the function is finished
    // checkout this video
    // https://www.youtube.com/watch?v=RWNM7CzDNyY
    double* p = new double[4];
    p[0] = CGAL::to_double(wp.x());
    p[1] = CGAL::to_double(wp.y());
    p[2] = CGAL::to_double(wp.z());
    // p[3] = 0; // weight
    p[3] = CGAL::to_double(wp.weight());  // weight
    return p;
  };

  // need to delete[] after calling
  inline double* get_double_data(const Vector3& wp) const {
    // new operator creates variable p in the HEAP
    // instead of StACK which willl be destroyed
    // when the function is finished
    // checkout this video
    // https://www.youtube.com/watch?v=RWNM7CzDNyY
    double* p = new double[3];
    p[0] = wp[0];
    p[1] = wp[1];
    p[2] = wp[2];
    return p;
  }

  // need to delete[] after calling
  inline double* get_double_data(const GEO::index_t tag) const {
    Vertex_handle_rt vh = get_vh(tag);
    return get_double_data(vh->point());
  };
  /////////////////////////////////////////////////////////

  inline double get_weight(const Weighted_point& wp) const {
    return CGAL::to_double(wp.weight());
  }
  inline double get_weight(const GEO::index_t tag) const {
    const Vertex_handle_rt& vh = get_vh(tag);
    return CGAL::to_double(vh->point().weight());
  }

 protected:
  // vertex tag -> vertex handle
  std::map<int, Vertex_handle_rt> tag_to_vh;

  // nb_vertices >= number_of_vertices()
  int nb_vertices;
};

/**
 * \brief Smart pointer that refers to a Regular Triangulation object
 */
typedef GEO::SmartPointer<RegularTriangulationNN> RegularTriangulationNN_var;

void generate_RT_CGAL_and_purge_spheres(
    const Parameter& params, std::vector<MedialSphere>& all_medial_spheres,
    RegularTriangulationNN& rt);

int get_RT_vertex_neighbors(const RegularTriangulationNN& rt, const int& n_site,
                            std::vector<int>& site_knn);

void get_PC_vertices(const RegularTriangulationNN& rt,
                     std::vector<Vector3>& pc_vertices);