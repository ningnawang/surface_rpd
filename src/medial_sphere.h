#ifndef H_medial_sphere_H
#define H_medial_sphere_H

#include <geogram/mesh/mesh.h>

#include <vector>

#include "common_cxx.h"
#include "input_types.h"

// for sphere shrinking
struct ss_params {
  Vector3 p, p_normal, q, q_normal;
  // matching GEO::Mesh facets
  int q_fid = -1;
  int p_fid = -1;
};

class TangentPlane {
 public:
  // TangentPlane();
  // TangentPlane(const Vector3& _normal);
  TangentPlane(const Vector3& _normal, const Vector3& _point, const int _fid);
  ~TangentPlane(){};

 public:
  void print_info() const;
  bool is_same_normal(const Vector3& bnormal,
                      const double eps_degree = EPS_DEGREE_10) const;
  bool is_same_tan_pl(const TangentPlane& tan_pl2) const;
  static double get_energy_value(const Vector3& theta, const double& radius,
                                 const double alpha1, const double alpha2,
                                 const Vector3& tan_point,
                                 const Vector3& normal);
  double update_energy_value(const Vector3& theta, const double& radius,
                             const double alpha1, const double alpha2);

 public:
  Vector3 normal;               // normal of tangent plane
  std::vector<Vector3> points;  // points on plane
  int fid;                      // corresponding fid from sf_mesh, -1 as default
  double energy;                // sphere to minimize this energy function,
                                // DBL_MAX default
  double energy_over_sq_radius;  // energy / sq_radius, used as break threshold

  bool is_deleted;
};

class TangentConcaveLine {
 public:
  // TangentConcaveLine(const int _id, const int _id_fe, const int
  // _num_ce_group)
  //     : id(_id), id_fe(_id_fe), num_ce_group(_num_ce_group){};
  TangentConcaveLine(const int _id, const FeatureEdge& fe);
  ~TangentConcaveLine(){};

  // TangentConcaveLine(const int _id, const aint3 _t2vs_group,
  //                    const avec2 _t2vs_pos, const aint2& _sf_fs_pair,
  //                    const avec2& _adj_normals);

 public:
  void print_info() const;
  // if concave line is a curve
  // then each line should cover more adjacent faces
  bool is_normal_covered_by_adj_fs(
      const Vector3& n, double esp_degree_given = EPS_DEGREE_10) const;
  // given two points of a segment
  void update_tangent_point(const GEO::vec3& v1_p, const GEO::vec3& v2_p);

  static double get_energy_value(const Vector3& theta, const double& radius,
                                 const double alpha3, const Vector3& tan_point,
                                 const Vector3& normal);
  double update_energy_value(const Vector3& theta, const double& radius,
                             const double alpha3);

 public:
  bool operator==(TangentConcaveLine const& b) const;

 public:
  int id;
  int id_fe;                  // matching to TetMesh::feature_edges::id
  int num_ce_group;           // matching TetMesh::feature_edges::t2vs_group[3]
  bool is_tan_point_updated;  // point M will be updated to tangent point X
                              // after each RPD

  Vector3 direction;  // (ref_vs_pos[1] - ref_vs_pos[0]).normalize()
  Vector3 tan_point;  // tangent point, init as center of ref_vs_pos
  Vector3 normal;     // a random vector inside the range of 2 adjacent normals

  double energy;  // distance from sphere to concave line, updated by alpha_3
  double energy_over_sq_radius;  // used as break threshold

  // aint3 t2vs_group;      // edges in format <tvid_min, tvid_max,
  // num_se_group> avec2 t2vs_pos;        // 2 tet_vs indices aint2
  // adj_sf_fs_pair;  // facet id of input GEO::Mesh for RPD avec2 adj_normals;
  // // two normals of adjacent reference plane
};

enum Topo_Status {
  low_edge_euler = -4,
  low_facet_euler = -3,
  low_cell_euler = -2,  // ninwang: site_knn __K__ might be small
  unkown = -1,          // newly added sphere to fix
  ok = 0,
  high_cell_cc = 1,  // ninwang: tet_knn tet_k might be small
  high_facet_cc = 2,
  high_edge_cc = 3
};

// each PowerCell contains multiple ConvexCellTransfer
struct PowerCell {
  PowerCell(){};

  std::set<int> cell_ids;  // matching ConvexCellTransfer::id

  /** For PowerCell Cells **/
  // cell_id -> {neighboring cell ids}
  std::map<int, std::set<int>> cell_neighbors;
  // original tet fid -> 1 or 2 cell ids
  std::map<int, std::set<int>> f_id1_to_cells;
  // cell id -> orignal tet fids
  // first #sf_mesh.facets.nb()-1 matches GEO::Mesh, later are unique fids
  // from orignal tet
  std::map<int, std::set<int>> cell_to_f_id1s;
  // pc_face (not surface triangle) centroids on surface mesh in
  // <pc_face_centroid, surf_fid> pair (<Vector3, int>)
  // (fid defines the same as cell_to_f_id1s, but we only care about the sf_mesh
  // fid that matches GEO::Mesh here).
  // cell id -> vector of <pc_face_centroid, surf_fid>
  std::map<int, std::vector<v2int>> cell_to_surfv2fid;

  std::vector<std::set<int>> cc_cells;  // grouped for each CC
  std::vector<std::set<int>> cc_bfids;  // boundary fids, matching GEO::Mesh or
                                        // just unique id from original tet mesh
  std::vector<std::vector<v2int>>
      cc_surf_v2fids;  // surface vertex2fids, matching GEO::Mesh only

  /** For PowerCell Facets **/
  // halfplane seed_neigh_id -> list of cell ids
  // (all cells that clipped by the halfplane)
  std::map<int, std::set<int>> f_id2_to_cells;
  // store seed_neigh_id that needs to add new spheres around
  // neigh_id -> fixed (true) or not fixed (false)
  std::map<int, bool> f_id2_is_fixed;
  // neigh_id -> { set of cell_ids in one facet CC }
  std::map<int, std::vector<std::set<int>>> facet_cc_cells;
  // neigh_id -> { set of v2fids in one facet CC }
  // relates to cell_to_surfv2fid
  // here we only care about the sf_mesh fid that matches GEO::Mesh
  std::map<int, std::vector<std::vector<v2int>>> facet_cc_surf_v2fids;

  /** For PowerCell Edges **/
  // Each edge in powercell cur_id is uniquely defined by 3 halfplanes [cur_id,
  // neigh_id_min, neigh_id_max] and we use [neigh_id_min, neigh_id_max] to
  // represent this unique edge
  //
  // [neigh_id_min, neigh_id_max] -> list of cell ids
  // (all cells that clipped by 2 halfplanes)
  std::map<aint2, std::set<int>> e_to_cells;
  // store shared edge that needs to add new spheres around
  // [neigh_id_min, neigh_id_max] -> fixed (true) or not fixed (false)
  std::map<aint2, bool> e_is_fixed;
  // [neigh_id_min, neigh_id_max] -> { set of cell_ids in one edge CC }
  std::map<aint2, std::vector<std::set<int>>> edge_cc_cells;
  // [neigh_id_min, neigh_id_max] -> { set of edge endpoints <pos, sf_fid>}
  // each powercell edge is dual to a medial face of 3 spheres
  //
  // (not anymore) used for geometry check, not topology check, and only check
  // when there is only 2 endpoints for each [neigh_id_min, neigh_id_max]
  //
  // used for thinning, store dual medial face's importance
  // sf_fid = -1 if not on surface
  // std::map<aint2, std::vector<v2int>> edge_endpoints;
  //
  // [neigh_id_min, neigh_id_max] -> [aint2, aint2]
  // aint2 is [cell_id, lvid (from cc_trans.nb_v)]
  std::map<aint2, std::array<aint2, 2>> edge_2endvertices;

  /** For PowerCell Vertices **/
  // Each vertex in powercell cur_id is uniquely defined by
  // [cur_id, neigh_id1, neigh_id2, neigh_id3]
  // and we use sorted [neigh_id1, neigh_id2, neigh_id3] to represent this
  // unique vertex
  //
  // 1. If neigh_id1 is -1, then this powercell vertex is on surface
  // 2. If all neigh_id(i) > 0, then vid is inside the shape. If this vertex is
  // a powercell vertex (not just a convex cell vertex), then it is dual to a
  // medial tet
  //
  // [cell_id, lvid (from cc_trans.nb_v)] -> [neigh_id1, neigh_id2, neigh_id3]
  std::map<aint2, aint3> vertex_2id;
  // [cell_id, lvid (from cc_trans.nb_v)] -> <pos, sf_fid>
  // sf_fid can be -1 if vertex is not on surface
  std::map<aint2, v2int> vertex_2pos;

  /** For External Edge Features **/
  // all touched sharp edge, stored in [cell_id, lvid1, lvid2,
  // num_se_group]
  std::set<aint4> se_covered_lvids;

  /** For Internal Features **/
  // grouped sf_mesh <pc_face_centroid, surf_fid> pairs (<Vector3, int>)
  // not crossing sharp edges
  std::vector<std::vector<v2int>> surf_v2fid_in_groups;

  Topo_Status topo_status;
};

enum SphereType {
  T_UNK = -1,
  T_1_2 = -2,  // external feature, sharp edge
  T_1_N = -3,  // external feature corners
  T_2 = 2,
  T_N = 3,     // including T_3, T_4 ...
  T_c = 10,    // invisible pin sphere added on concave line to avoid RPD
               // degeneration
  T_2_c = 11,  // added for one concave line using sphere shrinking
  T_N_c = 12   // added through internal feature preservation
};

class MedialSphere {
 public:
  MedialSphere(int _id, Vector3 _pin, Vector3 _pin_normal,
               SphereType _type = SphereType::T_2);
  MedialSphere(int _id, Vector3 _center, double _radius, SphereType _type);
  ~MedialSphere(){};

 public:
  void print_info() const;
  void print_tan_planes() const;
  double update_all_energy_values(const double alpha1, const double alpha2);
  void save_old_center_radius(bool is_clear = true);
  void update_tan_planes_from_ss_params();
  void new_tan_plane_no_dup(const Vector3& _normal, const Vector3& _point,
                            const int _fid);
  void new_cc_line_no_dup(const FeatureEdge& one_fe);
  void apply_perturb();

  bool is_on_se() const;
  bool is_on_corner() const;
  bool is_on_ce_pin() const;
  bool is_on_ce() const;
  bool is_on_intf() const;
  bool is_on_extf() const;

  void topo_clear();
  void pcell_insert(int cell_id);
  void fcc_fixed(int neigh_id);

  bool operator==(const MedialSphere& m2) const;
  bool operator!=(const MedialSphere& m2) const;

 public:
  int id;
  bool is_deleted = false;
  Vector3 center;
  double radius;
  SphereType type;
  ss_params ss;
  std::vector<TangentPlane> tan_planes;
  std::vector<TangentConcaveLine> tan_cc_lines;

  /* For topology check */
  PowerCell pcell;  // one powercell per sphere
  // Euler characteristics for all cells
  double euler = 0;  // = euler_sum - num_cells
  double euler_sum = 0;
  uint num_cells = 0;

  /* For medial mesh */
  std::set<int> edges_;  // matches MedialMesh::edges
  std::set<int> faces_;  // matches MedialMesh::faces
  void clear_mm();

  // for relaxation (CVT-ish) and sphere iterating
  Vector3 old_center;
  double old_radius;

  /* For external feature */
  int num_se_group = -1;
};

bool validate_new_sphere(const std::vector<MedialSphere>& all_medial_spheres,
                         const MedialSphere& new_sphere);
bool add_new_sphere_validate(std::vector<MedialSphere>& all_medial_spheres,
                             MedialSphere& new_sphere);

void purge_deleted_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres);

bool is_two_mspheres_on_same_se(const MedialSphere& msphere1,
                                const MedialSphere& msphere2);

#endif  // __H_medial_sphere_H__