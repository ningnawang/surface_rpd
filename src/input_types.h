#pragma once
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>

#include <memory>

#include "common_cxx.h"
#include "params.h"

enum EdgeType {
  UE = -1,  // unkown edge
  SE = 1,   // convex sharp edge
  CE = 2    // concave edge
};

// either sharp edge or concave edge
class FeatureEdge {
 public:
  FeatureEdge(const int _id, const EdgeType &_t, const aint3 &_t2vs_group)
      : id(_id), type(_t), t2vs_group(_t2vs_group){};
  ~FeatureEdge(){};
  void print_info() const;

  int id;
  EdgeType type;
  // (store in tet mesh indices, matching TetMesh::tet_vertices)
  aint3 t2vs_group;      // edges in format <tvid_min, tvid_max, num_se_group>
  avec2 t2vs_pos;        // tet_vs indices of 2 endpoints
  aint2 adj_sf_fs_pair;  // mapping to GEO::Mesh sf_mesh
  avec2 adj_tan_points;  // centroid of adj_sf_fs_pair
  avec2 adj_normals;     // order maps adj_sf_fs_pair
};

class TetMesh {
 public:
  TetMesh(std::string path) : tet_path_with_ext(path){};

 public:
  std::string tet_path_with_ext;
  std::vector<float> tet_vertices;  // tet vertices
  std::vector<int> tet_indices;     // tet 4 indices of vertices
  std::map<int, std::set<int>> v2tets;

  // mapping to SurfaceMesh
  std::map<int, std::set<int>> tet_vs2sf_fids;

  // for both sharp edges and concave edges
  std::vector<FeatureEdge> feature_edges;

  // for sharp edge
  // format: <tid, lfid_min, lfid_max> -> <num_se_group>.
  // the mapping can be many-to-one (multiple tets share same sharp edge).
  std::map<aint3, int> se_lfs2group_map;
  // mapping vertex from <tid, lfid1, lfid2, lfid3> -> tvid
  std::map<aint4, int> tet_vs_lfs2tvs_map;

  // for concave edges, vector index matching
  // FeatureEdge::t2vs_group[3] (num_ce_group)
  std::vector<std::vector<aint3>> ce_tet_groups;

  // for corner
  std::set<int> corners_tet;
};

class AABBWrapper {
 private:
  std::shared_ptr<GEO::MeshFacetsAABB> sf_tree;

 public:
  AABBWrapper() {}

  // is_reorder == true => will reorder sf_mesh
  void init_sf_mesh_and_tree(GEO::Mesh &_sf_mesh, bool is_reorder = true) {
    sf_tree = std::make_shared<GEO::MeshFacetsAABB>(_sf_mesh, is_reorder);
  }

 public:
  inline int project_to_sf_get_nearest_face(Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet(p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];
    return fidx;
  }

  inline double project_to_sf(Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    sf_tree->nearest_facet(p, nearest_p, sq_dist);
    p = nearest_p;
    return sq_dist;
  }

  inline int get_nearest_face_sf(const Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    return sf_tree->nearest_facet(p, nearest_p, sq_dist);
  }

  inline int get_nearest_face_sf(const Vector3 &p, double &sq_dist) const {
    Vector3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    return sf_tree->nearest_facet(p, nearest_p, sq_dist);
  }

  inline int get_nearest_point_on_sf(const Vector3 &p, Vector3 &nearest_p,
                                     double &sq_dist) const {
    sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet(p, nearest_p, sq_dist);
    return fidx;
  }

  inline double get_sq_dist_to_sf(const Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    sf_tree->nearest_facet(p, nearest_p, sq_dist);
    return sq_dist;
  }
};

class SurfaceMesh : public GEO::Mesh {
 public:
  void reload_sf2tet_vs_mapping();

 public:
  AABBWrapper aabb_wrapper;
  std::vector<int> sf2tet_vs_mapping;  // matching TetMesh::tet_vertices

  // for feature edges
  // store <sf_fid_min, sf_fid_max> if on feature edge
  std::set<aint2> fe_sf_fs_pairs;
};

void load_sf_tet_mapping(const GEO::Mesh &sf_mesh,
                         std::map<int, int> &vs_tet2sf_mapping,
                         std::map<int, std::set<int>> &vs2fids,
                         std::map<int, std::set<int>> &tet_vs2sf_fids);

void get_k_ring_neighbors_no_cross(const GEO::Mesh &sf_mesh,
                                   const std::set<aint2> &fe_sf_pairs_not_cross,
                                   const int fid_given, const int k,
                                   std::set<int> &k_ring_fids, bool is_debug);

void store_special_edge(const TetMesh &tet_mesh, const SurfaceMesh &sf_mesh,
                        const EdgeType &fe_type, const aint3 &t2vs_group,
                        std::vector<FeatureEdge> &feature_edges);