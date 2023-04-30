#pragma once
#include <geogram/basic/geometry.h>
#include <geogram/basic/matrix.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <array>
#include <ctime>
#include <fstream>
#include <queue>
#include <string>
#include <unordered_set>

#include "assert.h"
#include "params.h"

typedef GEO::vec2 Vector2;
typedef GEO::vec3 Vector3;  // all coordinates are initialized to 0 (zero).
typedef GEO::vec4 Vector4;
typedef GEO::vec3i Vector3i;
typedef GEO::mat2 Matrix2;
typedef GEO::mat3 Matrix3;  //  initializes the matrix to the identity matrix
typedef GEO::mat4 Matrix4;
typedef unsigned char uchar;  // local indices with special values
typedef unsigned int uint;
typedef std::array<int, 2> aint2;  //  unique id of powercell edge, medial edge
typedef std::array<int, 3> aint3;  //  medial face
typedef std::array<int, 4> aint4;
typedef std::array<double, 3> adouble3;  // for gui
typedef std::pair<Vector3, int>
    v2int;  // mapping from vertex to one sf_mesh fid
typedef std::array<Vector3, 2> avec2;

// === Math utilities
const double PI = 3.14159265358979323;

constexpr int DIM = 3;
constexpr uchar END_OF_LIST = 255;
constexpr int UNK_INT = -1;
constexpr float UNK_FLOAT = -1.f;
constexpr uchar UNK_UCHAR = UCHAR_MAX;
// every tet face is shared by 2 tets
constexpr float F_TET_ADJ_DEFAULT = 2.f;
// every new clip of convex cell is only shared by 1
// bcs we only care about the Euler of a single cell
constexpr float F_CELL_ADJ_DEFAULT = 1.f;
constexpr float UNK_FACE_ID = -1.f;
constexpr float INIT_RADIUS =
    1000.f;  // for sphere shriking, matching Parameter::scale_max
constexpr double HALF_PI = M_PI / 180.;
constexpr double EPS_DEGREE_10 = 10;
constexpr double EPS_DEGREE_30 = 30;
constexpr double EPS_DEGREE_90 = 90;

// Feature Edges
constexpr int SHARP_EDGE = 1;
constexpr int CONCAVE_EDGE = 2;

// 4 faces of tet abcd: cbd acd bad abc (in vids)
constexpr int tet_faces_lvid_host[4][3] = {
    {2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
// 6 edges of tet (in vids)
constexpr int tet_edges_lvid_host[6][2] = {{2, 3}, {1, 3}, {1, 2},
                                           {0, 3}, {0, 2}, {0, 1}};
// 4 vertices of tet (in lfids)
constexpr int tet_vs_lfid_host[4][3] = {
    {1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};
// 6 edges of tet (in lfids)
constexpr int tet_edges_lfid_host[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                           {1, 2}, {1, 3}, {2, 3}};

template <typename T>
struct is_array_or_vector {
  enum { value = false };
};

template <typename T, typename A>
struct is_array_or_vector<std::vector<T, A>> {
  enum { value = true };
};

template <typename T, std::size_t N>
struct is_array_or_vector<std::array<T, N>> {
  enum { value = true };
};

/**
 * @brief Compute A * transpose(B)
 *
 * @param A 3x3 vector
 * @param B 3x3 vector
 * @return Matrix3
 */
inline Matrix3 vec_vec_trans(Vector3 A, Vector3 B) {
  std::array<double, 9> result = {{A[0] * B[0], A[0] * B[1], A[0] * B[2],
                                   A[1] * B[0], A[1] * B[1], A[1] * B[2],
                                   A[2] * B[0], A[2] * B[1], A[2] * B[2]}};
  return Matrix3(result.data());
}

inline std::string get_timestamp() {
  auto now = std::time(nullptr);
  std::ostringstream os;
  // os << std::put_time(std::gmtime(&now),"%F  %T");
  os << std::put_time(std::localtime(&now), "%F_%T");
  return os.str();
}

inline bool is_file_exist(const std::string filePath) {
  std::ifstream infile(filePath);
  return infile.good();
}

inline std::string get_file_ext(const std::string filePath) {
  return filePath.substr(filePath.find_last_of(".") + 1);
}

inline std::string get_file_no_ext(std::string filePath) {
  std::string filename = filePath.substr(0, filePath.find_last_of("."));
  return filename;
}

inline std::string get_only_file_name(std::string filePath, bool withExtension,
                                      char seperator = '/') {
  std::string filename_ext =
      filePath.substr(filePath.find_last_of(seperator) + 1);
  if (withExtension) return filename_ext;
  size_t lastindex = filename_ext.find_last_of(".");
  return filename_ext.substr(0, lastindex);
}

// using .geogram for saving vertex attributes
//
// (NO) matching ftetwild output (https://github.com/wildmeshing/fTetWild)
//
// type:
// 0 -> .geogram
// 1 -> .obj
// 2 -> .xyz
inline std::string get_other_file_path(std::string filePath, int type) {
  std::string file_path = get_file_no_ext(filePath);
  switch (type) {
    case 0:
      return file_path + "_sf.geogram";
      break;
    case 1:
      return file_path + "_sf.obj";
    case 2:
      return file_path + "_pts.xyz";
    default:
      break;
  }
  assert(false);
}

inline double get_triangle_area(const Vector3& v1, const Vector3& v2,
                                const Vector3& v3) {
  // Heron's formula
  double l1 = GEO::length(v1 - v2);
  double l2 = GEO::length(v1 - v3);
  double l3 = GEO::length(v2 - v3);
  double p = (l1 + l2 + l3) / 2;
  double s = std::sqrt(p * (p - l1) * (p - l2) * (p - l3));
  return s;
}

inline Vector3 get_triangle_centroid(const Vector3& v1, const Vector3& v2,
                                     const Vector3& v3) {
  return (v1 + v2 + v3) / 3.;
};

// Return unit normal given 3 points
inline Vector3 get_normal(const Vector3& a, const Vector3& b,
                          const Vector3& c) {
  // facet abc is clockwise/counter-clockwise
  // (a-c) x (b-c) is clockwise/counter-clockwise
  return GEO::normalize(GEO::cross(a - c, b - c));
}

// unit vector from a to b
inline Vector3 get_direction(const Vector3& a, const Vector3& b) {
  return GEO::normalize(b - a);
}

inline void get_random_vector_between_two_vectors(const Vector3& nA,
                                                  const Vector3& nB,
                                                  Vector3& nX) {
  // double rand_k = ((double)std::rand() / (RAND_MAX));
  std::srand(RAN_SEED);
  Scalar rand_k = RANDOM_01();
  nX = GEO::normalize(rand_k * nA + (1. - rand_k) * nB);
  // logger().debug("rand_k: {}, nA ({},{},{}), nB: ({},{},{}), nX: ({},{},{})",
  // 	rand_k, nA[0], nA[1], nA[2], nB[0], nB[1], nB[2], nX[0], nX[1], nX[2]);
}

inline double angle_between_two_vectors_in_degrees(const Vector3& a,
                                                   const Vector3& b) {
  // a and b should be normalized, so a.dot(b) in range [-1, 1]
  // However, with floating-point math, this is not necessarily true
  // so we need to clamp into [-1, 1]
  Vector3 a_normalized = GEO::normalize(a);
  Vector3 b_normalized = GEO::normalize(b);
  double ab_dot = GEO::dot(a_normalized, b_normalized);
  if (ab_dot <= -1.0) {
    return 180.;
  } else if (ab_dot >= 1.0) {
    return 0.;
  }
  // ab_dot in (-1, 1)
  double diff_angle = std::acos(ab_dot);
  // logger().debug("calculate angle ({},{},{}), ({},{},{}), diff_angle {}",
  //     a[0], a[1], a[2], b[0], b[1], b[2], diff_angle
  // );
  diff_angle *= (180. / PI);
  return diff_angle;
}

// total k, k >= 2 includes 2 boundary normals
inline void sample_k_vectors_given_two_vectors(const Vector3& nA,
                                               const Vector3& nB, const int k,
                                               std::vector<Vector3>& nXs) {
  nXs.clear();
  if (k < 2) return;
  double step = 1. / std::max(k - 1, 1);
  for (double t = 0; t <= 1.; t += step) {
    Vector3 nX = GEO::normalize((1. - t) * nA + t * nB);
    nXs.push_back(nX);
  }
  if (nXs.size() == k - 1) {
    double t = 1;
    Vector3 nX = GEO::normalize((1. - t) * nA + t * nB);
    nXs.push_back(nX);
  }
  if (k != nXs.size()) {
    printf("[ERROR] did not created k: %d normlas with step %f, created %ld \n",
           k, step, nXs.size());
    assert(false);
  }
}

// Project p onto line defined by [v0, v1]
inline void project_point_onto_line(const Vector3& p, const Vector3& v0,
                                    const Vector3& v1, Vector3& p_proj,
                                    double& dist) {
  Vector3 v0p(p - v0), v0v1(v1 - v0);
  p_proj = GEO::dot(v0 + v0p, v0v1) / GEO::dot(v0v1, v0v1) * v0v1;
  dist = (p - p_proj).length();
}

// Return normalized vector
inline Vector3 get_mesh_facet_normal(const GEO::Mesh& mesh, const int fidx) {
  Vector3 fn = GEO::normalize(GEO::Geom::mesh_facet_normal(mesh, fidx));
  return fn;
}

inline Vector3 get_mesh_facet_centroid(const GEO::Mesh& mesh, const int fid) {
  assert(fid >= 0 || fid < mesh.facets.nb());
  return get_triangle_centroid(mesh.vertices.point(mesh.facets.vertex(fid, 0)),
                               mesh.vertices.point(mesh.facets.vertex(fid, 1)),
                               mesh.vertices.point(mesh.facets.vertex(fid, 2)));
}

// Return normalized vector
inline Vector3 get_mesh_vertex_normal(const GEO::Mesh& mesh, const int vidx) {
  Vector3 vn = GEO::normalize(GEO::Geom::mesh_vertex_normal(mesh, vidx));
  return vn;
}

inline Vector3 get_random_point_given_facet(const Vector3& p1,
                                            const Vector3& p2,
                                            const Vector3& p3) {
  std::srand(RAN_SEED);
  // https://stackoverflow.com/a/21722167
  Scalar r1 = RANDOM_01();
  Scalar r2 = RANDOM_01();
  if (r1 + r2 > 1.f) {
    r1 = 1.f - r1;
    r2 = 1.f - r2;
  }
  std::array<Scalar, 3> w = {{1 - r1 - r2, r1, r2}};
  return w[0] * p1 + w[1] * p2 + w[2] * p3;
}

inline Vector3 get_random_point_given_edge(const Vector3& p1,
                                           const Vector3& p2) {
  std::srand(RAN_SEED);
  // https://stackoverflow.com/a/21722167
  Scalar r = RANDOM_01();
  std::array<Scalar, 2> w = {{1 - r, r}};
  return w[0] * p1 + w[1] * p2;
}

inline bool all_finite(const Vector3& p) {
  for (uint i = 0; i < 3; i++) {
    if (std::isnan(p[i]) || !std::isfinite(p[i])) return false;
  }
  return true;
}

// make sure f1_normal and f2_normal are normalized!!
inline bool is_share_a_sharp_edge(const Vector3& f1_normal,
                                  const Vector3& f2_normal) {
  double cos = GEO::dot(f1_normal, f2_normal);
  // sharp edge: theta in [90, 180), cos in (-1, 0]
  if (cos > -1 && cos < SCALAR_ZERO_3) {
    return true;  // skip checking its neighbor
  }
  return false;
}

inline bool is_vector_same_direction(const Vector3& a, const Vector3& b,
                                     const double degree) {
  // angle betwen a and b is in [0, degree]
  // const double halfC = M_PI / 180;HALF_PI
  return GEO::dot(a, b) / (a.length() * b.length()) >=
         std::cos(degree * HALF_PI);
}

inline double get_distance_between_two_vectors(const Vector3& a,
                                               const Vector3& b) {
  return (a - b).length();
}

template <typename T>
inline void vector_unique(std::vector<T>& v) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

template <typename T>
inline T get_sorted(const T& in) {
  T out = in;
  std::sort(out.begin(), out.end());
  return out;
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      std::set<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  for (uint x : s1) {
    if (s2.count(x)) {
      v.insert(x);
    }
  }
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      std::vector<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.push_back(x);
    }
  }
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      const std::set<T>& s3, std::set<T>& v) {
  v.clear();
  std::set<T> neigh_cells_tmp;
  set_intersection<T>(s1, s2, neigh_cells_tmp);
  if (neigh_cells_tmp.empty()) return;
  set_intersection<T>(neigh_cells_tmp, s3, v);
}

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2, std::vector<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.push_back(x);
    }
  }
}

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2,
                      const std::unordered_set<T>& s3, std::vector<T>& v) {
  if (s2.size() < s1.size() && s2.size() < s1.size()) {
    set_intersection(s2, s1, s3, v);
    return;
  }

  if (s3.size() < s1.size() && s3.size() < s2.size()) {
    set_intersection(s3, s1, s2, v);
    return;
  }

  assert(s1.size() <= s2.size());
  assert(s1.size() <= s3.size());

  v.clear();
  v.reserve(s1.size());
  for (int x : s1) {
    if (s2.count(x) && s3.count(x)) {
      v.push_back(x);
      if (v.size() == 2) break;
    }
  }
}

template <typename T>
std::set<T> to_set(const std::vector<T>& in) {
  std::set<T> out;
  for (int x : in) {
    out.insert(x);
  }
  return out;
}

inline int mod3(int j) { return j % 3; }
inline int mod4(int j) { return j % 4; }
inline int modk(int j, int k) { return j % k; }

/**
 * @brief Given cells/others to trace, and return the number of connected
 * components
 *
 * @param cell_to_visit
 * @param cell_neighbors
 * @param cc_cells return, grouped in each CC
 * @return int number of CC
 */
template <typename T>
inline int get_CC_given_neighbors(
    const std::set<T>& cell_to_visit,
    const std::map<T, std::set<T>>& cell_neighbors,
    std::vector<std::set<T>>& cc_cells) {
  std::set<T> cell_unvisited = cell_to_visit;  // copy
  cc_cells.clear();

  // printf("cell_unvisited: %ld, cell_neighbors size: %ld \n",
  //        cell_unvisited.size(), cell_neighbors.size());
  int num_cc = 0;
  std::set<T> cell_visited, one_cc_visited;
  // calculate the number of CCs
  std::queue<T> cc_queue;
  while (!cell_unvisited.empty()) {
    cc_queue.push(*cell_unvisited.begin());
    while (!cc_queue.empty()) {
      T cell_id = cc_queue.front();
      cc_queue.pop();

      if (cell_visited.find(cell_id) != cell_visited.end()) continue;
      cell_visited.insert(cell_id);    // all visited cells
      cell_unvisited.erase(cell_id);   // all unvisited cells
      one_cc_visited.insert(cell_id);  // visited cells in one CC
      // push all neighbors
      if (cell_neighbors.find(cell_id) == cell_neighbors.end()) {
        // assert(false);
        continue;
      }
      auto& neighbors = cell_neighbors.at(cell_id);
      for (auto& neigh_id : neighbors) {
        // if not in cell_to_visit, do not add to queue
        if (cell_to_visit.find(neigh_id) != cell_to_visit.end())
          cc_queue.push(neigh_id);
      }
    }
    // printf("found one CC with cell_visited: %zu/%zu, cell_unvisited:
    // %zu/%zu\n",
    //        cell_visited.size(), cell_to_visit.size(), cell_unvisited.size(),
    //        cell_to_visit.size());
    num_cc++;
    cc_cells.push_back(one_cc_visited);  // store cell ids in one CC
    one_cc_visited.clear();
  }

  return num_cc;
}