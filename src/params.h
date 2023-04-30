#pragma once

// user-defined
struct Parameter {
  float bbox_diag_l = 0.f;
  std::vector<float> bb_points;  // 8 bbox

  float scale_max = 1000.f;  // normalize models between [0,1000]^3

  // For finding sphere neighbors (using RT)
  // will be updated by function get_RT_vertex_neighbors()
  int site_k = 90;

  // For adding spheres on feature edges if too long
  float fe_rel_len = 1. / 50.;

  // For measuring the Hausdorff distance
  // between medial slab and the surface
  float hd_rel_slab = 1. / 60.;

  // For adding medial spheres that tangent to concave lines using
  // ball-shrinking algo. We need to specify two parameters:
  // 1. cc_len_eps_rel:
  //    defines the length between to pin points on concave lines, relative to
  //    bbox_diag_l
  // 2. cc_normal_eps:
  //    define the angle between two normals given a pin point, scaled in
  //    [0,360]
  // Given a pin point and a normal, we can add a sphere using ball-shrinking
  // algo. Concave lines can be better preserved when cc_len_eps is smaller
  // (more medial spheres inserted), but it also means more processing time.
  double cc_len_eps_rel = 1. / 5.;
  double cc_normal_eps = 90;

  // Threshold for detecting sharp/concave edges and corners
  double thres_concave = 0.3;  // smaller more sensative
  double thres_convex = 30.;   // smaller more sensative
};

#define KNEARESTS_MIN_N 14000

// random seed
#define RAN_SEED 200
#define RANDOM_01() ((double)std::rand() / RAND_MAX)  // random double in [0, 1]
#define RANDOM_INT(l, h) \
  (l + std::rand() % (h - l))  // random int in [l, h)
                               // (https://stackoverflow.com/a/7560146)

#define FOR(I, UPPERBND) for (int I = 0; I < int(UPPERBND); ++I)
#define SCALAR_ZERO 1e-3
#define SCALAR_ZERO_6 1e-6  // float, 1e-8 if double
#define SCALAR_ZERO_5 1e-5
#define SCALAR_ZERO_4 1e-4
#define SCALAR_ZERO_3 1e-3
#define SCALAR_ZERO_2 1e-2
#define SCALAR_ZERO_1 1e-1
#define SCALAR_FEATURE_RADIUS 1e-2
#define SCALAR_SE_MERGE_RADIUS 5 /* scaled to [0,1000]^3 */
#define SCALAR_CE_PIN_RADIUS 30

#define PRESET 3

#if 0 == PRESET  // conservative settings (white noise)
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 32
#define _K_ 90  // k_site nearest sites of each site
// #define _TET_K_         70
#define _MAX_P_ 64
#define _MAX_T_ 96
#define _MAX_E_ 152
#elif 1 == PRESET  // perturbed grid settings
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 32
#define _K_ 90
#define _MAX_P_ 50
#define _MAX_T_ 96
#elif 2 == PRESET  // blue noise settings
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 64
#define _K_ 35
#define _MAX_P_ 32
#define _MAX_T_ 96
#elif 3 == PRESET  // ninwang setting
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 32
#define _K_ 150  // for testing
#define _MAX_P_ 64
#define _MAX_T_ 96
#define _MAX_E_ 152
typedef double Scalar;  // use for calculating Euler
#endif

// Uncomment to activate arithmetic filters.
//   If arithmetic filters are activated,
//   status is set to needs_exact_predicates
//   whenever predicates could not be evaluated
//   using floating points on the GPU
// #define USE_ARITHMETIC_FILTER

#define IF_VERBOSE(x) x

// #define USE_NINWANG
