#include "triangulation.h"

/**
 * @brief Given sphere + 8 bbox, we generate RT. Since some spheres may not
 * exist in RT (bcs of weights), we purge the all_medial_spheres to those valid
 * in RT.
 *
 * @param params provides 8 bbox info
 * @param all_medial_spheres will be updated after RT
 * @param rt
 */
void generate_RT_CGAL_and_purge_spheres(
    const Parameter& params, std::vector<MedialSphere>& all_medial_spheres,
    RegularTriangulationNN& rt) {
  int num_spheres = all_medial_spheres.size();
  printf("[RT] generate RT for %d spheres\n", num_spheres);
  rt.clean();
  // add all medial spheres
  for (int mid = 0; mid < num_spheres; mid++) {
    const MedialSphere& msphere = all_medial_spheres.at(mid);
    // do not add deleted sphere in RT
    // so later we can purge them after RT
    if (msphere.is_deleted) continue;
    Point_rt p(msphere.center[0], msphere.center[1], msphere.center[2]);
    Weight weight = std::pow(msphere.radius, 2);
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
    // update RT
    vh->info().all_id = msphere.id;
  }

  // add 8 bbox, vh->info().all_id = -1
  assert(params.bb_points.size() / 3 == 8);
  for (int i = 0; i < 8; i++) {
    Point_rt p(params.bb_points[i * 3], params.bb_points[i * 3 + 1],
               params.bb_points[i * 3 + 2]);
    Weight weight = SCALAR_FEATURE_RADIUS;
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
  }
  printf("[RT] number_of_vertices - 8: %ld, number_of_finite_edges: %ld\n",
         rt.number_of_vertices() - 8, rt.number_of_finite_edges());
  // no need to purge
  if (num_spheres == rt.number_of_vertices() - 8) return;

  // purge non-exist RT vertices (spheres)
  std::vector<MedialSphere> valid_medial_spheres;
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); vit++) {
    // skip 8 bbox points
    if (vit->info().all_id == -1) continue;
    // start purging
    int valid_id = valid_medial_spheres.size();
    Vertex_handle_rt vh = vit;
    int all_id = vh->info().all_id;
    assert(all_id != -1);
    MedialSphere new_msphere = all_medial_spheres.at(all_id);  // copy
    vh->info().all_id = valid_id;
    new_msphere.id = valid_id;
    valid_medial_spheres.push_back(new_msphere);
  }
  all_medial_spheres.clear();
  all_medial_spheres = valid_medial_spheres;
  printf("[RT] purged spheres %d->%ld, rt.number_of_vertices: %ld\n",
         num_spheres, valid_medial_spheres.size(), rt.number_of_vertices());
  assert(all_medial_spheres.size() == rt.number_of_vertices() - 8);
}

/**
 * @brief
 *
 * @param rt
 * @param n_site
 * @param site_knn
 * @return int return the site_k, maximum size of sphere neighbors
 */
int get_RT_vertex_neighbors(const RegularTriangulationNN& rt, const int& n_site,
                            std::vector<int>& site_knn) {
  // n_site == all_medial_spheres.size();
  // get sphere neighbors
  // sphere -> a set of neighbors
  std::vector<std::vector<int>> sphere_neighbors(n_site);
  std::vector<Vertex_handle_rt> one_neighs;
  int num_neigh_max = 0;
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); ++vit) {
    // skip 8 bbox points
    if (vit->info().all_id == -1) continue;
    one_neighs.clear();
    rt.finite_adjacent_vertices(vit, std::back_inserter(one_neighs));
    if (one_neighs.size() > num_neigh_max) {
      num_neigh_max = one_neighs.size();
    }

    int all_id = vit->info().all_id;
    for (auto& neigh_handle : one_neighs) {
      // skip 8 bbox points
      if (neigh_handle->info().all_id == -1) continue;
      sphere_neighbors[all_id].push_back(neigh_handle->info().all_id);
    }
  }

  // convert to flat 2D matrix
  // 1. dim: (num_neigh_max+1) x n_site
  // 2. each column j store all neighbors of sphere all_medial_spheres.at(j)
  // 3. init as -1
  site_knn.clear();
  site_knn.resize((num_neigh_max + 1) * n_site, -1);
  for (int all_id = 0; all_id < sphere_neighbors.size(); all_id++) {
    const auto& neighbors = sphere_neighbors.at(all_id);
    for (int neigh_id = 0; neigh_id < neighbors.size(); neigh_id++) {
      assert(neigh_id <= num_neigh_max);
      site_knn[neigh_id * n_site + all_id] = neighbors[neigh_id];
    }
  }

  // // for debug
  // printf("site_knn matrix: \n\t ");
  // for (uint i = 0; i < (num_neigh_max + 1); i++) {
  //   printf("\n\t ");
  //   // for (uint seed = 0; seed < n_site; seed++) {
  //   for (uint seed = 274; seed < 275; seed++) {
  //     printf("%d ", site_knn[i * n_site + seed]);
  //   }
  // }
  // printf("\n");

  // printf("maximum %d neighbors per sphere\n", num_neigh_max);
  return num_neigh_max;
}

// Each (non-restricted) powercell is dual to a RT tet
void get_PC_vertices(const RegularTriangulationNN& rt,
                     std::vector<Vector3>& pc_vertices) {
  pc_vertices.clear();
  for (Finite_cells_iterator_rt fci = rt.finite_cells_begin();
       fci != rt.finite_cells_end(); fci++) {
    RegularTriangulationNN::Bare_point bp =
        rt.dual(fci);  // dunno why we do not need *fci here
    pc_vertices.push_back(Vector3(bp[0], bp[1], bp[2]));
  }
  printf("found pc_vertices: %ld\n", pc_vertices.size());
}