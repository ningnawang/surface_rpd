#include "medial_sphere.h"

#include "assert.h"

/////////////////////////////////////////////////////////////////////////////////////
// TangentPlane
/////////////////////////////////////////////////////////////////////////////////////
TangentPlane::TangentPlane(const Vector3& _normal, const Vector3& _point,
                           const int _fid) {
  normal = _normal;
  points.push_back(_point);
  fid = _fid;
  energy = DBL_MAX;
  energy_over_sq_radius = DBL_MAX;
  is_deleted = false;
}

void TangentPlane::print_info() const {
  printf(
      "TangentPlane: is_deleted %d, fid: %d energy: %.10e, "
      "energy_over_sq_radius: "
      "%.10e, points[0] (%f,%f,%f), normal: (%f,%f,%f) , points size %ld\n",
      is_deleted, fid, energy, energy_over_sq_radius, points[0][0],
      points[0][1], points[0][2], normal[0], normal[1], normal[2],
      points.size());
}

bool TangentPlane::is_same_normal(const Vector3& bnormal,
                                  const double eps_degree) const {
  return is_vector_same_direction(normal, bnormal, eps_degree);
}

// 1. normals are simialr
// 2. tangent point is close
bool TangentPlane::is_same_tan_pl(const TangentPlane& tan_pl2) const {
  if (!is_same_normal(tan_pl2.normal)) {
    return false;
  }
  double dist = get_distance_between_two_vectors(points[0], tan_pl2.points[0]);
  if (dist > SCALAR_ZERO_1) return false;
  return true;
}

// static function
double TangentPlane::get_energy_value(const Vector3& theta,
                                      const double& radius, const double alpha1,
                                      const double alpha2,
                                      const Vector3& tan_point,
                                      const Vector3& normal) {
  Vector3 T = theta + radius * normal - tan_point;
  double K = GEO::dot(tan_point - theta, normal) - radius;
  return alpha1 * GEO::dot(T, T) + alpha2 * std::pow(K, 2);
}
// sphere (theta, radiu)
double TangentPlane::update_energy_value(const Vector3& theta,
                                         const double& radius,
                                         const double alpha1,
                                         const double alpha2) {
  energy = get_energy_value(theta, radius, alpha1, alpha2, points[0], normal);
  energy_over_sq_radius = energy / std::pow(radius, 2);
  return energy;
}

/////////////////////////////////////////////////////////////////////////////////////
// TangentConcaveLine
/////////////////////////////////////////////////////////////////////////////////////
TangentConcaveLine::TangentConcaveLine(const int _id, const FeatureEdge& fe) {
  id = _id;
  id_fe = fe.id;
  num_ce_group = fe.t2vs_group[3];
  energy = DBL_MAX;
}
// TangentConcaveLine::TangentConcaveLine(const int _id, const aint3
// _t2vs_group,
//                                        const avec2 _t2vs_pos,
//                                        const aint2& _sf_fs_pair,
//                                        const avec2& _adj_normals) {
//   id = _id;
//   t2vs_group = _t2vs_group;
//   t2vs_pos = _t2vs_pos;
//   direction = GEO::normalize(t2vs_pos[1] - t2vs_pos[0]);
//   tan_point = (t2vs_pos[0] + t2vs_pos[1]) / 2.;  // TODO: perturb?
//   is_tan_point_updated = false;
//   adj_sf_fs_pair = _sf_fs_pair;
//   adj_normals = _adj_normals;
//   get_random_vector_between_two_vectors(adj_normals[0], adj_normals[1],
//   normal); energy = DBL_MAX;
// }

// bool TangentConcaveLine::operator==(TangentConcaveLine const& b) const {
//   // belongs to the same grouped concave line
//   if (t2vs_group[3] == b.t2vs_group[3]) return true;
//   return false;

//   // Vector3 vec = (b.tan_point - tan_point).normalized();
//   // if (adj_sf_fs_pair == b.adj_sf_fs_pair) return true;

//   // // more likely to be different
//   // // adj_sf_fs_pair might be different, but represent the same line
//   // if (!is_vector_same_direction(direction, b.direction, esp_degree_5) &&
//   //     !is_vector_oppo_direction(direction, b.direction, esp_degree_5))
//   //   return false;

//   // if (!is_vector_same_direction(vec, direction, esp_degree_5) &&
//   //     !is_vector_oppo_direction(vec, direction, esp_degree_5))
//   //   return false;

//   // if (!is_vector_same_direction(vec, b.direction, esp_degree_5) &&
//   //     !is_vector_oppo_direction(vec, b.direction, esp_degree_5))
//   //   return false;

//   // // vec should parallel to both direction and b.direction
//   // return true;
// }

// void TangentConcaveLine::print_info() const {
//   printf(
//       "TangentConcaveLine: energy: %f, energy_over_sq_radius: %f, "
//       "adj_sf_fs_pair: (%d,%d), tan_point (%f,%f,%f), normal: "
//       "(%f,%f,%f) \n",
//       energy, energy_over_sq_radius, adj_sf_fs_pair[0], adj_sf_fs_pair[1],
//       tan_point[0], tan_point[1], tan_point[2], normal[0], normal[1],
//       normal[2]);
//   printf("adj_normals (%f,%f,%f) - (%f,%f,%f)", adj_normals[0][0],
//          adj_normals[0][1], adj_normals[0][2], adj_normals[1][0],
//          adj_normals[1][1], adj_normals[1][2]);
// }

// // if concave line is a curve
// // then each line should cover more adjacent faces
// bool TangentConcaveLine::is_normal_covered_by_adj_fs(
//     const Vector3& n, double esp_degree_given) const {
//   if (is_vector_same_direction(adj_normals[0], n, esp_degree_given))
//     return true;
//   if (is_vector_same_direction(adj_normals[1], n, esp_degree_given))
//     return true;
//   return false;
// }

// given two points of a segment
void TangentConcaveLine::update_tangent_point(const GEO::vec3& v1_p,
                                              const GEO::vec3& v2_p) {
  GEO::vec3 v_mid = 1. / 2. * (v1_p + v2_p);
  Vector3 p2(v_mid[0], v_mid[1], v_mid[2]);
  if (!is_tan_point_updated)
    tan_point = p2;
  else
    tan_point = 1. / 2. * (p2 + tan_point);
}

// static function
double TangentConcaveLine::get_energy_value(const Vector3& theta,
                                            const double& radius,
                                            const double alpha3,
                                            const Vector3& tan_point,
                                            const Vector3& normal) {
  Vector3 T = theta + radius * normal - tan_point;
  return alpha3 * GEO::dot(T, T);
}
// sphere (theta, radiu)
double TangentConcaveLine::update_energy_value(const Vector3& theta,
                                               const double& radius,
                                               const double alpha3) {
  // sanity_check();
  energy = get_energy_value(theta, radius, alpha3, tan_point, normal);
  energy_over_sq_radius = energy / std::pow(radius, 2);
  return energy;
}

/////////////////////////////////////////////////////////////////////////////////////
// MedialSphere
/////////////////////////////////////////////////////////////////////////////////////
MedialSphere::MedialSphere(int _id, Vector3 _pin, Vector3 _pin_normal,
                           SphereType _type) {
  id = _id;
  ss.p = _pin;
  ss.p_normal = _pin_normal;
  radius = INIT_RADIUS;
  type = _type ? _type : SphereType::T_2;
  center = ss.p - ss.p_normal * radius;
  pcell.topo_status = Topo_Status::unkown;
}
MedialSphere::MedialSphere(int _id, Vector3 _center, double _radius,
                           SphereType _type) {
  id = _id;
  center = _center;
  radius = _radius;
  type = _type;
  pcell.topo_status = Topo_Status::unkown;
}

void MedialSphere::clear_mm() {
  this->edges_.clear();
  this->faces_.clear();
}

void MedialSphere::print_info() const {
  printf("------ MedialSphere Info ------\n");
  printf(
      "id: %d, is_deleted: %d, is_on_se: %d, is_on_corner: %d, is_on_ce_pin: "
      "%d, is_on_ce: %d \n",
      id, is_deleted, is_on_se(), is_on_corner(), is_on_ce_pin(), is_on_ce());
  printf("center: (%f,%f,%f), radius: %f, type: %d \n", center[0], center[1],
         center[2], radius, type);
  printf("ss_params: p: (%f,%f,%f), q: (%f,%f,%f), q_fid: %d, p_fid: %d \n",
         ss.p[0], ss.p[1], ss.p[2], ss.q[0], ss.q[1], ss.q[2], ss.q_fid,
         ss.p_fid);
  printf("tan_planes: %ld \n", tan_planes.size());
  printf("pcell: topo_status: %d\n", pcell.topo_status);
  print_tan_planes();
}

double MedialSphere::update_all_energy_values(const double alpha1,
                                              const double alpha2) {
  double sum = 0.;
  for (auto& tan_pl : tan_planes) {
    if (tan_pl.is_deleted) continue;
    tan_pl.update_energy_value(center, radius, alpha1, alpha2);
    sum += tan_pl.energy;
  }
  return sum;
}

void MedialSphere::print_tan_planes() const {
  printf("------- MVertex %d Tangent Planes: %ld\n", id, tan_planes.size());
  for (const auto& tan_pl : tan_planes) {
    tan_pl.print_info();
  }
}

void MedialSphere::save_old_center_radius(bool is_clear) {
  old_center = center;
  old_radius = radius;
  if (is_clear) {
    center = Vector3(0.0, 0.0, 0.0);
    radius = 0.;
  }
}

void MedialSphere::update_tan_planes_from_ss_params() {
  TangentPlane tan_pl1(ss.p_normal, ss.p, ss.p_fid);
  TangentPlane tan_pl2(ss.q_normal, ss.q, ss.q_fid);
  tan_planes.push_back(tan_pl1);
  tan_planes.push_back(tan_pl2);
}

void MedialSphere::new_tan_plane_no_dup(const Vector3& _normal,
                                        const Vector3& _point, const int _fid) {
  // if (id == 11) print_tan_planes();
  TangentPlane tan_pl(_normal, _point, _fid);
  // if (id == 11) {
  //   printf("+++");
  //   tan_pl.print_info();
  // }
  bool is_add = true;
  for (const auto& curr_tan_pl : tan_planes) {
    if (curr_tan_pl.is_same_tan_pl(tan_pl)) {
      // if (id == 11) printf("duplicated!! not add \n");
      is_add = false;
      break;
    }
  }
  if (!is_add) return;
  tan_planes.push_back(tan_pl);
}

void MedialSphere::new_cc_line_no_dup(const FeatureEdge& one_fe) {
  TangentConcaveLine one_cc_line(tan_cc_lines.size(), one_fe);
  bool is_added = true;
  for (const auto& curr_cc_line : tan_cc_lines) {
    if (curr_cc_line.num_ce_group == one_cc_line.num_ce_group ||
        curr_cc_line.id_fe == one_cc_line.id_fe) {
      is_added = false;
      break;
    }
  }

  if (!is_added) return;
  tan_cc_lines.push_back(one_cc_line);
}

// perturbation of the coordinates
// to avoid boundary error (degenerations) while calculating RPD
void MedialSphere::apply_perturb() {
  std::srand(RAN_SEED);
  int i = RANDOM_INT(0, 3);
  center[i] += SCALAR_ZERO_2;
}

void MedialSphere::pcell_insert(int cell_id) { pcell.cell_ids.insert(cell_id); }

void MedialSphere::topo_clear() {
  // cells
  pcell.cell_ids.clear();
  pcell.cell_neighbors.clear();
  pcell.cell_to_surfv2fid.clear();
  pcell.cc_cells.clear();
  pcell.cc_bfids.clear();
  pcell.cc_surf_v2fids.clear();
  // facets
  pcell.f_id1_to_cells.clear();
  pcell.f_id2_to_cells.clear();
  pcell.cell_to_f_id1s.clear();
  pcell.facet_cc_cells.clear();
  pcell.facet_cc_surf_v2fids.clear();
  pcell.f_id2_is_fixed.clear();
  // powercell edges
  pcell.e_to_cells.clear();
  pcell.e_is_fixed.clear();
  pcell.edge_cc_cells.clear();
  pcell.edge_2endvertices.clear();
  // powercell vertices
  pcell.vertex_2id.clear();
  pcell.vertex_2pos.clear();
  // sharp edges
  pcell.se_covered_lvids.clear();
  // surface faces
  pcell.surf_v2fid_in_groups.clear();

  pcell.topo_status = Topo_Status::ok;
  euler = 0;
  euler_sum = 0;
  num_cells = 0;
}

void MedialSphere::fcc_fixed(int neigh_id) {
  // only update when the topo_type is the same
  // if (pcell.topo_status == Topo_Status::high_facet_cc)
  // assert(pcell.f_id2_is_fixed.find(neigh_id) != pcell.f_id2_is_fixed.end());
  pcell.f_id2_is_fixed[neigh_id] = true;
}

// for removing medial spheres that are too close
// using absolute value is not reliable for all models
// let's use ratio to measure how deep two spheres intersect
bool MedialSphere::operator==(const MedialSphere& m2) const {
  // our mesh has been normalized through
  // MeshIO::normalize_mesh()
  double threshold = SCALAR_ZERO_2;
  double diff = (center - m2.center).length();
  if (diff > threshold) return false;
  diff = std::abs(radius - m2.radius);
  if (diff > threshold) return false;
  return true;
}

bool MedialSphere::operator!=(const MedialSphere& m2) const {
  // if (id != m2.id) return true;
  double diff = (center - m2.center).length();
  if (diff > SCALAR_ZERO_2) return true;
  diff = std::abs(radius - m2.radius);
  if (diff > SCALAR_ZERO_2) return true;
  return false;
}

bool MedialSphere::is_on_se() const {
  if (type == SphereType::T_1_2) return true;
  return false;
}
bool MedialSphere::is_on_corner() const {
  if (type == SphereType::T_1_N) return true;
  return false;
}
bool MedialSphere::is_on_ce_pin() const {
  if (type == SphereType::T_c) return true;
  return false;
}
bool MedialSphere::is_on_ce() const {
  if (type == SphereType::T_2_c || type == SphereType::T_N_c) return true;
  return false;
}
bool MedialSphere::is_on_intf() const {
  if (type == SphereType::T_N) return true;
  return false;
}
bool MedialSphere::is_on_extf() const {
  if (is_on_se()) return true;
  if (is_on_corner()) return true;
  return false;
}

/////////////////////////////////////////////////////////////////////////////////////
// Other Functions
/////////////////////////////////////////////////////////////////////////////////////
bool validate_new_sphere(const std::vector<MedialSphere>& all_medial_spheres,
                         const MedialSphere& new_sphere) {
  // make sure new_sphere is not duplicated with any sphere
  for (const auto& msphere : all_medial_spheres) {
    if (msphere == new_sphere) {
      printf(
          "[NewSphereAdd Failed] new_sphere "
          "(%f,%f,%f,%f) too close to sphere "
          "%d (%f,%f,%f,%f), not add\n",
          new_sphere.center[0], new_sphere.center[1], new_sphere.center[2],
          new_sphere.radius, msphere.id, msphere.center[0], msphere.center[1],
          msphere.center[2], msphere.radius);
      return false;
    }
  }
  return true;
}

bool add_new_sphere_validate(std::vector<MedialSphere>& all_medial_spheres,
                             MedialSphere& new_sphere) {
  if (!validate_new_sphere(all_medial_spheres, new_sphere)) return false;
  // save to add
  new_sphere.id = all_medial_spheres.size();
  printf("[NewSphereAdd Success] new_sphere added %d\n", new_sphere.id);
  all_medial_spheres.push_back(new_sphere);
  return true;
}

// This function will update MedialSphere::id
void purge_deleted_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres) {
  std::vector<MedialSphere> purged;
  int new_id = 0, num_deleted = 0;
  for (const auto& msphere : all_medial_spheres) {
    if (msphere.is_deleted) {
      num_deleted++;
      continue;
    }
    purged.push_back(msphere);
    purged.back().id = new_id++;
  }

  if (num_deleted == 0) return;
  all_medial_spheres.clear();
  all_medial_spheres = purged;
  printf("[Purge] purged %d/%ld deleted spheres \n", num_deleted,
         all_medial_spheres.size());
}

// check if two medial spheres on the same sharp edge
bool is_two_mspheres_on_same_se(const MedialSphere& msphere1,
                                const MedialSphere& msphere2) {
  if (!msphere1.is_on_se() || !msphere2.is_on_se()) return false;
  std::set<int> s1_se_group, s2_se_group;
  for (const aint4& se_lvds : msphere1.pcell.se_covered_lvids)
    s1_se_group.insert(se_lvds[3]);
  for (const aint4& se_lvds : msphere2.pcell.se_covered_lvids)
    s2_se_group.insert(se_lvds[3]);
  if (s1_se_group != s2_se_group) {
    return false;
  }
  return true;
}