#include "input_types.h"

void SurfaceMesh::reload_sf2tet_vs_mapping() {
  // load sf2tet_vs_mapping
  const GEO::Attribute<int> tet_vid_attr(this->vertices.attributes(),
                                         "tet_vid");
  sf2tet_vs_mapping.clear();
  sf2tet_vs_mapping.resize(this->vertices.nb());
  for (uint v = 0; v < this->vertices.nb(); v++) {
    sf2tet_vs_mapping[v] = tet_vid_attr[v];
    // if (v == 214 || v == 224) {
    //   printf(
    //       "!!!!! surface vid %d has sf2tet_vs_mapping[i] %d, tet_vid_attr[i]:
    //       "
    //       "%d \n",
    //       v, sf2tet_vs_mapping[v], tet_vid_attr[v]);
    // }
  }
}

/**
 * tet2sf_vs_mapping:   id mapping from tet_vertices to sf_mesh
 * sf_vs2fids:          sf_mesh vertices to adjacent facets
 * tet_vs2sf_fids:      tet vertex to sf_mesh fids
 */
void load_sf_tet_mapping(const GEO::Mesh& sf_mesh,
                         std::map<int, int>& tet2sf_vs_mapping,
                         std::map<int, std::set<int>>& sf_vs2fids,
                         std::map<int, std::set<int>>& tet_vs2sf_fids) {
  tet2sf_vs_mapping.clear();
  sf_vs2fids.clear();
  const GEO::Attribute<int> tet_vid_attr(sf_mesh.vertices.attributes(),
                                         "tet_vid");
  for (uint v = 0; v < sf_mesh.vertices.nb(); v++) {
    tet2sf_vs_mapping[tet_vid_attr[v]] = v;
  }

  for (uint f = 0; f < sf_mesh.facets.nb(); f++) {
    int f_nb_v = sf_mesh.facets.nb_vertices(f);
    assert(f_nb_v == 3);
    for (uint lv = 0; lv < f_nb_v; lv++) {
      uint vid = sf_mesh.facets.vertex(f, lv);
      uint tvid = tet_vid_attr[vid];
      sf_vs2fids[vid].insert(f);
      tet_vs2sf_fids[tvid].insert(f);
    }
  }
}

// Given fid_given, we fetch the k_ring neighboring faces
// not crossing sharp edges
void get_k_ring_neighbors_no_cross(const GEO::Mesh& sf_mesh,
                                   const std::set<aint2>& fe_sf_pairs_not_cross,
                                   const int fid_given, const int k,
                                   std::set<int>& k_ring_fids, bool is_debug) {
  auto is_skip_se_neighbor = [&](const int f, const int nf) {
    std::array<int, 2> ref_fs_pair = {{f, nf}};
    std::sort(ref_fs_pair.begin(), ref_fs_pair.end());
    if (fe_sf_pairs_not_cross.find(ref_fs_pair) !=
        fe_sf_pairs_not_cross.end()) {
      if (is_debug)
        printf("[K_RING] face %d skip nf %d since sharing a sharp edge \n", f,
               nf);
      return true;  // skip checking its neighbor
    }
    return false;
  };

  assert(fid_given >= 0);
  if (is_debug)
    printf("calling get_k_ring_neighbor_no_se for fid_given: %d \n", fid_given);

  k_ring_fids.clear();
  k_ring_fids.insert(fid_given);
  std::set<int> new_added_fids, new_added_fids_copy;
  new_added_fids.insert(fid_given);

  int i = 0;
  while (i < k) {
    new_added_fids_copy = new_added_fids;
    new_added_fids.clear();
    for (const auto& fid : new_added_fids_copy) {
      for (GEO::index_t le = 0; le < sf_mesh.facets.nb_vertices(fid); le++) {
        GEO::index_t nfid = sf_mesh.facets.adjacent(fid, le);
        if (nfid == GEO::NO_FACET) continue;
        if (is_skip_se_neighbor(fid, nfid)) continue;
        if (is_debug) printf("fid %d has neighbor face nfid %d\n", fid, nfid);
        k_ring_fids.insert(nfid);
        new_added_fids.insert(nfid);
      }  // for facets.nb_vertices
    }    // for k_ring_fids
    i++;
  }  // while (i < k)
  if (is_debug) {
    printf("[K_RING] fid %d found k_ring %d neighbors: %ld\n", fid_given, k,
           k_ring_fids.size());
  }
}

void store_special_edge(const TetMesh& tet_mesh, const SurfaceMesh& sf_mesh,
                        const EdgeType& fe_type, const aint3& t2vs_group,
                        std::vector<FeatureEdge>& feature_edges) {
  int tv1 = t2vs_group[0];
  int tv2 = t2vs_group[1];
  assert(tet_mesh.tet_vs2sf_fids.find(tv1) != tet_mesh.tet_vs2sf_fids.end());
  assert(tet_mesh.tet_vs2sf_fids.find(tv2) != tet_mesh.tet_vs2sf_fids.end());

  Vector3 tv1_pos(tet_mesh.tet_vertices.at(tv1 * 3),
                  tet_mesh.tet_vertices.at(tv1 * 3 + 1),
                  tet_mesh.tet_vertices.at(tv1 * 3 + 2));
  Vector3 tv2_pos(tet_mesh.tet_vertices.at(tv2 * 3),
                  tet_mesh.tet_vertices.at(tv2 * 3 + 1),
                  tet_mesh.tet_vertices.at(tv2 * 3 + 2));

  std::vector<int> adj_sf_fids;
  set_intersection<int>(tet_mesh.tet_vs2sf_fids.at(tv1),
                        tet_mesh.tet_vs2sf_fids.at(tv2), adj_sf_fids);
  assert(adj_sf_fids.size() == 2);
  std::sort(adj_sf_fids.begin(), adj_sf_fids.end());

  FeatureEdge fe(feature_edges.size(), fe_type, t2vs_group);
  fe.t2vs_pos = {{tv1_pos, tv2_pos}};
  fe.adj_sf_fs_pair = {{adj_sf_fids[0], adj_sf_fids[1]}};
  fe.adj_normals = {{get_mesh_facet_normal(sf_mesh, adj_sf_fids[0]),
                     get_mesh_facet_normal(sf_mesh, adj_sf_fids[1])}};
  fe.adj_tan_points = {{get_mesh_facet_centroid(sf_mesh, adj_sf_fids[0]),
                        get_mesh_facet_centroid(sf_mesh, adj_sf_fids[1])}};
  feature_edges.push_back(fe);
}