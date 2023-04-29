#ifndef H_IO_H
#define H_IO_H

#include <assert.h>
#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "common_cxx.h"
#include "medial_sphere.h"
#include "params.h"

void get_bbox(const std::vector<float>& vertices, float& xmin, float& ymin,
              float& zmin, float& xmax, float& ymax, float& zmax,
              float& bbox_diag_l);

bool load_tet(const std::string& filename, std::vector<float>& vertices,
              std::vector<int>& indices, bool normalize, Parameter& params);

void load_spheres_to_sites_given(std::vector<MedialSphere>& all_medial_spheres,
                                 bool& site_is_transposed,
                                 std::vector<float>& site,
                                 std::vector<float>& site_weights, int& n_site);

void load_spheres_from_file(const char* filename,
                            std::vector<MedialSphere>& all_medial_spheres,
                            bool is_load_type);

void load_sites_from_file(bool& site_is_transposed, std::vector<float>& site,
                          std::vector<float>& site_weights, int& n_site,
                          const char* filename);
void save_spheres_file(const std::vector<MedialSphere>& all_medial_spheres,
                       const std::string filename, bool is_save_type);

void load_v2tets(const std::vector<float>& vertices,
                 const std::vector<int>& indices,
                 std::map<int, std::set<int>>& v2tets);

void load_surface_vertices(const std::vector<float>& vertices,
                           const std::vector<int>& indices,
                           std::vector<float>& surf_vertices);

bool load_surface_mesh(const std::string& path, GEO::Mesh& input);
bool load_surface_mesh_geogram(const std::string& path, GEO::Mesh& input);

void load_sf_mesh_from_internal(const std::vector<std::array<float, 3>>& points,
                                const std::vector<std::array<int, 3>>& faces,
                                const std::vector<int>& input_vs_id_attr,
                                GEO::Mesh& input);

// void write_convex_cells(std::vector<float3>& voro_points,
//                         std::vector<std::array<int, 4>>& voro_tets,
//                         std::vector<int>& voro_tets_sites);

void get_surface_from_tet(const std::vector<float>& tet_vertices,
                          const std::vector<int>& tet_indices,
                          std::vector<std::array<float, 3>>& surf_vertices,
                          std::vector<std::array<int, 3>>& surf_faces,
                          std::vector<int>& sf_vs_2_tet_vs);

bool save_sf_mesh(const std::string sf_path, const GEO::Mesh& sf_mesh);
bool save_sf_mesh_geogram(const std::string sf_path, GEO::Mesh& sf_mesh);

#endif  // __IO_H__
