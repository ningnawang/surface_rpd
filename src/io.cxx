
#include "io.h"

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>

#include <bitset>

#include "../extern/mshloader/MshLoader.h"

void get_bbox(const std::vector<float>& vertices, float& xmin, float& ymin,
              float& zmin, float& xmax, float& ymax, float& zmax,
              float& bbox_diag_l) {
  int nb_v = vertices.size() / 3;
  xmin = xmax = vertices[0];
  ymin = ymax = vertices[1];
  zmin = zmax = vertices[2];
  for (int i = 1; i < nb_v; ++i) {
    xmin = std::min(xmin, vertices[3 * i]);
    ymin = std::min(ymin, vertices[3 * i + 1]);
    zmin = std::min(zmin, vertices[3 * i + 2]);
    xmax = std::max(xmax, vertices[3 * i]);
    ymax = std::max(ymax, vertices[3 * i + 1]);
    zmax = std::max(zmax, vertices[3 * i + 2]);
  }
  float d = xmax - xmin;
  d = std::max(d, ymax - ymin);
  d = std::max(d, zmax - zmin);
  d = 0.001f * d;
  xmin -= d;
  ymin -= d;
  zmin -= d;
  xmax += d;
  ymax += d;
  zmax += d;

  bbox_diag_l = std::sqrt(std::pow(xmax - xmin, 2) + std::pow(ymax - ymin, 2) +
                          std::pow(zmax - zmin, 2));
}

// vertices: tet vertices
// indices: tet 4 indices of vertices
using namespace PyMesh;
bool load_tet(const std::string& filename, std::vector<float>& vertices,
              std::vector<int>& indices, bool normalize, Parameter& params) {
  std::string s;
  int n_vertex, n_tet, temp;

  std::ifstream input(filename);
  if (input.fail()) return false;

  std::string ext = filename.substr(filename.find_last_of('.') + 1);
  if (ext == "msh") {
    MshLoader msh_loader(filename);
    vertices = msh_loader.get_nodes();
    indices = msh_loader.get_elements();
    n_vertex = vertices.size() / 3;
  } else if (ext == "tet") {
    input >> n_vertex >> n_tet;
    vertices.resize(3 * n_vertex);
    indices.resize(n_tet << 2);

    for (int i = 0; i < n_vertex; ++i)
      input >> vertices[3 * i] >> vertices[3 * i + 1] >> vertices[3 * i + 2];

    for (int i = 0; i < n_tet; ++i) {
      input >> temp >> indices[(i << 2)] >> indices[(i << 2) + 1] >>
          indices[(i << 2) + 2] >> indices[(i << 2) + 3];
      assert(temp == 4);
    }
  } else if (ext == "vtk") {
    for (int i = 0; i < 4; ++i) std::getline(input, s);  // skip first 4 lines

    input >> s >> n_vertex >> s;
    vertices.resize(3 * n_vertex);
    for (int i = 0; i < n_vertex; ++i)
      input >> vertices[3 * i] >> vertices[3 * i + 1] >> vertices[3 * i + 2];

    input >> s >> n_tet >> s;
    indices.resize(n_tet << 2);
    for (int i = 0; i < n_tet; ++i) {
      // A single left shift multiplies a binary number by 2:
      input >> temp >> indices[(i << 2)] >> indices[(i << 2) + 1] >>
          indices[(i << 2) + 2] >> indices[(i << 2) + 3];
      assert(temp == 4);
      for (uint j = 0; j < 4; ++j) --indices[(i << 2) + j];
    }
  } else {
    input.close();
    return false;
  }
  input.close();
  std::cout << "loaded tet_mesh #v: " << n_vertex
            << ", #t: " << indices.size() / 4 << std::endl;

  // normalize vertices between [0,1000]^3
  float xmin, ymin, zmin, xmax, ymax, zmax;
  get_bbox(vertices, xmin, ymin, zmin, xmax, ymax, zmax, params.bbox_diag_l);
  if (normalize) {
    float maxside = std::max(std::max(xmax - xmin, ymax - ymin), zmax - zmin);
    // #pragma omp parallel for
    for (int i = 0; i < n_vertex; i++) {
      vertices[3 * i] = params.scale_max * (vertices[3 * i] - xmin) / maxside;
      vertices[3 * i + 1] =
          params.scale_max * (vertices[3 * i + 1] - ymin) / maxside;
      vertices[3 * i + 2] =
          params.scale_max * (vertices[3 * i + 2] - zmin) / maxside;
    }
    get_bbox(vertices, xmin, ymin, zmin, xmax, ymax, zmax, params.bbox_diag_l);
    std::cerr << "bbox [" << xmin << ":" << xmax << "], [" << ymin << ":"
              << ymax << "], [" << zmin << ":" << zmax
              << "], bbox_diag_l: " << params.bbox_diag_l << std::endl;
  }

  // update 8 bbox points
  Vector3 pmin(xmin, ymin, zmin);
  Vector3 pmax(xmax, ymax, zmax);
  for (int i = 0; i < 8; i++) {
    std::array<double, 3> p;
    std::bitset<sizeof(int) * 8> flag(i);
    for (int j = 0; j < 3; j++) {
      if (flag.test(j))
        p[j] = pmax[j];
      else
        p[j] = pmin[j];
    }
    params.bb_points.push_back(p[0]);
    params.bb_points.push_back(p[1]);
    params.bb_points.push_back(p[2]);
  }

  return true;
}

void load_spheres_to_sites_given(std::vector<MedialSphere>& all_medial_spheres,
                                 bool& site_is_transposed,
                                 std::vector<float>& site,
                                 std::vector<float>& site_weights,
                                 int& n_site) {
  std::vector<std::array<float, 4>> spheres;
  // spheres.push_back({{-2.251256, -1.251256, 1.543551, 1.543575}});
  // spheres.push_back({{0.000000, -3.000000, 1.684660, 1.684660}});
  // spheres.push_back({{0.000000, 3.000000, 2.552344, 2.552344}});

  // spheres.push_back({{-2.5, -1.5, 1.5, 1.0}});
  // spheres.push_back({{0.0, -3.0, 1.5, 1.0}});
  // spheres.push_back({{0.0, 3.0, 2.5, 1.0}});

  // spheres.push_back({{-2.251262, -1.251262, 1.543560, 1.543562}});
  // spheres.push_back({{0.000000, -3.000000, 1.684660, 1.684660}});
  // spheres.push_back({{2.251262, -1.251262, 1.543560, 1.643562}});

  spheres.push_back({{1.944161, -12.897811, -2.092603, 0.528340}});
  spheres.push_back({{-3.450441, 13.444244, -0.001805, 0.891773}});

  n_site = spheres.size();
  int dim = 4;
  std::cout << "n_site: " << n_site << ", dim: " << dim << std::endl;

  all_medial_spheres.clear();
  for (int i = 0; i < n_site; i++) {
    MedialSphere msphere(all_medial_spheres.size(), Vector3(0, 0, 0),
                         Vector3(0, 0, 0));
    msphere.center = Vector3(spheres[i][0], spheres[i][1], spheres[i][2]);
    msphere.radius = spheres[i][3];
    all_medial_spheres.push_back(msphere);
  }

  // use knn_cuda_global_dev if true (n_site is small), else knearests
  site_is_transposed = n_site < KNEARESTS_MIN_N;
  site.resize(n_site * 3);
  site_weights.resize(n_site);
  for (int i = 0; i < n_site; ++i) {
    const auto& msphere = all_medial_spheres.at(i);
    if (site_is_transposed) {
      site[i] = msphere.center[0];
      site[i + n_site] = msphere.center[1];
      site[i + (n_site << 1)] = msphere.center[2];
    } else {
      site[3 * i] = msphere.center[0];
      site[3 * i + 1] = msphere.center[1];
      site[3 * i + 2] = msphere.center[2];
    }
    // add weight (sq_radius) info
    if (dim == 4) {
      site_weights[i] = std::pow(float(msphere.radius), 2);
    } else
      site_weights[i] = 1;
    printf("site %d has sphere: (%lf %lf %lf %lf) \n", i, msphere.center[0],
           msphere.center[1], msphere.center[2], msphere.radius);
  }
}

void load_spheres_from_file(const char* filename,
                            std::vector<MedialSphere>& all_medial_spheres,
                            bool is_load_type) {
  std::ifstream file(filename);
  int dim, n_site;
  int type;
  file >> dim >> n_site;
  assert(dim == 3 || dim == 4);
  std::cout << "n_site: " << n_site << ", dim: " << dim << std::endl;

  all_medial_spheres.clear();
  for (int i = 0; i < n_site; i++) {
    MedialSphere msphere(all_medial_spheres.size(), Vector3(0, 0, 0),
                         Vector3(0, 0, 0));
    file >> msphere.center[0] >> msphere.center[1] >> msphere.center[2];
    if (dim == 4)
      file >> msphere.radius;
    else
      msphere.radius = 1;
    if (is_load_type) {
      file >> type;
      if (type != SphereType::T_UNK)
        msphere.type = SphereType(type);  // else as T2 sphere
    }
    all_medial_spheres.push_back(msphere);
    printf("site %d has sphere: (%lf %lf %lf %lf), type: %d\n", i,
           msphere.center[0], msphere.center[1], msphere.center[2],
           msphere.radius, msphere.type);
  }
  file.close();
}

void load_sites_from_file(bool& site_is_transposed, std::vector<float>& site,
                          std::vector<float>& site_weights, int& n_site,
                          const char* filename) {
  std::ifstream file(filename);
  int dim = -1;
  file >> dim >> n_site;
  assert(dim == 3 || dim == 4);
  std::cout << "n_site: " << n_site << ", dim: " << dim << std::endl;

  // use knn_cuda_global_dev if true (n_site is small), else knearests
  site_is_transposed = n_site < KNEARESTS_MIN_N;
  site.resize(n_site * 3);
  site_weights.resize(n_site);
  for (int i = 0; i < n_site; ++i) {
    if (site_is_transposed)
      file >> site[i] >> site[i + n_site] >> site[i + (n_site << 1)];
    else
      file >> site[3 * i] >> site[3 * i + 1] >> site[3 * i + 2];
    // add weight (sq_radius) info
    if (dim == 4) {
      float radius;
      file >> radius;
      site_weights[i] = std::pow(radius, 2);
    } else
      site_weights[i] = 1;

    std::cout << "site " << i << " has weight: " << site_weights[i]
              << std::endl;
  }
  file.close();
}

void save_spheres_file(const std::vector<MedialSphere>& all_medial_spheres,
                       const std::string filename, bool is_save_type) {
  std::string sphere_path =
      "../out/sph/sph_" + filename + "_" + get_timestamp() + ".sph";
  int n_site = all_medial_spheres.size();
  std::fstream file;
  file.open(sphere_path, std::ios_base::out);
  file << 4 << " " << n_site << std::endl;
  for (int i = 0; i < n_site; i++) {
    const auto& msphere = all_medial_spheres.at(i);
    file << msphere.center[0] << " " << msphere.center[1] << " "
         << msphere.center[2] << " " << msphere.radius;
    if (is_save_type) file << " " << msphere.type;
    file << std::endl;
  }
  file.close();
  printf("saved .sph file %s\n", sphere_path.c_str());
}

void load_v2tets(const std::vector<float>& vertices,
                 const std::vector<int>& indices,
                 std::map<int, std::set<int>>& v2tets) {
  v2tets.clear();
  int n_tets = indices.size() / 4;
  for (int t = 0; t < n_tets; t++) {
    for (uint lv = 0; lv < 4; lv++) {
      v2tets[indices[t * 4 + lv]].insert(t);
    }
  }
  // sanity
  if (v2tets.size() != vertices.size() / 3) {
    std::cerr << "ERROR: vertices size / 3: " << vertices.size() / 3
              << " not equal to v2tets size: " << v2tets.size() << std::endl;
    std::cout << "indices size / 4: " << indices.size() / 4 << std::endl;
    exit(1);
  }
};

void load_surface_vertices(const std::vector<float>& vertices,
                           const std::vector<int>& indices,
                           std::vector<float>& surf_vertices) {
  uint dim = 4;
  assert(indices.size() % dim == 0);

  // facet -> count, surface facet only count 1, otherwise 2
  std::map<std::array<int, 3>, int> map_fcnt;  // facet count
  for (uint i = 0; i < indices.size() / dim; i++) {
    uint index = i * dim;
    for (uint j = 0; j < dim; j++) {
      std::array<int, 3> facet = {{indices.at(index + j),
                                   indices.at(index + (j + 1) % dim),
                                   indices.at(index + (j + 2) % dim)}};

      std::sort(facet.begin(), facet.end());
      if (map_fcnt.find(facet) == map_fcnt.end()) {
        map_fcnt[facet] = 0;
      }
      map_fcnt[facet]++;
    }
  }

  std::set<int> surf_vset;
  for (auto& pair : map_fcnt) {
    if (pair.second != 1) continue;  // not on surface
    auto& facet = pair.first;
    for (int i = 0; i < 3; i++) surf_vset.insert(facet[i]);
  }
  printf("[Surface] found %lu surface vertices\n", surf_vset.size());

  surf_vertices.clear();
  for (uint vid : surf_vset) {
    surf_vertices.push_back(vertices[vid]);
  }
}

bool load_surface_mesh(const std::string& path, GEO::Mesh& input) {
  std::cout << "Loading mesh at" << path << std::endl;
  input.clear(false, false);
  const bool ok = GEO::mesh_load(path, input);
  std::cout << "ok: " << ok << std::endl;
  return ok;
}

/**
 * Here we save the vertex's old id in tet_vertices as attributes
 */
bool load_surface_mesh_geogram(const std::string& path, GEO::Mesh& input) {
  std::cout << "Loading mesh at" << path << std::endl;
  if (get_file_ext(path) != "geogram") {
    printf("Please use mesh format as .geogram");
    return false;
  }
  input.clear(false, false);
  GEO::InputGeoFile geo_file(path);
  GEO::MeshIOFlags flags;
  flags.set_element(GEO::MeshElementsFlags::MESH_ALL_ELEMENTS);
  flags.set_attributes(GEO::MeshAttributesFlags::MESH_ALL_ATTRIBUTES);
  const bool ok = GEO::mesh_load(geo_file, input, flags);
  if (!ok) return false;
  return ok;
}

/**
 * sf2tet_vs_mapping: mapping to tet_vertices indices
 */
void load_sf_mesh_from_internal(const std::vector<std::array<float, 3>>& points,
                                const std::vector<std::array<int, 3>>& faces,
                                const std::vector<int>& sf2tet_vs_mapping,
                                GEO::Mesh& input) {
  std::cout << "Loading mesh from internal data ..." << std::endl;
  input.clear(false, false);
  GEO::Attribute<int> tet_vid_attr(input.vertices.attributes(), "tet_vid");

  // Setup vertices
  input.vertices.create_vertices(points.size());
  for (uint i = 0; i < input.vertices.nb(); ++i) {
    GEO::vec3& p = input.vertices.point(i);
    p[0] = points[i][0];
    p[1] = points[i][1];
    p[2] = points[i][2];
    tet_vid_attr[i] = sf2tet_vs_mapping[i];
  }

  // Setup faces
  input.facets.create_triangles(faces.size());
  for (uint c = 0; c < input.facets.nb(); ++c) {
    for (uint lv = 0; lv < 3; ++lv) {
      input.facets.set_vertex(c, lv, faces[c][lv]);
    }
  }

  // Setup edges
  std::vector<aint2> edges;
  for (uint i = 0; i < faces.size(); i++) {
    const auto& f = faces[i];
    for (uint j = 0; j < 3; j++) {
      aint2 e = {{f[j], f[(j + 1) % 3]}};
      if (e[0] > e[1]) std::swap(e[0], e[1]);
      edges.push_back(e);
    }
  }
  vector_unique(edges);
  input.edges.create_edges(edges.size());
  for (uint e = 0; e < edges.size(); e++) {
    for (uint lv = 0; lv < 2; ++lv) {
      input.edges.set_vertex(e, lv, edges[e][lv]);
    }
  }

  GEO::mesh_reorder(input, GEO::MESH_ORDER_MORTON);
  // we did not setup the adjacent info till now
  input.facets.connect();
}

bool save_sf_mesh(const std::string sf_path, const GEO::Mesh& sf_mesh) {
  bool ok = GEO::mesh_save(sf_mesh, sf_path);
  std::cout << "saving sf_mesh ok: " << ok << std::endl;
  return ok;
}

bool save_sf_mesh_geogram(const std::string sf_path, GEO::Mesh& sf_mesh) {
  GEO::OutputGeoFile geo_file(sf_path);
  GEO::MeshIOFlags flags;
  flags.set_element(GEO::MeshElementsFlags::MESH_ALL_ELEMENTS);
  flags.set_attributes(GEO::MeshAttributesFlags::MESH_ALL_ATTRIBUTES);
  if (!GEO::mesh_save(sf_mesh, geo_file, flags)) {
    std::cout << "Unable to save file at " << sf_path << std::endl;
    return false;
  }
  std::cout << "Done saving file  " << sf_path << std::endl;
  return true;
}