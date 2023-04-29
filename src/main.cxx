#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "io.h"
#include "main_gui_cxx.h"
#include "params.h"

int main(int argc, char** argv) {
  if (argc < 1) {
    std::cerr << "Usage: " << argv[0]
              << " <tet_mesh.tet/vtk> (e.g.: " << argv[0]
              << " ../data/joint.tet)" << std::endl;
    return 1;
  }

  std::string surface_path = argv[1];
  // read surface file from tetwild/ftetwild
  // this will make sure all geogram algoriths can be used
  GEO::initialize();
  GEO::CmdLine::import_arg_group("algo");

  GEO::Mesh sf_mesh;
  if (!load_surface_mesh(surface_path, sf_mesh)) {
    std::cerr << surface_path << ": could not load surface file" << std::endl;
    return 1;
  }
  std::cout << "Done loading surface mesh" << std::endl;

  std::vector<MedialSphere> all_medial_spheres;

  // show Gui
  MainGuiWindow main_gui;
  main_gui.set_sf_mesh(sf_mesh);
  main_gui.set_all_medial_spheres(all_medial_spheres);
  main_gui.show();

  return 0;
}
