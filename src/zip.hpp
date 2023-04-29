#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Simplified version of partio
namespace ZIP {
// Create streams thath write .gz
std::ostream* Gzip_Out(const std::string& filename, std::ios::openmode mode);
}  // namespace ZIP