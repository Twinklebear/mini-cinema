#pragma once

#include <string>
#include <vector>
#include <ospcommon/math/vec.h>

using namespace ospcommon::math;

// Read the contents of a file into the string
std::string get_file_content(const std::string &fname);

std::string get_file_extension(const std::string &fname);

std::string get_file_basename(const std::string &path);

std::string get_file_basepath(const std::string &path);

bool starts_with(const std::string &str, const std::string &prefix);

std::vector<vec3f> generate_fibonacci_sphere(const float radius, const size_t n_points);

