#pragma once

#include <string>
#include <vector>
#include <ospcommon/math/vec.h>
#include "json.hpp"

using namespace ospcommon::math;
using json = nlohmann::json;

// Read the contents of a file into the string
std::string get_file_content(const std::string &fname);

std::string get_file_extension(const std::string &fname);

std::string get_file_basename(const std::string &path);

std::string get_file_basepath(const std::string &path);

bool starts_with(const std::string &str, const std::string &prefix);

std::vector<vec3f> generate_fibonacci_sphere(const size_t n_points, const float radius);

template <typename T, size_t N>
inline vec_t<T, N> get_vec(const json &j)
{
    vec_t<T, N> v;
    for (size_t i = 0; i < N; ++i) {
        v[i] = j[i].get<T>();
    }
    return v;
}

