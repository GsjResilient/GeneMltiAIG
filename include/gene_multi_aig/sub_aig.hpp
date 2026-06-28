#pragma once

#include "gene_multi_aig/aag.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace gene_multi_aig {

std::vector<Lit> read_sub_variables(const Aag& base, const std::filesystem::path& path);
Aag constrain_output_with_units(const Aag& base, const std::vector<Lit>& unit_literals);
std::string assignment_bits(std::uint64_t mask, std::size_t width);
std::filesystem::path sub_aig_path(const std::filesystem::path& input,
                                   const std::filesystem::path& out_dir,
                                   const std::string& bits);

} // namespace gene_multi_aig
