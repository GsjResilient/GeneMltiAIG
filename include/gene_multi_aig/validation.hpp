#pragma once

#include "gene_multi_aig/aag.hpp"

#include <filesystem>
#include <string>

namespace gene_multi_aig {

struct ValidationResult {
    bool ok = true;
    std::string message;
};

ValidationResult validate_input_single_eq_pair(const std::filesystem::path& aig_path,
                                               Lit lhs,
                                               Lit rhs);

ValidationResult validate_output_single_eq_pair(const std::filesystem::path& aig_path,
                                                Lit expected_lhs,
                                                Lit expected_rhs);

} // namespace gene_multi_aig
