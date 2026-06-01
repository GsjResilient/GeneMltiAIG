#pragma once

#include "gene_multi_aig/aag.hpp"

#include <filesystem>
#include <string>

namespace gene_multi_aig {

struct ValidationResult {
    bool ok = true;
    std::string message;
};

struct FastLecValidationConfig {
    std::filesystem::path executable;
    std::filesystem::path log_path;
    unsigned cores = 8;
    unsigned timeout_seconds = 100;
    bool verbose = false;
};

ValidationResult validate_input_single_eq_pair(const std::filesystem::path& aig_path,
                                               Lit lhs,
                                               Lit rhs,
                                               const FastLecValidationConfig* fastlec);

ValidationResult validate_output_single_eq_pair(const std::filesystem::path& aig_path,
                                                Lit expected_lhs,
                                                Lit expected_rhs,
                                                const FastLecValidationConfig* fastlec);

} // namespace gene_multi_aig
