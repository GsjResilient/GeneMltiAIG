#include "gene_multi_aig/validation.hpp"

namespace gene_multi_aig {

ValidationResult validate_input_single_eq_pair(const std::filesystem::path&,
                                               Lit,
                                               Lit) {
    return ValidationResult{true, "input validation placeholder"};
}

ValidationResult validate_output_single_eq_pair(const std::filesystem::path&,
                                                Lit,
                                                Lit) {
    return ValidationResult{true, "output validation placeholder"};
}

} // namespace gene_multi_aig
