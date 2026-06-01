#include "gene_multi_aig/validation.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

namespace gene_multi_aig {

namespace {

std::string shell_quote(const std::filesystem::path& path) {
#ifdef _WIN32
    std::string value = path.string();
    std::string out = "\"";
    for (const char ch : value) {
        if (ch == '"') {
            out += "\\\"";
        } else {
            out += ch;
        }
    }
    out += "\"";
    return out;
#else
    std::string value = path.string();
    std::string out = "'";
    for (const char ch : value) {
        if (ch == '\'') {
            out += "'\\''";
        } else {
            out += ch;
        }
    }
    out += "'";
    return out;
#endif
}

ValidationResult validate_unique_equiv_pair(const std::filesystem::path& aig_path,
                                            Lit lhs,
                                            Lit rhs,
                                            const FastLecValidationConfig* fastlec) {
    if (fastlec == nullptr || fastlec->executable.empty()) {
        return ValidationResult{false, "fastLEC validation is required"};
    }

    if (fastlec->log_path.empty()) {
        return ValidationResult{false, "fastLEC log path is empty"};
    }

    if (!std::filesystem::exists(fastlec->executable)) {
        return ValidationResult{false, "fastLEC executable does not exist: " +
                                           fastlec->executable.string()};
    }

    const std::string log_path = shell_quote(fastlec->log_path);
    const std::string fastlec_invocation =
        shell_quote(fastlec->executable) +
        " -i " + shell_quote(aig_path) +
        " -m SAT_sweeping" +
        " -c " + std::to_string(fastlec->cores) +
        " -v 2" +
        " -t " + std::to_string(fastlec->timeout_seconds);
    const std::string command =
        "rm -f " + log_path + "; " +
        "( " + fastlec_invocation + " > " + log_path + " 2>&1 & "
        "pid=$!; "
        "while kill -0 \"$pid\" 2>/dev/null; do "
        "if grep -Eq 'remain[[:space:]]+[0-9]+[[:space:]]+potent(ial|ail)-eql classes' " +
        log_path + " 2>/dev/null; then "
        "kill \"$pid\" 2>/dev/null; "
        "wait \"$pid\" 2>/dev/null; "
        "exit 0; "
        "fi; "
        "sleep 0.05; "
        "done; "
        "wait \"$pid\" 2>/dev/null; "
        "exit 0; "
        ")";
    if (fastlec->verbose) {
        std::cerr << "[validation] " << command << '\n';
    }

    std::system(command.c_str());

    std::ifstream log(fastlec->log_path);
    if (!log) {
        return ValidationResult{false, "could not read fastLEC log: " +
                                           fastlec->log_path.string()};
    }

    bool saw_class_count = false;
    unsigned remaining_classes = 0;
    const std::regex remain_re(R"(remain\s+([0-9]+)\s+potent(?:ial|ail)-eql classes)");

    std::string line;
    while (std::getline(log, line)) {
        std::smatch match;
        if (std::regex_search(line, match, remain_re)) {
            saw_class_count = true;
            remaining_classes = static_cast<unsigned>(std::stoul(match[1].str()));
        }
    }

    if (!saw_class_count) {
        return ValidationResult{false,
                                "fastLEC log did not contain potential-equivalence class count; see " +
                                    fastlec->log_path.string()};
    }
    if (remaining_classes != 1U) {
        return ValidationResult{false,
                                "expected exactly 1 potential-equivalence class, found " +
                                    std::to_string(remaining_classes) +
                                    "; see " + fastlec->log_path.string()};
    }

    return ValidationResult{true,
                            "fastLEC found exactly one potential-equivalence class for roots " +
                                std::to_string(lhs) + " and " + std::to_string(rhs)};
}

} // namespace

ValidationResult validate_input_single_eq_pair(const std::filesystem::path& aig_path,
                                               Lit lhs,
                                               Lit rhs,
                                               const FastLecValidationConfig* fastlec) {
    return validate_unique_equiv_pair(aig_path, lhs, rhs, fastlec);
}

ValidationResult validate_output_single_eq_pair(const std::filesystem::path& aig_path,
                                                Lit expected_lhs,
                                                Lit expected_rhs,
                                                const FastLecValidationConfig* fastlec) {
    return validate_unique_equiv_pair(aig_path, expected_lhs, expected_rhs, fastlec);
}

} // namespace gene_multi_aig
