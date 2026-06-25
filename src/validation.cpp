#include "gene_multi_aig/validation.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

namespace gene_multi_aig {

namespace {

enum class FastLecRunMode {
    StopAtClassCount,
    RunToCompletion
};

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

PotentialEquivClassCount count_potential_equiv_classes(
    const std::filesystem::path& aig_path,
    const FastLecValidationConfig* fastlec,
    FastLecRunMode mode) {
    if (fastlec == nullptr || fastlec->executable.empty()) {
        return PotentialEquivClassCount{false, 0, "fastLEC validation is required"};
    }

    if (fastlec->log_path.empty()) {
        return PotentialEquivClassCount{false, 0, "fastLEC log path is empty"};
    }

    if (!std::filesystem::exists(fastlec->executable)) {
        return PotentialEquivClassCount{false, 0,
                                        "fastLEC executable does not exist: " +
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

    std::string command = "rm -f " + log_path + "; ";
    if (mode == FastLecRunMode::StopAtClassCount) {
        command +=
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
    } else {
        command += fastlec_invocation + " > " + log_path + " 2>&1";
    }
    if (fastlec->verbose) {
        std::cerr << "[validation] " << command << '\n';
    }

    std::system(command.c_str());

    std::ifstream log(fastlec->log_path);
    if (!log) {
        return PotentialEquivClassCount{false, 0,
                                        "could not read fastLEC log: " +
                                            fastlec->log_path.string()};
    }

    bool saw_class_count = false;
    bool saw_bug = false;
    bool not_equivalent = false;
    bool equivalent = false;
    bool unknown = false;
    unsigned remaining_classes = 0;
    const std::regex remain_re(R"(remain\s+([0-9]+)\s+potent(?:ial|ail)-eql classes)");

    std::string line;
    while (std::getline(log, line)) {
        std::smatch match;
        if (std::regex_search(line, match, remain_re)) {
            saw_class_count = true;
            remaining_classes = static_cast<unsigned>(std::stoul(match[1].str()));
        }
        if (line.find("Find bugs") != std::string::npos) {
            saw_bug = true;
        }
        if (line.find("s Not Equivalent.") != std::string::npos) {
            not_equivalent = true;
        }
        if (line.find("s Equivalent.") != std::string::npos) {
            equivalent = true;
        }
        if (line.find("s Unknown.") != std::string::npos) {
            unknown = true;
        }
    }

    if (!saw_class_count) {
        std::string detail = "fastLEC log did not contain potential-equivalence class count";
        if (saw_bug || not_equivalent) {
            detail += " and reported a satisfiable counterexample";
        } else if (unknown) {
            detail += " and ended with Unknown";
        } else if (equivalent) {
            detail += " and ended with Equivalent";
        }
        return PotentialEquivClassCount{
            false,
            0,
            detail + "; see " + fastlec->log_path.string(),
            saw_bug,
            not_equivalent,
            equivalent,
            unknown};
    }

    return PotentialEquivClassCount{
        true,
        remaining_classes,
        "fastLEC found " + std::to_string(remaining_classes) +
            " potential-equivalence classes; see " + fastlec->log_path.string(),
        saw_bug,
        not_equivalent,
        equivalent,
        unknown};
}

ValidationResult validate_unique_equiv_pair(const std::filesystem::path& aig_path,
                                            Lit lhs,
                                            Lit rhs,
                                            const FastLecValidationConfig* fastlec) {
    const PotentialEquivClassCount count =
        count_potential_equiv_classes(aig_path, fastlec, FastLecRunMode::StopAtClassCount);
    if (!count.ok) {
        return ValidationResult{false, count.message};
    }

    const unsigned remaining_classes = count.count;
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

PotentialEquivClassCount count_output_potential_equiv_classes(
    const std::filesystem::path& aig_path,
    const FastLecValidationConfig* fastlec,
    bool stop_at_class_count) {
    return count_potential_equiv_classes(
        aig_path,
        fastlec,
        stop_at_class_count ? FastLecRunMode::StopAtClassCount
                            : FastLecRunMode::RunToCompletion);
}

} // namespace gene_multi_aig
