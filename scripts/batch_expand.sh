#!/usr/bin/env bash
set -u

usage() {
    cat <<'EOF'
Usage: scripts/batch_expand.sh --input-dir DIR --list FILE --out-dir DIR --aigtoaig PATH --fastlec PATH [options]

Required:
  --input-dir DIR      Directory containing input AIG files
  --list FILE          Text file listing AIG names or paths relative to --input-dir
  --out-dir DIR        Directory for expanded AIG files
  --aigtoaig PATH      Path to aigtoaig
  --fastlec PATH       Path to fastLEC

Options:
  --ec-expand PATH     Path to ec_expand [default: build/ec_expand]
  --tmp-dir DIR        Temporary directory [default: OUT_DIR/tmp]
  --sub-out-dir DIR    Generate sub-AIGs under DIR/<input-stem>/
  --keep-temp          Keep ec_expand temporary files
  --verbose            Print ec_expand progress details
  --help               Show this help
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"

input_dir=""
list_file=""
out_dir=""
aigtoaig=""
fastlec=""
ec_expand="${repo_dir}/build/ec_expand"
tmp_dir=""
sub_out_dir=""
keep_temp=0
verbose=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            input_dir="$2"
            shift 2
            ;;
        --input-dir=*)
            input_dir="${1#*=}"
            shift
            ;;
        --list)
            list_file="$2"
            shift 2
            ;;
        --list=*)
            list_file="${1#*=}"
            shift
            ;;
        --out-dir)
            out_dir="$2"
            shift 2
            ;;
        --out-dir=*)
            out_dir="${1#*=}"
            shift
            ;;
        --aigtoaig)
            aigtoaig="$2"
            shift 2
            ;;
        --aigtoaig=*)
            aigtoaig="${1#*=}"
            shift
            ;;
        --fastlec)
            fastlec="$2"
            shift 2
            ;;
        --fastlec=*)
            fastlec="${1#*=}"
            shift
            ;;
        --ec-expand)
            ec_expand="$2"
            shift 2
            ;;
        --ec-expand=*)
            ec_expand="${1#*=}"
            shift
            ;;
        --tmp-dir)
            tmp_dir="$2"
            shift 2
            ;;
        --tmp-dir=*)
            tmp_dir="${1#*=}"
            shift
            ;;
        --sub-out-dir)
            sub_out_dir="$2"
            shift 2
            ;;
        --sub-out-dir=*)
            sub_out_dir="${1#*=}"
            shift
            ;;
        --keep-temp)
            keep_temp=1
            shift
            ;;
        --verbose)
            verbose=1
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "batch_expand: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

require_value() {
    local name="$1"
    local value="$2"
    if [[ -z "$value" ]]; then
        echo "batch_expand: missing required ${name}" >&2
        exit 2
    fi
}

require_value "--input-dir" "$input_dir"
require_value "--list" "$list_file"
require_value "--out-dir" "$out_dir"
require_value "--aigtoaig" "$aigtoaig"
require_value "--fastlec" "$fastlec"

if [[ ! -d "$input_dir" ]]; then
    echo "batch_expand: input directory does not exist: $input_dir" >&2
    exit 2
fi
if [[ ! -f "$list_file" ]]; then
    echo "batch_expand: list file does not exist: $list_file" >&2
    exit 2
fi
if [[ ! -x "$aigtoaig" ]]; then
    echo "batch_expand: aigtoaig is not executable: $aigtoaig" >&2
    exit 2
fi
if [[ ! -x "$fastlec" ]]; then
    echo "batch_expand: fastLEC is not executable: $fastlec" >&2
    exit 2
fi

if [[ -z "$tmp_dir" ]]; then
    tmp_dir="${out_dir}/tmp"
fi

mkdir -p "$out_dir" "$tmp_dir" "${out_dir}/logs"
if [[ -n "$sub_out_dir" ]]; then
    mkdir -p "$sub_out_dir"
fi

if [[ ! -x "$ec_expand" ]]; then
    if [[ "$ec_expand" != "${repo_dir}/build/ec_expand" ]]; then
        echo "batch_expand: ec_expand is not executable: $ec_expand" >&2
        exit 2
    fi
    cmake -S "$repo_dir" -B "${repo_dir}/build"
    cmake --build "${repo_dir}/build"
fi

success_count=0
failure_count=0
total_count=0
failures=()

while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
    line="${raw_line#"${raw_line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    if [[ -z "$line" || "${line:0:1}" == "#" ]]; then
        continue
    fi

    total_count=$((total_count + 1))
    input_path="${input_dir}/${line}"
    rel_dir="$(dirname "$line")"
    filename="$(basename "$line")"
    stem="${filename%.*}"

    if [[ "$rel_dir" == "." ]]; then
        output_dir="$out_dir"
        log_dir="${out_dir}/logs"
        sub_dir="${sub_out_dir}/${stem}"
    else
        output_dir="${out_dir}/${rel_dir}"
        log_dir="${out_dir}/logs/${rel_dir}"
        sub_dir="${sub_out_dir}/${rel_dir}/${stem}"
    fi

    output_path="${output_dir}/${stem}_expanded.aig"
    log_path="${log_dir}/${stem}.log"

    mkdir -p "$output_dir" "$log_dir"
    if [[ -n "$sub_out_dir" ]]; then
        mkdir -p "$sub_dir"
    fi

    if [[ ! -f "$input_path" ]]; then
        echo "[batch_expand] missing input: $input_path" | tee "$log_path"
        failure_count=$((failure_count + 1))
        failures+=("$line")
        continue
    fi

    cmd=(
        "$ec_expand"
        --in "$input_path"
        --out "$output_path"
        --tmp-dir "$tmp_dir"
        --aigtoaig "$aigtoaig"
        --fastlec "$fastlec"
    )
    if [[ -n "$sub_out_dir" ]]; then
        cmd+=(--sub-out-dir "$sub_dir")
    fi
    if [[ "$keep_temp" -eq 1 ]]; then
        cmd+=(--keep-temp)
    fi
    if [[ "$verbose" -eq 1 ]]; then
        cmd+=(--verbose)
    fi

    echo "[batch_expand] expanding $line -> $output_path"
    if "${cmd[@]}" >"$log_path" 2>&1; then
        success_count=$((success_count + 1))
    else
        failure_count=$((failure_count + 1))
        failures+=("$line")
        echo "[batch_expand] failed: $line (see $log_path)" >&2
    fi
done < "$list_file"

echo "[batch_expand] total: ${total_count}, succeeded: ${success_count}, failed: ${failure_count}"

if [[ "$failure_count" -gt 0 ]]; then
    echo "[batch_expand] failed files:" >&2
    for failure in "${failures[@]}"; do
        echo "  $failure" >&2
    done
    exit 1
fi

exit 0
