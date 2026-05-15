#!/usr/bin/env bash
set -euo pipefail

trim_spaces() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  echo "$s"
}

normalize_compact_age() {
  local age="$1"
  age="$(trim_spaces "$age")"
  age="${age,,}"
  if [[ -z "$age" ]]; then
    echo ""
  elif [[ "$age" =~ ^[0-9]+$ ]]; then
    echo "${age}y"
  else
    echo "$age"
  fi
}

is_valid_compact_age() {
  local age="$1"
  age="$(trim_spaces "$age")"
  [[ -z "$age" || "$age" =~ ^[0-9]+([dDmMyY])?$ ]]
}

reject_csv_field() {
  local label="$1"
  local value="$2"
  if [[ "$value" == *","* ]]; then
    echo "ERROR: ${label} contains a comma, which is not supported by this CSV generator: ${value}" >&2
    exit 1
  fi
}

csv_join_row() {
  local first="yes"
  local value
  for value in "$@"; do
    reject_csv_field "CSV field" "$value"
    if [[ "$first" == "yes" ]]; then
      printf '%s' "$value"
      first="no"
    else
      printf ',%s' "$value"
    fi
  done
  printf '\n'
}

normalize_prefix_from_suffix() {
  local filename="$1"
  local suffix="$2"
  if [[ "$filename" == *"$suffix" ]]; then
    echo "${filename%"$suffix"}"
  else
    echo "$filename"
  fi
}

normalize_for_match() {
  local s="$1"
  s="${s,,}"
  s="${s//[^a-z0-9]/}"
  echo "$s"
}

match_score() {
  local sample_prefix="$1"
  local candidate_prefix="$2"
  local s_norm c_norm
  s_norm="$(normalize_for_match "$sample_prefix")"
  c_norm="$(normalize_for_match "$candidate_prefix")"
  if [[ -z "$s_norm" || -z "$c_norm" ]]; then
    echo 0
  elif [[ "$sample_prefix" == "$candidate_prefix" ]]; then
    echo 1000
  elif [[ "$s_norm" == "$c_norm" ]]; then
    echo 900
  elif [[ "$s_norm" == *"$c_norm"* || "$c_norm" == *"$s_norm"* ]]; then
    local len_diff
    if (( ${#c_norm} > ${#s_norm} )); then
      len_diff=$(( ${#c_norm} - ${#s_norm} ))
    else
      len_diff=$(( ${#s_norm} - ${#c_norm} ))
    fi
    echo $((700 - len_diff))
  else
    local s_token c_token
    s_token="$(echo "$sample_prefix" | awk -F'[_ .-]' '{print tolower($1)}')"
    c_token="$(echo "$candidate_prefix" | awk -F'[_ .-]' '{print tolower($1)}')"
    [[ -n "$s_token" && "$s_token" == "$c_token" ]] && echo 400 || echo 0
  fi
}

find_best_match() {
  local sample_prefix="$1"
  shift
  local candidates=("$@")
  local best_score=0
  local best=()
  local candidate score
  for candidate in "${candidates[@]}"; do
    score="$(match_score "$sample_prefix" "$candidate")"
    if (( score > best_score )); then
      best_score=$score
      best=("$candidate")
    elif (( score > 0 && score == best_score )); then
      best+=("$candidate")
    fi
  done
  if (( best_score == 0 )); then
    echo ""
  elif (( ${#best[@]} > 1 )); then
    echo "AMBIGUOUS:$(IFS='|'; echo "${best[*]}")"
  else
    echo "${best[0]}"
  fi
}

load_prefix_map() {
  local dir="$1"
  local suffix="$2"
  local map_name="$3"
  local list_name="$4"
  declare -n out_map="$map_name"
  declare -n out_list="$list_name"
  local path name prefix
  while IFS= read -r path; do
    name="$(basename "$path")"
    prefix="$(normalize_prefix_from_suffix "$name" "$suffix")"
    reject_csv_field "sample prefix" "$prefix"
    reject_csv_field "path" "$path"
    if [[ -n "$prefix" ]]; then
      if [[ -n "${out_map[$prefix]:-}" ]]; then
        echo "WARNING: duplicate prefix '$prefix' found; keeping first:" >&2
        echo "         ${out_map[$prefix]}" >&2
        echo "         skipped: $path" >&2
      else
        out_map["$prefix"]="$path"
        out_list+=("$prefix")
      fi
    fi
  done < <(find "$dir" -maxdepth 1 -type f -name "*${suffix}" | sort)
}

prompt_age() {
  local sample="$1"
  local age_in
  while true; do
    printf "Age for sample '%s' [empty|e.g. 10d,9m,7y,7]: " "$sample" > /dev/tty
    read -r age_in < /dev/tty
    age_in="$(trim_spaces "$age_in")"
    if is_valid_compact_age "$age_in"; then
      normalize_compact_age "$age_in"
      return
    fi
    echo "Invalid age format. Use empty, <int>, or <int><d|m|y>." > /dev/tty
  done
}

select_age_prompting() {
  INCLUDE_AGE="no"
  PROMPT_AGE="no"
  read -r -p "Include optional age_of_onset column? [y/N]: " INCLUDE_AGE
  INCLUDE_AGE="$(trim_spaces "$INCLUDE_AGE")"
  INCLUDE_AGE="${INCLUDE_AGE,,}"
  if [[ "$INCLUDE_AGE" == "y" || "$INCLUDE_AGE" == "yes" ]]; then
    INCLUDE_AGE="yes"
    read -r -p "Prompt for age per sample? [y/N]: " PROMPT_AGE
    PROMPT_AGE="$(trim_spaces "$PROMPT_AGE")"
    PROMPT_AGE="${PROMPT_AGE,,}"
    [[ "$PROMPT_AGE" == "y" || "$PROMPT_AGE" == "yes" ]] && PROMPT_AGE="yes" || PROMPT_AGE="no"
    if [[ "$PROMPT_AGE" == "yes" && ! -r /dev/tty ]]; then
      echo "WARNING: no interactive TTY available; disabling per-sample age prompts." >&2
      PROMPT_AGE="no"
    fi
  fi
}

resolve_match_or_skip() {
  local sample="$1"
  local label="$2"
  local map_name="$3"
  local list_name="$4"
  declare -n ref_map="$map_name"
  declare -n ref_list="$list_name"
  local path="${ref_map[$sample]:-}"
  if [[ -n "$path" ]]; then
    echo "$path"
    return
  fi
  local matched_prefix
  matched_prefix="$(find_best_match "$sample" "${ref_list[@]}")"
  if [[ "$matched_prefix" == AMBIGUOUS:* || -z "$matched_prefix" ]]; then
    echo "WARNING: no unambiguous ${label} match for '$sample'; candidates/result: ${matched_prefix:-none}; skipping." >&2
    echo ""
    return
  fi
  echo "${ref_map[$matched_prefix]}"
}

discover_alignment_index() {
  local alignment_path="$1"
  local index_path=""
  if [[ "$alignment_path" == *.bam ]]; then
    [[ -f "${alignment_path}.bai" ]] && index_path="${alignment_path}.bai"
    [[ -z "$index_path" && -f "${alignment_path%.bam}.bai" ]] && index_path="${alignment_path%.bam}.bai"
  elif [[ "$alignment_path" == *.cram ]]; then
    [[ -f "${alignment_path}.crai" ]] && index_path="${alignment_path}.crai"
    [[ -z "$index_path" && -f "${alignment_path%.cram}.crai" ]] && index_path="${alignment_path%.cram}.crai"
  fi
  echo "$index_path"
}

generate_legacy_csv() {
  local data_dir data_suffix note_dir note_suffix out_csv format_choice
  echo
  read -r -p "Enter data-file directory (bam/cram/vcf/vcf.gz): " data_dir
  data_dir="$(trim_spaces "$data_dir")"
  [[ -d "$data_dir" ]] || { echo "ERROR: directory not found: $data_dir" >&2; exit 1; }
  cat <<'MSG'
Choose data file format:
  1) bam
  2) cram
  3) vcf
  4) vcf.gz
MSG
  read -r -p "Select format [1-4]: " format_choice
  format_choice="$(trim_spaces "$format_choice")"
  case "$format_choice" in
    1) data_suffix=".bam" ;;
    2) data_suffix=".cram" ;;
    3) data_suffix=".vcf" ;;
    4) data_suffix=".vcf.gz" ;;
    *) echo "ERROR: invalid selection: $format_choice" >&2; exit 1 ;;
  esac
  read -r -p "Enter medical-note text directory: " note_dir
  note_dir="$(trim_spaces "$note_dir")"
  [[ -d "$note_dir" ]] || { echo "ERROR: directory not found: $note_dir" >&2; exit 1; }
  read -r -p "Enter note filename suffix (example: _note.txt or .txt): " note_suffix
  note_suffix="$(trim_spaces "$note_suffix")"
  [[ -n "$note_suffix" ]] || { echo "ERROR: note suffix cannot be empty" >&2; exit 1; }
  read -r -p "Output CSV path [input_csv_generated.csv]: " out_csv
  out_csv="$(trim_spaces "$out_csv")"
  out_csv="${out_csv:-input_csv_generated.csv}"
  reject_csv_field "output path" "$out_csv"
  select_age_prompting

  declare -A note_map=()
  note_prefixes=()
  load_prefix_map "$note_dir" "$note_suffix" note_map note_prefixes
  [[ "${#note_map[@]}" -gt 0 ]] || { echo "ERROR: no note files matched suffix '$note_suffix' in $note_dir" >&2; exit 1; }

  local data_count=0 match_count=0 missing_count=0
  {
    if [[ "$INCLUDE_AGE" == "yes" ]]; then
      csv_join_row "sample" "file_path" "note_path" "age_of_onset"
    else
      csv_join_row "sample" "file_path" "note_path"
    fi
    while IFS= read -r data_path; do
      data_count=$((data_count + 1))
      local data_name sample note_path matched_prefix age_value
      data_name="$(basename "$data_path")"
      sample="$(normalize_prefix_from_suffix "$data_name" "$data_suffix")"
      reject_csv_field "sample" "$sample"
      reject_csv_field "path" "$data_path"
      note_path="${note_map[$sample]:-}"
      if [[ -z "$note_path" ]]; then
        matched_prefix="$(find_best_match "$sample" "${note_prefixes[@]}")"
        if [[ "$matched_prefix" == AMBIGUOUS:* ]]; then
          echo "WARNING: ambiguous note match for '$sample'; candidates: ${matched_prefix#AMBIGUOUS:}; skipping." >&2
          missing_count=$((missing_count + 1))
          continue
        elif [[ -n "$matched_prefix" ]]; then
          echo "INFO: rough matched '$sample' -> '$matched_prefix'" >&2
          note_path="${note_map[$matched_prefix]}"
        else
          echo "WARNING: no note match for sample prefix '$sample' (${data_name})" >&2
          missing_count=$((missing_count + 1))
          continue
        fi
      fi
      match_count=$((match_count + 1))
      age_value=""
      [[ "$PROMPT_AGE" == "yes" ]] && age_value="$(prompt_age "$sample")"
      if [[ "$INCLUDE_AGE" == "yes" ]]; then
        csv_join_row "$sample" "$data_path" "$note_path" "$age_value"
      else
        csv_join_row "$sample" "$data_path" "$note_path"
      fi
    done < <(find "$data_dir" -maxdepth 1 -type f -name "*${data_suffix}" | sort)
  } > "$out_csv"
  echo
  echo "Done."
  echo "  Data files scanned : $data_count"
  echo "  Matched rows       : $match_count"
  echo "  Missing note match : $missing_count"
  echo "  Output CSV         : $out_csv"
  [[ "$match_count" -gt 0 ]] || { echo "ERROR: zero matched samples. Check prefixes and note suffix." >&2; exit 1; }
}

choose_input_kind() {
  cat <<'MSG'
Choose input kind:
  1) annotated_snv
  2) vcf_snv
  3) vcf_sv
  4) bam_ngs
  5) cram_ngs
MSG
  read -r -p "Select input kind [1-5]: " input_choice
  input_choice="$(trim_spaces "$input_choice")"
  case "$input_choice" in
    1) INPUT_KIND="annotated_snv" ;;
    2) INPUT_KIND="vcf_snv" ;;
    3) INPUT_KIND="vcf_sv" ;;
    4) INPUT_KIND="bam_ngs" ;;
    5) INPUT_KIND="cram_ngs" ;;
    *) echo "ERROR: invalid selection: $input_choice" >&2; exit 1 ;;
  esac
}

choose_phenotype() {
  cat <<'MSG'
Choose phenotype format:
  1) clinical_note
  2) hpo
MSG
  read -r -p "Select phenotype format [1-2]: " phenotype_choice
  phenotype_choice="$(trim_spaces "$phenotype_choice")"
  case "$phenotype_choice" in
    1) PHENOTYPE_FORMAT="clinical_note" ;;
    2) PHENOTYPE_FORMAT="hpo" ;;
    *) echo "ERROR: invalid selection: $phenotype_choice" >&2; exit 1 ;;
  esac
}

generate_unified_csv() {
  echo
  choose_input_kind
  choose_phenotype
  local phenotype_dir phenotype_suffix out_csv
  read -r -p "Enter phenotype directory: " phenotype_dir
  phenotype_dir="$(trim_spaces "$phenotype_dir")"
  [[ -d "$phenotype_dir" ]] || { echo "ERROR: directory not found: $phenotype_dir" >&2; exit 1; }
  read -r -p "Enter phenotype filename suffix (example: _note.txt or .hpo.txt): " phenotype_suffix
  phenotype_suffix="$(trim_spaces "$phenotype_suffix")"
  [[ -n "$phenotype_suffix" ]] || { echo "ERROR: phenotype suffix cannot be empty" >&2; exit 1; }
  read -r -p "Output CSV path [input_csv_generated.csv]: " out_csv
  out_csv="$(trim_spaces "$out_csv")"
  out_csv="${out_csv:-input_csv_generated.csv}"
  reject_csv_field "output path" "$out_csv"
  select_age_prompting

  declare -A phenotype_map=()
  phenotype_prefixes=()
  load_prefix_map "$phenotype_dir" "$phenotype_suffix" phenotype_map phenotype_prefixes
  [[ "${#phenotype_map[@]}" -gt 0 ]] || { echo "ERROR: no phenotype files matched suffix '$phenotype_suffix' in $phenotype_dir" >&2; exit 1; }

  declare -A primary_map=()
  primary_prefixes=()
  declare -A secondary_map=()
  secondary_prefixes=()
  declare -A alignment_map=()
  alignment_prefixes=()
  declare -A sv_vcf_map=()
  sv_vcf_prefixes=()
  local primary_dir primary_suffix secondary_dir secondary_suffix include_alignment="no"

  if [[ "$INPUT_KIND" == "annotated_snv" ]]; then
    read -r -p "Enter ANNOVAR TXT directory: " primary_dir
    primary_dir="$(trim_spaces "$primary_dir")"
    [[ -d "$primary_dir" ]] || { echo "ERROR: directory not found: $primary_dir" >&2; exit 1; }
    read -r -p "Enter ANNOVAR TXT suffix [.hg38_multianno.txt]: " primary_suffix
    primary_suffix="$(trim_spaces "$primary_suffix")"
    primary_suffix="${primary_suffix:-.hg38_multianno.txt}"
    read -r -p "Enter ANNOVAR VCF directory: " secondary_dir
    secondary_dir="$(trim_spaces "$secondary_dir")"
    [[ -d "$secondary_dir" ]] || { echo "ERROR: directory not found: $secondary_dir" >&2; exit 1; }
    read -r -p "Enter ANNOVAR VCF suffix [.hg38_multianno.vcf]: " secondary_suffix
    secondary_suffix="$(trim_spaces "$secondary_suffix")"
    secondary_suffix="${secondary_suffix:-.hg38_multianno.vcf}"
    read -r -p "Also include BAM/CRAM paths for SV/STR/mito analysis? [y/N]: " include_alignment
    include_alignment="$(trim_spaces "$include_alignment")"
    include_alignment="${include_alignment,,}"
    [[ "$include_alignment" == "y" || "$include_alignment" == "yes" ]] && include_alignment="yes" || include_alignment="no"
  else
    read -r -p "Enter input-data directory: " primary_dir
    primary_dir="$(trim_spaces "$primary_dir")"
    [[ -d "$primary_dir" ]] || { echo "ERROR: directory not found: $primary_dir" >&2; exit 1; }
    case "$INPUT_KIND" in
      vcf_snv|vcf_sv)
        read -r -p "Enter VCF suffix [.vcf]: " primary_suffix
        primary_suffix="$(trim_spaces "$primary_suffix")"
        primary_suffix="${primary_suffix:-.vcf}"
        ;;
      bam_ngs) primary_suffix=".bam" ;;
      cram_ngs) primary_suffix=".cram" ;;
    esac
  fi
  load_prefix_map "$primary_dir" "$primary_suffix" primary_map primary_prefixes
  [[ "${#primary_map[@]}" -gt 0 ]] || { echo "ERROR: no input files matched suffix '$primary_suffix' in $primary_dir" >&2; exit 1; }
  if [[ "$INPUT_KIND" == "annotated_snv" ]]; then
    load_prefix_map "$secondary_dir" "$secondary_suffix" secondary_map secondary_prefixes
    [[ "${#secondary_map[@]}" -gt 0 ]] || { echo "ERROR: no ANNOVAR VCF files matched suffix '$secondary_suffix' in $secondary_dir" >&2; exit 1; }
  fi
  if [[ "$include_alignment" == "yes" ]]; then
    local alignment_dir alignment_suffix sv_vcf_dir sv_vcf_suffix alignment_choice
    read -r -p "Enter alignment directory: " alignment_dir
    alignment_dir="$(trim_spaces "$alignment_dir")"
    [[ -d "$alignment_dir" ]] || { echo "ERROR: directory not found: $alignment_dir" >&2; exit 1; }
    cat <<'MSG'
Choose alignment format:
  1) bam
  2) cram
MSG
    read -r -p "Select alignment format [1-2]: " alignment_choice
    alignment_choice="$(trim_spaces "$alignment_choice")"
    case "$alignment_choice" in
      1) alignment_suffix=".bam" ;;
      2) alignment_suffix=".cram" ;;
      *) echo "ERROR: invalid selection: $alignment_choice" >&2; exit 1 ;;
    esac
    load_prefix_map "$alignment_dir" "$alignment_suffix" alignment_map alignment_prefixes
    [[ "${#alignment_map[@]}" -gt 0 ]] || { echo "ERROR: no alignment files matched suffix '$alignment_suffix' in $alignment_dir" >&2; exit 1; }
    read -r -p "Enter annotated SV VCF directory: " sv_vcf_dir
    sv_vcf_dir="$(trim_spaces "$sv_vcf_dir")"
    [[ -d "$sv_vcf_dir" ]] || { echo "ERROR: directory not found: $sv_vcf_dir" >&2; exit 1; }
    read -r -p "Enter annotated SV VCF suffix [.hg38_multianno.vcf]: " sv_vcf_suffix
    sv_vcf_suffix="$(trim_spaces "$sv_vcf_suffix")"
    sv_vcf_suffix="${sv_vcf_suffix:-.hg38_multianno.vcf}"
    load_prefix_map "$sv_vcf_dir" "$sv_vcf_suffix" sv_vcf_map sv_vcf_prefixes
    [[ "${#sv_vcf_map[@]}" -gt 0 ]] || { echo "ERROR: no annotated SV VCF files matched suffix '$sv_vcf_suffix' in $sv_vcf_dir" >&2; exit 1; }
  fi

  local data_count=0 match_count=0 missing_count=0
  {
    csv_join_row sample input_kind phenotype_path phenotype_format age_of_onset snv_txt_path snv_vcf_path sv_vcf_path vcf_path alignment_path alignment_index_path
    local sample primary_path phenotype_path snv_txt_path snv_vcf_path sv_vcf_path vcf_path alignment_path alignment_index_path age_value
    for sample in "${primary_prefixes[@]}"; do
      data_count=$((data_count + 1))
      primary_path="${primary_map[$sample]}"
      phenotype_path="$(resolve_match_or_skip "$sample" "phenotype" phenotype_map phenotype_prefixes)"
      [[ -n "$phenotype_path" ]] || { missing_count=$((missing_count + 1)); continue; }
      snv_txt_path=""
      snv_vcf_path=""
      sv_vcf_path=""
      vcf_path=""
      alignment_path=""
      alignment_index_path=""
      if [[ "$INPUT_KIND" == "annotated_snv" ]]; then
        snv_txt_path="$primary_path"
        snv_vcf_path="$(resolve_match_or_skip "$sample" "ANNOVAR VCF" secondary_map secondary_prefixes)"
        [[ -n "$snv_vcf_path" ]] || { missing_count=$((missing_count + 1)); continue; }
        if [[ "$include_alignment" == "yes" ]]; then
          alignment_path="$(resolve_match_or_skip "$sample" "alignment" alignment_map alignment_prefixes)"
          [[ -n "$alignment_path" ]] || { missing_count=$((missing_count + 1)); continue; }
          sv_vcf_path="$(resolve_match_or_skip "$sample" "annotated SV VCF" sv_vcf_map sv_vcf_prefixes)"
          [[ -n "$sv_vcf_path" ]] || { missing_count=$((missing_count + 1)); continue; }
          alignment_index_path="$(discover_alignment_index "$alignment_path")"
        fi
      elif [[ "$INPUT_KIND" == "vcf_snv" || "$INPUT_KIND" == "vcf_sv" ]]; then
        vcf_path="$primary_path"
      else
        alignment_path="$primary_path"
        alignment_index_path="$(discover_alignment_index "$alignment_path")"
      fi
      age_value=""
      [[ "$PROMPT_AGE" == "yes" ]] && age_value="$(prompt_age "$sample")"
      csv_join_row "$sample" "$INPUT_KIND" "$phenotype_path" "$PHENOTYPE_FORMAT" "$age_value" "$snv_txt_path" "$snv_vcf_path" "$sv_vcf_path" "$vcf_path" "$alignment_path" "$alignment_index_path"
      match_count=$((match_count + 1))
    done
  } > "$out_csv"
  echo
  echo "Done."
  echo "  Input files scanned : $data_count"
  echo "  Matched rows        : $match_count"
  echo "  Missing pair/match  : $missing_count"
  echo "  Output CSV          : $out_csv"
  [[ "$match_count" -gt 0 ]] || { echo "ERROR: zero matched samples. Check suffixes and prefixes." >&2; exit 1; }
}

update_existing_annotated_snv_csv() {
  local input_csv out_csv sv_vcf_dir sv_vcf_suffix overwrite_existing
  echo
  read -r -p "Existing annotated_snv unified CSV path: " input_csv
  input_csv="$(trim_spaces "$input_csv")"
  [[ -f "$input_csv" ]] || { echo "ERROR: CSV file not found: $input_csv" >&2; exit 1; }
  read -r -p "Output updated CSV path [updated_annotated_snv.csv]: " out_csv
  out_csv="$(trim_spaces "$out_csv")"
  out_csv="${out_csv:-updated_annotated_snv.csv}"
  reject_csv_field "output path" "$out_csv"
  read -r -p "Enter annotated SV VCF directory: " sv_vcf_dir
  sv_vcf_dir="$(trim_spaces "$sv_vcf_dir")"
  [[ -d "$sv_vcf_dir" ]] || { echo "ERROR: directory not found: $sv_vcf_dir" >&2; exit 1; }
  read -r -p "Enter annotated SV VCF suffix [.hg38_multianno.vcf]: " sv_vcf_suffix
  sv_vcf_suffix="$(trim_spaces "$sv_vcf_suffix")"
  sv_vcf_suffix="${sv_vcf_suffix:-.hg38_multianno.vcf}"
  read -r -p "Overwrite existing non-empty sv_vcf_path values? [y/N]: " overwrite_existing
  overwrite_existing="$(trim_spaces "$overwrite_existing")"
  overwrite_existing="${overwrite_existing,,}"
  [[ "$overwrite_existing" == "y" || "$overwrite_existing" == "yes" ]] && overwrite_existing="yes" || overwrite_existing="no"
  declare -A sv_update_map=()
  sv_update_prefixes=()
  load_prefix_map "$sv_vcf_dir" "$sv_vcf_suffix" sv_update_map sv_update_prefixes
  [[ "${#sv_update_map[@]}" -gt 0 ]] || { echo "ERROR: no annotated SV VCF files matched suffix '$sv_vcf_suffix' in $sv_vcf_dir" >&2; exit 1; }

  local header
  IFS= read -r header < "$input_csv" || { echo "ERROR: CSV is empty: $input_csv" >&2; exit 1; }
  local headers=()
  IFS=',' read -r -a headers <<< "$header"
  local sample_idx=-1 input_kind_idx=-1 snv_vcf_idx=-1 sv_idx=-1 alignment_idx=-1 added_sv_column="no"
  local idx column
  for idx in "${!headers[@]}"; do
    column="$(trim_spaces "${headers[$idx]}")"
    headers[$idx]="$column"
    case "$column" in
      sample) sample_idx=$idx ;;
      input_kind) input_kind_idx=$idx ;;
      snv_vcf_path) snv_vcf_idx=$idx ;;
      sv_vcf_path) sv_idx=$idx ;;
      alignment_path) alignment_idx=$idx ;;
    esac
  done
  (( sample_idx >= 0 && input_kind_idx >= 0 && alignment_idx >= 0 )) || { echo "ERROR: CSV must include sample, input_kind, and alignment_path columns." >&2; exit 1; }
  if (( sv_idx < 0 )); then
    if (( snv_vcf_idx >= 0 )); then
      sv_idx=$((snv_vcf_idx + 1))
      headers=("${headers[@]:0:sv_idx}" "sv_vcf_path" "${headers[@]:sv_idx}")
    else
      headers+=("sv_vcf_path")
      sv_idx=$((${#headers[@]} - 1))
    fi
    added_sv_column="yes"
  fi
  local row_count=0 updated_count=0 skipped_count=0 missing_count=0
  {
    csv_join_row "${headers[@]}"
    local line
    while IFS= read -r line || [[ -n "$line" ]]; do
      [[ -n "$(trim_spaces "$line")" ]] || continue
      local values=()
      IFS=',' read -r -a values <<< "$line"
      if [[ "$added_sv_column" == "yes" ]]; then
        values=("${values[@]:0:sv_idx}" "" "${values[@]:sv_idx}")
      fi
      while (( ${#values[@]} < ${#headers[@]} )); do
        values+=("")
      done
      row_count=$((row_count + 1))
      local sample input_kind current_sv alignment_path matched_prefix matched_sv
      sample="$(trim_spaces "${values[$sample_idx]}")"
      input_kind="$(trim_spaces "${values[$input_kind_idx]}")"
      current_sv="$(trim_spaces "${values[$sv_idx]}")"
      alignment_path="$(trim_spaces "${values[$alignment_idx]}")"
      if [[ "$input_kind" != "annotated_snv" || -z "$alignment_path" ]]; then
        skipped_count=$((skipped_count + 1))
        csv_join_row "${values[@]}"
        continue
      fi
      if [[ -n "$current_sv" && "$overwrite_existing" != "yes" ]]; then
        skipped_count=$((skipped_count + 1))
        csv_join_row "${values[@]}"
        continue
      fi
      matched_sv="${sv_update_map[$sample]:-}"
      if [[ -z "$matched_sv" ]]; then
        matched_prefix="$(find_best_match "$sample" "${sv_update_prefixes[@]}")"
        if [[ "$matched_prefix" == AMBIGUOUS:* || -z "$matched_prefix" ]]; then
          echo "WARNING: no unambiguous annotated SV VCF match for '$sample'; candidates/result: ${matched_prefix:-none}; leaving sv_vcf_path unchanged." >&2
          missing_count=$((missing_count + 1))
          csv_join_row "${values[@]}"
          continue
        fi
        matched_sv="${sv_update_map[$matched_prefix]}"
      fi
      values[$sv_idx]="$matched_sv"
      updated_count=$((updated_count + 1))
      csv_join_row "${values[@]}"
    done
  } < <(tail -n +2 "$input_csv") > "$out_csv"
  echo
  echo "Done."
  echo "  Rows scanned     : $row_count"
  echo "  SV paths updated : $updated_count"
  echo "  Rows skipped     : $skipped_count"
  echo "  Missing matches  : $missing_count"
  echo "  Output CSV       : $out_csv"
}

main() {
  echo "PipeVar_mito input CSV generator"
  echo
  cat <<'MSG'
Choose action:
  1) Generate legacy CSV
  2) Generate unified CSV
  3) Update annotated-SNV unified CSV with annotated SV VCF paths
MSG
  read -r -p "Select action [1-3]: " action_choice
  action_choice="$(trim_spaces "$action_choice")"
  case "$action_choice" in
    1|"") generate_legacy_csv ;;
    2) generate_unified_csv ;;
    3) update_existing_annotated_snv_csv ;;
    *) echo "ERROR: invalid selection: $action_choice" >&2; exit 1 ;;
  esac
}

main "$@"
