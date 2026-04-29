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
    return
  fi
  if [[ "$age" =~ ^[0-9]+$ ]]; then
    echo "${age}y"
    return
  fi
  echo "$age"
}

is_valid_compact_age() {
  local age="$1"
  age="$(trim_spaces "$age")"
  if [[ -z "$age" ]]; then
    return 0
  fi
  [[ "$age" =~ ^[0-9]+([dDmMyY])?$ ]]
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
    return
  fi
  if [[ "$sample_prefix" == "$candidate_prefix" ]]; then
    echo 1000
    return
  fi
  if [[ "$s_norm" == "$c_norm" ]]; then
    echo 900
    return
  fi
  if [[ "$s_norm" == *"$c_norm"* || "$c_norm" == *"$s_norm"* ]]; then
    local len_diff=${#s_norm}
    if (( ${#c_norm} > len_diff )); then
      len_diff=$(( ${#c_norm} - ${#s_norm} ))
    else
      len_diff=$(( ${#s_norm} - ${#c_norm} ))
    fi
    echo $((700 - len_diff))
    return
  fi
  echo 0
}

find_best_match() {
  local sample_prefix="$1"
  shift
  local candidates=("$@")
  local best_prefix=""
  local best_score=0
  local best_tie=0
  local candidate score

  for candidate in "${candidates[@]}"; do
    score="$(match_score "$sample_prefix" "$candidate")"
    if (( score > best_score )); then
      best_score=$score
      best_prefix="$candidate"
      best_tie=0
    elif (( score > 0 && score == best_score )); then
      best_tie=1
    fi
  done

  if (( best_score == 0 )); then
    echo ""
    return
  fi
  if (( best_tie == 1 )); then
    echo "AMBIGUOUS:${best_prefix}"
    return
  fi
  echo "$best_prefix"
}

collect_by_suffix() {
  local dir="$1"
  local suffix="$2"
  find "$dir" -maxdepth 1 -type f -name "*${suffix}" | sort
}

load_prefix_map() {
  local dir="$1"
  local suffix="$2"
  local map_name="$3"
  local list_name="$4"
  declare -n out_map="$map_name"
  declare -n out_list="$list_name"

  while IFS= read -r path; do
    local name prefix
    name="$(basename "$path")"
    prefix="$(normalize_prefix_from_suffix "$name" "$suffix")"
    if [[ -n "$prefix" && -z "${out_map[$prefix]:-}" ]]; then
      out_map["$prefix"]="$path"
      out_list+=("$prefix")
    fi
  done < <(collect_by_suffix "$dir" "$suffix")
}

csv_split_line() {
  local line="$1"
  local array_name="$2"
  declare -n out_array="$array_name"
  IFS=',' read -r -a out_array <<< "$line"
}

csv_join_row() {
  local first="yes"
  local value
  for value in "$@"; do
    if [[ "$first" == "yes" ]]; then
      printf '%s' "$value"
      first="no"
    else
      printf ',%s' "$value"
    fi
  done
  printf '\n'
}

array_insert_value() {
  local array_name="$1"
  local insert_idx="$2"
  local insert_value="$3"
  declare -n target_array="$array_name"
  local rebuilt=()
  local idx
  for idx in "${!target_array[@]}"; do
    if (( idx == insert_idx )); then
      rebuilt+=("$insert_value")
    fi
    rebuilt+=("${target_array[$idx]}")
  done
  if (( insert_idx >= ${#target_array[@]} )); then
    rebuilt+=("$insert_value")
  fi
  target_array=("${rebuilt[@]}")
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

update_existing_annotated_snv_csv() {
  local input_csv out_csv sv_vcf_dir sv_vcf_suffix overwrite_existing
  echo
  read -r -p "Existing annotated_snv CSV path: " input_csv
  input_csv="$(trim_spaces "$input_csv")"
  [[ -f "$input_csv" ]] || { echo "ERROR: CSV file not found: $input_csv" >&2; exit 1; }

  read -r -p "Output updated CSV path [updated_annotated_snv.csv]: " out_csv
  out_csv="$(trim_spaces "$out_csv")"
  out_csv="${out_csv:-updated_annotated_snv.csv}"

  read -r -p "Enter annotated SV VCF directory: " sv_vcf_dir
  sv_vcf_dir="$(trim_spaces "$sv_vcf_dir")"
  [[ -d "$sv_vcf_dir" ]] || { echo "ERROR: directory not found: $sv_vcf_dir" >&2; exit 1; }

  read -r -p "Enter annotated SV VCF suffix [.hg38_multianno.vcf]: " sv_vcf_suffix
  sv_vcf_suffix="$(trim_spaces "$sv_vcf_suffix")"
  sv_vcf_suffix="${sv_vcf_suffix:-.hg38_multianno.vcf}"

  read -r -p "Overwrite existing non-empty sv_vcf_path values? [y/N]: " overwrite_existing
  overwrite_existing="$(trim_spaces "$overwrite_existing")"
  overwrite_existing="${overwrite_existing,,}"
  if [[ "$overwrite_existing" == "y" || "$overwrite_existing" == "yes" ]]; then
    overwrite_existing="yes"
  else
    overwrite_existing="no"
  fi

  declare -A update_sv_vcf_map=()
  update_sv_vcf_prefixes=()
  load_prefix_map "$sv_vcf_dir" "$sv_vcf_suffix" update_sv_vcf_map update_sv_vcf_prefixes
  if [[ "${#update_sv_vcf_map[@]}" -eq 0 ]]; then
    echo "ERROR: no annotated SV VCF files matched suffix '$sv_vcf_suffix' in $sv_vcf_dir" >&2
    exit 1
  fi

  local header
  IFS= read -r header < "$input_csv" || { echo "ERROR: CSV is empty: $input_csv" >&2; exit 1; }
  local headers=()
  csv_split_line "$header" headers

  local sample_idx=-1 input_kind_idx=-1 snv_vcf_idx=-1 sv_idx=-1 added_sv_column="no"
  local idx column
  for idx in "${!headers[@]}"; do
    column="$(trim_spaces "${headers[$idx]}")"
    headers[$idx]="$column"
    case "$column" in
      sample) sample_idx=$idx ;;
      input_kind) input_kind_idx=$idx ;;
      snv_vcf_path) snv_vcf_idx=$idx ;;
      sv_vcf_path) sv_idx=$idx ;;
    esac
  done

  if (( sample_idx < 0 || input_kind_idx < 0 )); then
    echo "ERROR: CSV must include sample and input_kind columns." >&2
    exit 1
  fi
  if (( sv_idx < 0 )); then
    if (( snv_vcf_idx >= 0 )); then
      sv_idx=$((snv_vcf_idx + 1))
      array_insert_value headers "$sv_idx" "sv_vcf_path"
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
      csv_split_line "$line" values
      if [[ "$added_sv_column" == "yes" ]]; then
        array_insert_value values "$sv_idx" ""
      fi
      while (( ${#values[@]} < ${#headers[@]} )); do
        values+=("")
      done

      row_count=$((row_count + 1))
      local sample input_kind current_sv matched_sv sv_path
      sample="$(trim_spaces "${values[$sample_idx]}")"
      input_kind="$(trim_spaces "${values[$input_kind_idx]}")"
      current_sv="$(trim_spaces "${values[$sv_idx]}")"

      if [[ "$input_kind" != "annotated_snv" ]]; then
        skipped_count=$((skipped_count + 1))
        csv_join_row "${values[@]}"
        continue
      fi
      if [[ -n "$current_sv" && "$overwrite_existing" != "yes" ]]; then
        skipped_count=$((skipped_count + 1))
        csv_join_row "${values[@]}"
        continue
      fi

      matched_sv="${update_sv_vcf_map[$sample]:-}"
      if [[ -z "$matched_sv" ]]; then
        local matched_prefix
        matched_prefix="$(find_best_match "$sample" "${update_sv_vcf_prefixes[@]}")"
        if [[ "$matched_prefix" == AMBIGUOUS:* || -z "$matched_prefix" ]]; then
          echo "WARNING: no unambiguous annotated SV VCF match for '$sample'; leaving sv_vcf_path empty." >&2
          missing_count=$((missing_count + 1))
          csv_join_row "${values[@]}"
          continue
        fi
        matched_sv="${update_sv_vcf_map[$matched_prefix]}"
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

generate_new_csv() {
  echo
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
    *)
      echo "ERROR: invalid selection: $input_choice" >&2
      exit 1
      ;;
  esac

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
    *)
      echo "ERROR: invalid selection: $phenotype_choice" >&2
      exit 1
      ;;
  esac

  read -r -p "Enter phenotype directory: " PHENOTYPE_DIR
  PHENOTYPE_DIR="$(trim_spaces "$PHENOTYPE_DIR")"
  [[ -d "$PHENOTYPE_DIR" ]] || { echo "ERROR: directory not found: $PHENOTYPE_DIR" >&2; exit 1; }

  read -r -p "Enter phenotype filename suffix (example: _note.txt or .hpo.txt): " PHENOTYPE_SUFFIX
  PHENOTYPE_SUFFIX="$(trim_spaces "$PHENOTYPE_SUFFIX")"
  [[ -n "$PHENOTYPE_SUFFIX" ]] || { echo "ERROR: phenotype suffix cannot be empty" >&2; exit 1; }

  read -r -p "Output CSV path [input_csv_generated.csv]: " OUT_CSV
  OUT_CSV="$(trim_spaces "$OUT_CSV")"
  OUT_CSV="${OUT_CSV:-input_csv_generated.csv}"

  read -r -p "Include optional age_of_onset column? [y/N]: " INCLUDE_AGE
  INCLUDE_AGE="$(trim_spaces "$INCLUDE_AGE")"
  INCLUDE_AGE="${INCLUDE_AGE,,}"
  if [[ "$INCLUDE_AGE" == "y" || "$INCLUDE_AGE" == "yes" ]]; then
    INCLUDE_AGE="yes"
  else
    INCLUDE_AGE="no"
  fi

  PROMPT_AGE="no"
  if [[ "$INCLUDE_AGE" == "yes" ]]; then
    read -r -p "Prompt for age per sample? [y/N]: " PROMPT_AGE
    PROMPT_AGE="$(trim_spaces "$PROMPT_AGE")"
    PROMPT_AGE="${PROMPT_AGE,,}"
    if [[ "$PROMPT_AGE" == "y" || "$PROMPT_AGE" == "yes" ]]; then
      PROMPT_AGE="yes"
    else
      PROMPT_AGE="no"
    fi
    if [[ "$PROMPT_AGE" == "yes" && ! -r /dev/tty ]]; then
      echo "WARNING: no interactive TTY available; disabling per-sample age prompts." >&2
      PROMPT_AGE="no"
    fi
  fi

  declare -A phenotype_map=()
  phenotype_prefixes=()
  load_prefix_map "$PHENOTYPE_DIR" "$PHENOTYPE_SUFFIX" phenotype_map phenotype_prefixes
  if [[ "${#phenotype_map[@]}" -eq 0 ]]; then
    echo "ERROR: no phenotype files matched suffix '$PHENOTYPE_SUFFIX' in $PHENOTYPE_DIR" >&2
    exit 1
  fi

  declare -A primary_map=()
  primary_prefixes=()
  declare -A secondary_map=()
  secondary_prefixes=()
  PRIMARY_HEADER=""
  PRIMARY_KIND_LABEL=""

  if [[ "$INPUT_KIND" == "annotated_snv" ]]; then
    read -r -p "Enter ANNOVAR TXT directory: " PRIMARY_DIR
    PRIMARY_DIR="$(trim_spaces "$PRIMARY_DIR")"
    [[ -d "$PRIMARY_DIR" ]] || { echo "ERROR: directory not found: $PRIMARY_DIR" >&2; exit 1; }
    read -r -p "Enter ANNOVAR TXT suffix [.hg38_multianno.txt]: " PRIMARY_SUFFIX
    PRIMARY_SUFFIX="$(trim_spaces "$PRIMARY_SUFFIX")"
    PRIMARY_SUFFIX="${PRIMARY_SUFFIX:-.hg38_multianno.txt}"

    read -r -p "Enter ANNOVAR VCF directory: " SECONDARY_DIR
    SECONDARY_DIR="$(trim_spaces "$SECONDARY_DIR")"
    [[ -d "$SECONDARY_DIR" ]] || { echo "ERROR: directory not found: $SECONDARY_DIR" >&2; exit 1; }
    read -r -p "Enter ANNOVAR VCF suffix [.hg38_multianno.vcf]: " SECONDARY_SUFFIX
    SECONDARY_SUFFIX="$(trim_spaces "$SECONDARY_SUFFIX")"
    SECONDARY_SUFFIX="${SECONDARY_SUFFIX:-.hg38_multianno.vcf}"
    read -r -p "Also include BAM/CRAM paths for SV/CNV/mito analysis? [y/N]: " INCLUDE_ALIGNMENT
    INCLUDE_ALIGNMENT="$(trim_spaces "$INCLUDE_ALIGNMENT")"
    INCLUDE_ALIGNMENT="${INCLUDE_ALIGNMENT,,}"
    if [[ "$INCLUDE_ALIGNMENT" == "y" || "$INCLUDE_ALIGNMENT" == "yes" ]]; then
      INCLUDE_ALIGNMENT="yes"
      read -r -p "Enter alignment directory: " ALIGNMENT_DIR
      ALIGNMENT_DIR="$(trim_spaces "$ALIGNMENT_DIR")"
      [[ -d "$ALIGNMENT_DIR" ]] || { echo "ERROR: directory not found: $ALIGNMENT_DIR" >&2; exit 1; }
      cat <<'MSG'
Choose alignment format:
  1) bam
  2) cram
MSG
      read -r -p "Select alignment format [1-2]: " alignment_choice
      alignment_choice="$(trim_spaces "$alignment_choice")"
      case "$alignment_choice" in
        1) ALIGNMENT_SUFFIX=".bam" ;;
        2) ALIGNMENT_SUFFIX=".cram" ;;
        *)
          echo "ERROR: invalid selection: $alignment_choice" >&2
          exit 1
          ;;
      esac
      declare -gA alignment_map=()
      declare -ga alignment_prefixes=()
      load_prefix_map "$ALIGNMENT_DIR" "$ALIGNMENT_SUFFIX" alignment_map alignment_prefixes
      if [[ "${#alignment_map[@]}" -eq 0 ]]; then
        echo "ERROR: no alignment files matched suffix '$ALIGNMENT_SUFFIX' in $ALIGNMENT_DIR" >&2
        exit 1
      fi
      read -r -p "Enter annotated SV VCF directory: " SV_VCF_DIR
      SV_VCF_DIR="$(trim_spaces "$SV_VCF_DIR")"
      [[ -d "$SV_VCF_DIR" ]] || { echo "ERROR: directory not found: $SV_VCF_DIR" >&2; exit 1; }
      read -r -p "Enter annotated SV VCF suffix [.hg38_multianno.vcf]: " SV_VCF_SUFFIX
      SV_VCF_SUFFIX="$(trim_spaces "$SV_VCF_SUFFIX")"
      SV_VCF_SUFFIX="${SV_VCF_SUFFIX:-.hg38_multianno.vcf}"
      declare -gA sv_vcf_map=()
      declare -ga sv_vcf_prefixes=()
      load_prefix_map "$SV_VCF_DIR" "$SV_VCF_SUFFIX" sv_vcf_map sv_vcf_prefixes
      if [[ "${#sv_vcf_map[@]}" -eq 0 ]]; then
        echo "ERROR: no annotated SV VCF files matched suffix '$SV_VCF_SUFFIX' in $SV_VCF_DIR" >&2
        exit 1
      fi
    else
      INCLUDE_ALIGNMENT="no"
    fi
    PRIMARY_KIND_LABEL="ANNOVAR TXT"
  else
    read -r -p "Enter input-data directory: " PRIMARY_DIR
    PRIMARY_DIR="$(trim_spaces "$PRIMARY_DIR")"
    [[ -d "$PRIMARY_DIR" ]] || { echo "ERROR: directory not found: $PRIMARY_DIR" >&2; exit 1; }
    case "$INPUT_KIND" in
      vcf_snv|vcf_sv)
        read -r -p "Enter VCF suffix [.vcf]: " PRIMARY_SUFFIX
        PRIMARY_SUFFIX="$(trim_spaces "$PRIMARY_SUFFIX")"
        PRIMARY_SUFFIX="${PRIMARY_SUFFIX:-.vcf}"
        ;;
      bam_ngs)
        PRIMARY_SUFFIX=".bam"
        ;;
      cram_ngs)
        PRIMARY_SUFFIX=".cram"
        ;;
    esac
    PRIMARY_KIND_LABEL="primary input"
  fi

  load_prefix_map "$PRIMARY_DIR" "$PRIMARY_SUFFIX" primary_map primary_prefixes
  if [[ "${#primary_map[@]}" -eq 0 ]]; then
    echo "ERROR: no ${PRIMARY_KIND_LABEL} files matched suffix '$PRIMARY_SUFFIX' in $PRIMARY_DIR" >&2
    exit 1
  fi

  if [[ "$INPUT_KIND" == "annotated_snv" ]]; then
    load_prefix_map "$SECONDARY_DIR" "$SECONDARY_SUFFIX" secondary_map secondary_prefixes
    if [[ "${#secondary_map[@]}" -eq 0 ]]; then
      echo "ERROR: no ANNOVAR VCF files matched suffix '$SECONDARY_SUFFIX' in $SECONDARY_DIR" >&2
      exit 1
    fi
  fi

  data_count=0
  match_count=0
  missing_count=0

  {
    echo "sample,input_kind,phenotype_path,phenotype_format,age_of_onset,snv_txt_path,snv_vcf_path,sv_vcf_path,vcf_path,alignment_path,alignment_index_path"

    for sample in "${primary_prefixes[@]}"; do
      data_count=$((data_count + 1))
      primary_path="${primary_map[$sample]}"
      phenotype_path="${phenotype_map[$sample]:-}"
      if [[ -z "$phenotype_path" ]]; then
        matched_phenotype="$(find_best_match "$sample" "${phenotype_prefixes[@]}")"
        if [[ "$matched_phenotype" == AMBIGUOUS:* || -z "$matched_phenotype" ]]; then
          echo "WARNING: no unambiguous phenotype match for '$sample'; skipping." >&2
          missing_count=$((missing_count + 1))
          continue
        fi
        phenotype_path="${phenotype_map[$matched_phenotype]}"
      fi

      snv_txt_path=""
      snv_vcf_path=""
      sv_vcf_path=""
      vcf_path=""
      alignment_path=""
      alignment_index_path=""

      if [[ "$INPUT_KIND" == "annotated_snv" ]]; then
        snv_txt_path="$primary_path"
        paired_vcf="${secondary_map[$sample]:-}"
        if [[ -z "$paired_vcf" ]]; then
          matched_vcf="$(find_best_match "$sample" "${secondary_prefixes[@]}")"
          if [[ "$matched_vcf" == AMBIGUOUS:* || -z "$matched_vcf" ]]; then
            echo "WARNING: no unambiguous ANNOVAR VCF match for '$sample'; skipping." >&2
            missing_count=$((missing_count + 1))
            continue
          fi
          paired_vcf="${secondary_map[$matched_vcf]}"
        fi
        snv_vcf_path="$paired_vcf"
        if [[ "${INCLUDE_ALIGNMENT:-no}" == "yes" ]]; then
          paired_alignment="${alignment_map[$sample]:-}"
          if [[ -z "$paired_alignment" ]]; then
            matched_alignment="$(find_best_match "$sample" "${alignment_prefixes[@]}")"
            if [[ "$matched_alignment" == AMBIGUOUS:* || -z "$matched_alignment" ]]; then
              echo "WARNING: no unambiguous alignment match for '$sample'; skipping." >&2
              missing_count=$((missing_count + 1))
              continue
            fi
            paired_alignment="${alignment_map[$matched_alignment]}"
          fi
          alignment_path="$paired_alignment"
          paired_sv_vcf="${sv_vcf_map[$sample]:-}"
          if [[ -z "$paired_sv_vcf" ]]; then
            matched_sv_vcf="$(find_best_match "$sample" "${sv_vcf_prefixes[@]}")"
            if [[ "$matched_sv_vcf" == AMBIGUOUS:* || -z "$matched_sv_vcf" ]]; then
              echo "WARNING: no unambiguous annotated SV VCF match for '$sample'; skipping." >&2
              missing_count=$((missing_count + 1))
              continue
            fi
            paired_sv_vcf="${sv_vcf_map[$matched_sv_vcf]}"
          fi
          sv_vcf_path="$paired_sv_vcf"
          if [[ "$alignment_path" == *.bam ]]; then
            if [[ -f "${alignment_path}.bai" ]]; then
              alignment_index_path="${alignment_path}.bai"
            else
              alt_index="${alignment_path%.bam}.bai"
              [[ -f "$alt_index" ]] && alignment_index_path="$alt_index"
            fi
          else
            if [[ -f "${alignment_path}.crai" ]]; then
              alignment_index_path="${alignment_path}.crai"
            else
              alt_index="${alignment_path%.cram}.crai"
              [[ -f "$alt_index" ]] && alignment_index_path="$alt_index"
            fi
          fi
        fi
      elif [[ "$INPUT_KIND" == "vcf_snv" || "$INPUT_KIND" == "vcf_sv" ]]; then
        vcf_path="$primary_path"
      else
        alignment_path="$primary_path"
        if [[ "$INPUT_KIND" == "bam_ngs" ]]; then
          if [[ -f "${alignment_path}.bai" ]]; then
            alignment_index_path="${alignment_path}.bai"
          else
            alt_index="${alignment_path%.bam}.bai"
            [[ -f "$alt_index" ]] && alignment_index_path="$alt_index"
          fi
        else
          if [[ -f "${alignment_path}.crai" ]]; then
            alignment_index_path="${alignment_path}.crai"
          else
            alt_index="${alignment_path%.cram}.crai"
            [[ -f "$alt_index" ]] && alignment_index_path="$alt_index"
          fi
        fi
      fi

      age_value=""
      if [[ "$PROMPT_AGE" == "yes" ]]; then
        age_value="$(prompt_age "$sample")"
      fi

      printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$sample" \
        "$INPUT_KIND" \
        "$phenotype_path" \
        "$PHENOTYPE_FORMAT" \
        "$age_value" \
        "$snv_txt_path" \
        "$snv_vcf_path" \
        "$sv_vcf_path" \
        "$vcf_path" \
        "$alignment_path" \
        "$alignment_index_path"
      match_count=$((match_count + 1))
    done
  } > "$OUT_CSV"

  echo
  echo "Done."
  echo "  Input files scanned : $data_count"
  echo "  Matched rows        : $match_count"
  echo "  Missing pair/match  : $missing_count"
  echo "  Output CSV          : $OUT_CSV"

  if [[ "$match_count" -eq 0 ]]; then
    echo "ERROR: zero matched samples. Check suffixes and prefixes." >&2
    exit 1
  fi
}

main() {
  echo "PipeVar_annotated_snv unified CSV generator"
  echo
  cat <<'MSG'
Choose action:
  1) Generate new CSV
  2) Update existing annotated_snv CSV with annotated SV VCF paths
MSG
  read -r -p "Select action [1-2]: " action_choice
  action_choice="$(trim_spaces "$action_choice")"
  case "$action_choice" in
    1|"") generate_new_csv ;;
    2) update_existing_annotated_snv_csv ;;
    *)
      echo "ERROR: invalid selection: $action_choice" >&2
      exit 1
      ;;
  esac
}

main "$@"
