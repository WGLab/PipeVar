#!/usr/bin/env bash
set -euo pipefail

# Build PipeVar input CSV by matching data files and note files on shared prefixes.

print_formats() {
  cat <<'MSG'
Choose data file format:
  1) bam
  2) cram
  3) vcf
  4) vcf.gz
MSG
}

normalize_data_prefix() {
  local filename="$1"
  case "$DATA_FORMAT" in
    bam)    echo "${filename%.bam}" ;;
    cram)   echo "${filename%.cram}" ;;
    vcf)    echo "${filename%.vcf}" ;;
    vcf.gz) echo "${filename%.vcf.gz}" ;;
    *)      echo "$filename" ;;
  esac
}

normalize_note_prefix() {
  local filename="$1"
  echo "${filename%$NOTE_SUFFIX}"
}

trim_spaces() {
  local s="$1"
  # trim leading/trailing whitespace
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

collect_data_files() {
  case "$DATA_FORMAT" in
    vcf.gz)
      find "$DATA_DIR" -maxdepth 1 -type f -name "*.vcf.gz" | sort
      ;;
    *)
      find "$DATA_DIR" -maxdepth 1 -type f -name "*.${DATA_FORMAT}" | sort
      ;;
  esac
}

normalize_for_match() {
  # Lowercase and remove non-alphanumeric characters for rough matching.
  local s="$1"
  s="${s,,}"
  s="${s//[^a-z0-9]/}"
  echo "$s"
}

match_score() {
  # Return a numeric score for how well note_prefix matches sample_prefix.
  # Higher score is better.
  local sample_prefix="$1"
  local note_prefix="$2"
  local s_norm n_norm
  s_norm="$(normalize_for_match "$sample_prefix")"
  n_norm="$(normalize_for_match "$note_prefix")"

  if [[ -z "$s_norm" || -z "$n_norm" ]]; then
    echo 0
    return
  fi

  if [[ "$sample_prefix" == "$note_prefix" ]]; then
    echo 1000
    return
  fi

  if [[ "$s_norm" == "$n_norm" ]]; then
    echo 900
    return
  fi

  if [[ "$s_norm" == *"$n_norm"* || "$n_norm" == *"$s_norm"* ]]; then
    local len_diff=${#s_norm}
    if (( ${#n_norm} > len_diff )); then
      len_diff=$(( ${#n_norm} - ${#s_norm} ))
    else
      len_diff=$(( ${#s_norm} - ${#n_norm} ))
    fi
    echo $((700 - len_diff))
    return
  fi

  # Fallback: shared leading token before common separators.
  local s_token n_token
  s_token="$(echo "$sample_prefix" | awk -F'[_ .-]' '{print tolower($1)}')"
  n_token="$(echo "$note_prefix" | awk -F'[_ .-]' '{print tolower($1)}')"
  if [[ -n "$s_token" && "$s_token" == "$n_token" ]]; then
    echo 400
    return
  fi

  echo 0
}

find_best_note_match() {
  local sample_prefix="$1"
  local sample_norm prefix prefix_norm

  # First pass: normalized exact match (case/separator-insensitive).
  # This prevents false "ambiguous" for truly matching names like:
  #   SAMPLE-01.vcf.gz <-> sample_01.txt
  sample_norm="$(normalize_for_match "$sample_prefix")"
  local norm_hits=()
  for prefix in "${note_prefixes[@]}"; do
    prefix_norm="$(normalize_for_match "$prefix")"
    if [[ -n "$sample_norm" && "$sample_norm" == "$prefix_norm" ]]; then
      norm_hits+=("$prefix")
    fi
  done
  if (( ${#norm_hits[@]} == 1 )); then
    echo "${norm_hits[0]}"
    return
  fi
  if (( ${#norm_hits[@]} > 1 )); then
    # If one candidate is an exact string match, use it.
    for prefix in "${norm_hits[@]}"; do
      if [[ "$prefix" == "$sample_prefix" ]]; then
        echo "$prefix"
        return
      fi
    done

    # Then prefer case-insensitive exact raw-stem match.
    for prefix in "${norm_hits[@]}"; do
      if [[ "${prefix,,}" == "${sample_prefix,,}" ]]; then
        echo "$prefix"
        return
      fi
    done

    # Then prefer the candidate with smallest raw-length difference to sample.
    # This resolves common cases like sample=S1 with notes S1.txt and S-1.txt.
    local best_len_prefix=""
    local best_len_diff=999999
    local tie_len=0
    local diff
    for prefix in "${norm_hits[@]}"; do
      diff=${#prefix}
      if (( ${#sample_prefix} > diff )); then
        diff=$(( ${#sample_prefix} - ${#prefix} ))
      else
        diff=$(( ${#prefix} - ${#sample_prefix} ))
      fi
      if (( diff < best_len_diff )); then
        best_len_diff=$diff
        best_len_prefix="$prefix"
        tie_len=0
      elif (( diff == best_len_diff )); then
        tie_len=1
      fi
    done
    if (( tie_len == 0 )) && [[ -n "$best_len_prefix" ]]; then
      echo "$best_len_prefix"
      return
    fi

    echo "AMBIGUOUS:$(IFS='|'; echo "${norm_hits[*]}")"
    return
  fi

  local best_prefix=""
  local best_score=0
  local best_tie=0
  local prefix score

  for prefix in "${note_prefixes[@]}"; do
    score="$(match_score "$sample_prefix" "$prefix")"
    if (( score > best_score )); then
      best_score=$score
      best_prefix="$prefix"
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

main() {
  echo "PipeVar input CSV generator"
  echo

  read -r -p "Enter data-file directory (bam/cram/vcf/vcf.gz): " DATA_DIR
  DATA_DIR="$(trim_spaces "$DATA_DIR")"
  if [[ ! -d "$DATA_DIR" ]]; then
    echo "ERROR: directory not found: $DATA_DIR" >&2
    exit 1
  fi

  print_formats
  read -r -p "Select format [1-4]: " format_choice
  format_choice="$(trim_spaces "$format_choice")"
  case "$format_choice" in
    1) DATA_FORMAT="bam" ;;
    2) DATA_FORMAT="cram" ;;
    3) DATA_FORMAT="vcf" ;;
    4) DATA_FORMAT="vcf.gz" ;;
    *)
      echo "ERROR: invalid selection: $format_choice" >&2
      exit 1
      ;;
  esac

  read -r -p "Enter medical-note text directory: " NOTE_DIR
  NOTE_DIR="$(trim_spaces "$NOTE_DIR")"
  if [[ ! -d "$NOTE_DIR" ]]; then
    echo "ERROR: directory not found: $NOTE_DIR" >&2
    exit 1
  fi

  read -r -p "Enter note filename suffix (example: _note.txt or .txt): " NOTE_SUFFIX
  NOTE_SUFFIX="$(trim_spaces "$NOTE_SUFFIX")"
  if [[ -z "$NOTE_SUFFIX" ]]; then
    echo "ERROR: note suffix cannot be empty" >&2
    exit 1
  fi

  read -r -p "Output CSV path [input_csv_generated.csv]: " OUT_CSV
  OUT_CSV="$(trim_spaces "$OUT_CSV")"
  OUT_CSV="${OUT_CSV:-input_csv_generated.csv}"

  read -r -p "Include optional age_of_onset column for prioritization flows? (values: xd/xm/xy or integer years) [y/N]: " INCLUDE_AGE
  INCLUDE_AGE="$(trim_spaces "$INCLUDE_AGE")"
  INCLUDE_AGE="${INCLUDE_AGE,,}"
  if [[ "$INCLUDE_AGE" == "y" || "$INCLUDE_AGE" == "yes" ]]; then
    INCLUDE_AGE="yes"
  else
    INCLUDE_AGE="no"
  fi

  PROMPT_AGE_PER_SAMPLE="no"
  if [[ "$INCLUDE_AGE" == "yes" ]]; then
    read -r -p "Prompt for age per sample while building CSV? [y/N]: " PROMPT_AGE_PER_SAMPLE
    PROMPT_AGE_PER_SAMPLE="$(trim_spaces "$PROMPT_AGE_PER_SAMPLE")"
    PROMPT_AGE_PER_SAMPLE="${PROMPT_AGE_PER_SAMPLE,,}"
    if [[ "$PROMPT_AGE_PER_SAMPLE" == "y" || "$PROMPT_AGE_PER_SAMPLE" == "yes" ]]; then
      PROMPT_AGE_PER_SAMPLE="yes"
    else
      PROMPT_AGE_PER_SAMPLE="no"
    fi
    if [[ "$PROMPT_AGE_PER_SAMPLE" == "yes" && ! -r /dev/tty ]]; then
      echo "WARNING: no interactive TTY available; disabling per-sample age prompts." >&2
      PROMPT_AGE_PER_SAMPLE="no"
    fi
  fi

  # Build note prefix -> full path map
  declare -A note_map
  note_prefixes=()
  while IFS= read -r note_path; do
    note_name="$(basename "$note_path")"
    note_prefix="$(normalize_note_prefix "$note_name")"

    # only keep files that actually end with the requested suffix
    if [[ "$note_name" == *"$NOTE_SUFFIX" ]] && [[ -n "$note_prefix" ]]; then
      if [[ -n "${note_map[$note_prefix]:-}" ]]; then
        echo "WARNING: duplicate note prefix '$note_prefix' found; keeping first:" >&2
        echo "         ${note_map[$note_prefix]}" >&2
        echo "         skipped: $note_path" >&2
      else
        note_map[$note_prefix]="$note_path"
        note_prefixes+=("$note_prefix")
      fi
    fi
  done < <(find "$NOTE_DIR" -maxdepth 1 -type f -name "*${NOTE_SUFFIX}" | sort)

  if [[ "${#note_map[@]}" -eq 0 ]]; then
    echo "ERROR: no note files matched suffix '$NOTE_SUFFIX' in $NOTE_DIR" >&2
    exit 1
  fi

  data_count=0
  match_count=0
  missing_note_count=0

  {
    if [[ "$INCLUDE_AGE" == "yes" ]]; then
      echo "sample,file_path,note_path,age_of_onset"
    else
      echo "sample,file_path,note_path"
    fi

    while IFS= read -r data_path; do
      data_count=$((data_count + 1))
      data_name="$(basename "$data_path")"
      sample="$(normalize_data_prefix "$data_name")"

      if [[ -n "${note_map[$sample]:-}" ]]; then
        match_count=$((match_count + 1))
        if [[ "$INCLUDE_AGE" == "yes" ]]; then
          age_out=""
          if [[ "$PROMPT_AGE_PER_SAMPLE" == "yes" ]]; then
            while true; do
              printf "Age for sample '%s' [empty|e.g. 10d,9m,7y,7]: " "$sample" > /dev/tty
              read -r age_in < /dev/tty
              age_in="$(trim_spaces "$age_in")"
              if is_valid_compact_age "$age_in"; then
                age_out="$(normalize_compact_age "$age_in")"
                break
              fi
              echo "Invalid age format. Use empty, <int>, or <int><d|m|y>." > /dev/tty
            done
          fi
          echo "${sample},${data_path},${note_map[$sample]},${age_out}"
        else
          echo "${sample},${data_path},${note_map[$sample]}"
        fi
      else
        matched_prefix="$(find_best_note_match "$sample")"
        if [[ "$matched_prefix" == AMBIGUOUS:* ]]; then
          missing_note_count=$((missing_note_count + 1))
          candidates="${matched_prefix#AMBIGUOUS:}"
          echo "WARNING: ambiguous rough match for '$sample' (${data_name}); candidates: ${candidates}; skipping." >&2
        elif [[ -n "$matched_prefix" && -n "${note_map[$matched_prefix]:-}" ]]; then
          match_count=$((match_count + 1))
          echo "INFO: rough matched '$sample' -> '$matched_prefix'" >&2
          if [[ "$INCLUDE_AGE" == "yes" ]]; then
            age_out=""
            if [[ "$PROMPT_AGE_PER_SAMPLE" == "yes" ]]; then
              while true; do
                printf "Age for sample '%s' [empty|e.g. 10d,9m,7y,7]: " "$sample" > /dev/tty
                read -r age_in < /dev/tty
                age_in="$(trim_spaces "$age_in")"
                if is_valid_compact_age "$age_in"; then
                  age_out="$(normalize_compact_age "$age_in")"
                  break
                fi
                echo "Invalid age format. Use empty, <int>, or <int><d|m|y>." > /dev/tty
              done
            fi
            echo "${sample},${data_path},${note_map[$matched_prefix]},${age_out}"
          else
            echo "${sample},${data_path},${note_map[$matched_prefix]}"
          fi
        else
          missing_note_count=$((missing_note_count + 1))
          echo "WARNING: no note match for sample prefix '$sample' (${data_name})" >&2
        fi
      fi
    done < <(collect_data_files)
  } > "$OUT_CSV"

  echo
  echo "Done."
  echo "  Data files scanned : $data_count"
  echo "  Matched rows       : $match_count"
  echo "  Missing note match : $missing_note_count"
  echo "  Output CSV         : $OUT_CSV"

  if [[ "$match_count" -eq 0 ]]; then
    echo "ERROR: zero matched samples. Check prefixes and note suffix." >&2
    exit 1
  fi
}

main "$@"
