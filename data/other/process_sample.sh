#!/bin/bash

# Define paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGS_DIR="${SCRIPT_DIR}/logs"
SLURM_LOGS_DIR="${LOGS_DIR}/slurm"
RESULTS_DIR="${SCRIPT_DIR}/results"
ENV_DIR="/projects/microsplit/env/envSimpleTrimming"

# Function to get file size in GB
get_file_size_gb() {
    local file="$1"
    local size_bytes=$(stat -c%s "$file")
    local size_gb=$(echo "scale=2; $size_bytes/1024/1024/1024" | bc)
    echo "$size_gb"
}

# Function to estimate required memory based on input file size
estimate_required_memory() {
    local r1_input="$1"
    local r2_input="$2"
    local r1_size=$(get_file_size_gb "$r1_input")
    local r2_size=$(get_file_size_gb "$r2_input")
    local total_size=$(echo "$r1_size + $r2_size" | bc)
    
    # Estimate memory needed (2x the input size as a safety factor)
    local required_mem=$(echo "scale=0; $total_size * 2" | bc)
    
    # Round up to nearest GB and add 4GB base memory
    required_mem=$(echo "scale=0; ($required_mem + 4) / 1" | bc)
    
    # Ensure minimum of 8GB
    if [ "$required_mem" -lt 8 ]; then
        required_mem=8
    fi
    
    echo "${required_mem}G"
}

# Function to log messages
log() {
    local level=$1
    local message=$2
    local run_id="$3"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] [${level}] ${message}" | tee -a "${LOGS_DIR}/${run_id}_${sample_name}.log"
}

# Function to check if a command succeeded
check_status() {
    if [ $? -ne 0 ]; then
        log "ERROR" "$1 failed" "$RUN_ID"
        exit 1
    fi
    log "INFO" "$1 completed successfully" "$RUN_ID"
}

# Function to perform TSO trimming with resource management
trim_tso() {
    local r1_input="$1"
    local r2_input="$2"
    local r1_output="$3"
    local r2_output="$4"
    local sample_name="$5"
    local output_dir="$6"
    local run_id="$7"

    log "INFO" "Step 1: Trimming TSO sequences..." "$run_id"
    log "INFO" "Using ${SLURM_CPUS_PER_TASK} CPUs" "$run_id"
    
    # Monitor memory usage
    local mem_usage_file="${output_dir}/${sample_name}_tso_mem_usage.txt"
    (while true; do
        echo "$(date '+%Y-%m-%d %H:%M:%S') $(free -h | grep Mem | awk '{print $3 "/" $2}')"
        sleep 30
    done) > "$mem_usage_file" &
    local mem_monitor_pid=$!

    cutadapt -j ${SLURM_CPUS_PER_TASK} \
        -g "AAGCAGTGGTATCAACGCAGAGTGAATGGG; min_overlap=6; max_errors=0.2" \
        -g "CAGAGTGAATGGG; min_overlap=6; max_errors=0.2" \
        --pair-filter=both \
        -m 20: \
        --too-short-output "${output_dir}/${sample_name}_R1_too_short.fastq.gz" \
        --too-short-paired-output "${output_dir}/${sample_name}_R2_too_short.fastq.gz" \
        -o "${r1_output}" \
        -p "${r2_output}" \
        "${r1_input}" "${r2_input}" \
        --report=full \
        --json "${output_dir}/${sample_name}_stats.json"
    
    kill $mem_monitor_pid
    check_status "TSO trimming" "$run_id"
}

# Function to perform initial fastp trimming
trim_fastp_initial() {
    local r1_input="$1"
    local r2_input="$2"
    local r1_output="$3"
    local r2_output="$4"
    local unpaired1="$5"
    local unpaired2="$6"
    local sample_name="$7"
    local output_dir="$8"
    local run_id="$9"

    log "INFO" "Step 2: Performing initial fastp trimming..." "$run_id"
    fastp \
        -i "${r1_input}" \
        -I "${r2_input}" \
        -o "${r1_output}" \
        -O "${r2_output}" \
        --html "${output_dir}/${sample_name}_report.html" \
        --json "${output_dir}/${sample_name}_report.json" \
        --report_title "microSplit Initial Fastp Report - ${sample_name}" \
        --compression 4 \
        --verbose \
        --unpaired1 "${unpaired1}" \
        --unpaired2 "${unpaired2}" \
        --length_required 91 \
        --dont_overwrite \
        --trim_front1 0 \
        --trim_front2 0 \
        --trim_tail1 0 \
        --trim_tail2 0 \
        --trim_poly_g \
        --poly_g_min_len 10 \
        --trim_poly_x \
        --poly_x_min_len 12 \
        --detect_adapter_for_pe \
        --adapter_sequence=ATCTCGTATGCCGTCTTCTGCTTGA \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    check_status "Initial fastp trimming" "$run_id"
}

# Function to perform polyA trimming
trim_polya() {
    local r1_input="$1"
    local r2_input="$2"
    local r1_output="$3"
    local r2_output="$4"
    local sample_name="$5"
    local output_dir="$6"
    local run_id="$7"

    log "INFO" "Step 3: Trimming polyA sequences..." "$run_id"
    cutadapt -j ${SLURM_CPUS_PER_TASK} \
        -a "A{12}; min_overlap=12; max_errors=0.2" \
        --pair-filter=both \
        -m 20: \
        --too-short-output "${output_dir}/${sample_name}_R1_too_short.fastq.gz" \
        --too-short-paired-output "${output_dir}/${sample_name}_R2_too_short.fastq.gz" \
        -o "${r1_output}" \
        -p "${r2_output}" \
        "${r1_input}" "${r2_input}" \
        --report=full \
        --json "${output_dir}/${sample_name}_stats.json"
    check_status "PolyA trimming" "$run_id"
}

# Function to perform specific adapter trimming
trim_specific_adapter() {
    local r1_input="$1"
    local r2_input="$2"
    local r1_output="$3"
    local r2_output="$4"
    local sample_name="$5"
    local output_dir="$6"
    local run_id="$7"

    log "INFO" "Step 4: Trimming specific adapter CCACAGTCTCAAGCAC..." "$run_id"
    cutadapt -j ${SLURM_CPUS_PER_TASK} \
        -a "CCACAGTCTCAAGCAC; min_overlap=6; max_errors=0.1" \
        --pair-filter=both \
        -m 20: \
        --too-short-output "${output_dir}/${sample_name}_R1_too_short.fastq.gz" \
        --too-short-paired-output "${output_dir}/${sample_name}_R2_too_short.fastq.gz" \
        -o "${r1_output}" \
        -p "${r2_output}" \
        "${r1_input}" "${r2_input}" \
        --report=full \
        --json "${output_dir}/${sample_name}_stats.json"
    check_status "Specific adapter trimming" "$run_id"
}

# Function to perform linker trimming
trim_linkers() {
    local r1_input="$1"
    local r2_input="$2"
    local r1_output="$3"
    local r2_output="$4"
    local sample_name="$5"
    local output_dir="$6"
    local run_id="$7"

    log "INFO" "Step 5: Trimming linkers and specific adapters..." "$run_id"
    cutadapt -j ${SLURM_CPUS_PER_TASK} \
        -a "CCACAGTCTCAAGCACGTGGAT; min_overlap=6; max_errors=0.2" \
        -a "AGTCGTACGCCGATGCGAAACATCGGCCAC; min_overlap=6; max_errors=0.2" \
        -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA; min_overlap=6; max_errors=0.2" \
        --pair-filter=both \
        -m 20: \
        --too-short-output "${output_dir}/${sample_name}_R1_too_short.fastq.gz" \
        --too-short-paired-output "${output_dir}/${sample_name}_R2_too_short.fastq.gz" \
        -o "${r1_output}" \
        -p "${r2_output}" \
        "${r1_input}" "${r2_input}" \
        --report=full \
        --json "${output_dir}/${sample_name}_stats.json"
    check_status "Linker trimming" "$run_id"
}

# Function to perform final fastp trimming
trim_fastp_final() {
    local r1_input="$1"
    local r2_input="$2"
    local r1_output="$3"
    local r2_output="$4"
    local sample_name="$5"
    local output_dir="$6"
    local run_id="$7"

    log "INFO" "Step 6: Performing final fastp trimming..." "$run_id"
    fastp \
        -i "${r1_input}" \
        -I "${r2_input}" \
        -o "${r1_output}" \
        -O "${r2_output}" \
        --trim_front1 10 \
        --trim_front2 0 \
        --trim_tail1 16 \
        --trim_tail2 0 \
        --length_required 25 \
        --detect_adapter_for_pe \
        --adapter_sequence=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
        --adapter_sequence=CCACAGTCTCAAGCACGTGGAT \
        --adapter_sequence=AGTCGTACGCCGATGCGAAACATCGGCCAC \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --html "${output_dir}/${sample_name}_report.html" \
        --json "${output_dir}/${sample_name}_report.json" \
        --report_title "microSplit Final Fastp Report - ${sample_name}" \
        --compression 4 \
        --verbose
    check_status "Final fastp trimming" "$run_id"
}

# Main function to process a sample
main() {
    if [ $# -ne 5 ]; then
        echo "Usage: $0 <sample_name> <R1_input> <R2_input> <output_dir> <run_id>"
        exit 1
    fi

    local sample_name="$1"
    local r1_input="$2"
    local r2_input="$3"
    local output_dir="$4"
    local RUN_ID="$5"
    
    # Log system information
    log "INFO" "Starting processing for sample: ${sample_name}" "$RUN_ID"
    log "INFO" "Output directory: ${output_dir}" "$RUN_ID"
    log "INFO" "System information:" "$RUN_ID"
    log "INFO" "CPUs available: ${SLURM_CPUS_PER_TASK}" "$RUN_ID"
    log "INFO" "Memory available: ${SLURM_MEM_PER_NODE}" "$RUN_ID"
    log "INFO" "Node: ${SLURM_NODELIST}" "$RUN_ID"
    
    # Estimate required memory
    local estimated_mem=$(estimate_required_memory "$r1_input" "$r2_input")
    log "INFO" "Estimated required memory: ${estimated_mem}" "$RUN_ID"
    
    # Verify input files exist
    if [ ! -f "${r1_input}" ] || [ ! -f "${r2_input}" ]; then
        log "ERROR" "Input files not found for sample ${sample_name}" "$RUN_ID"
        exit 1
    fi
    
    # Create directory structure
    local dirs=(
        "${output_dir}"
        "${output_dir}/1_tso_trimming"
        "${output_dir}/2_fastp_initial"
        "${output_dir}/3_polya_trimming"
        "${output_dir}/4_specific_adapter"
        "${output_dir}/5_linker_trimming"
        "${output_dir}/6_fastp_final"
        "${output_dir}/resource_usage"  # Pour les fichiers de monitoring
    )
    
    for dir in "${dirs[@]}"; do
        mkdir -p "${dir}"
        log "DEBUG" "Created directory: ${dir}" "$RUN_ID"
    done
    
    # Create resource usage summary file
    local resource_summary="${output_dir}/resource_usage/${sample_name}_resource_usage.txt"
    echo "Resource Usage Summary for ${sample_name}" > "$resource_summary"
    echo "=======================================" >> "$resource_summary"
    echo "Start time: $(date)" >> "$resource_summary"
    echo "Input files:" >> "$resource_summary"
    echo "  R1: ${r1_input} ($(get_file_size_gb "$r1_input") GB)" >> "$resource_summary"
    echo "  R2: ${r2_input} ($(get_file_size_gb "$r2_input") GB)" >> "$resource_summary"
    echo "SLURM resources:" >> "$resource_summary"
    echo "  CPUs: ${SLURM_CPUS_PER_TASK}" >> "$resource_summary"
    echo "  Memory: ${SLURM_MEM_PER_NODE}" >> "$resource_summary"
    echo "  Node: ${SLURM_NODELIST}" >> "$resource_summary"
    echo "Estimated required memory: ${estimated_mem}" >> "$resource_summary"
    
    # Define intermediate files with new directory structure
    local r1_tso="${output_dir}/1_tso_trimming/${sample_name}_R1.fastq.gz"
    local r2_tso="${output_dir}/1_tso_trimming/${sample_name}_R2.fastq.gz"
    local r1_fastp="${output_dir}/2_fastp_initial/${sample_name}_R1.fastq.gz"
    local r2_fastp="${output_dir}/2_fastp_initial/${sample_name}_R2.fastq.gz"
    local r1_polya="${output_dir}/3_polya_trimming/${sample_name}_R1.fastq.gz"
    local r2_polya="${output_dir}/3_polya_trimming/${sample_name}_R2.fastq.gz"
    local r1_specific="${output_dir}/4_specific_adapter/${sample_name}_R1.fastq.gz"
    local r2_specific="${output_dir}/4_specific_adapter/${sample_name}_R2.fastq.gz"
    local r1_linker="${output_dir}/5_linker_trimming/${sample_name}_R1.fastq.gz"
    local r2_linker="${output_dir}/5_linker_trimming/${sample_name}_R2.fastq.gz"
    local r1_final="${output_dir}/6_fastp_final/${sample_name}_R1.fastq.gz"
    local r2_final="${output_dir}/6_fastp_final/${sample_name}_R2.fastq.gz"
    
    # Execute trimming steps in sequence with resource monitoring
    trim_tso "${r1_input}" "${r2_input}" "${r1_tso}" "${r2_tso}" "${sample_name}" "${output_dir}/1_tso_trimming" "$RUN_ID"
    trim_fastp_initial "${r1_tso}" "${r2_tso}" "${r1_fastp}" "${r2_fastp}" \
        "${output_dir}/2_fastp_initial/${sample_name}_R1_unpaired.fastq.gz" \
        "${output_dir}/2_fastp_initial/${sample_name}_R2_unpaired.fastq.gz" \
        "${sample_name}" "${output_dir}/2_fastp_initial" "$RUN_ID"
    trim_polya "${r1_fastp}" "${r2_fastp}" "${r1_polya}" "${r2_polya}" "${sample_name}" "${output_dir}/3_polya_trimming" "$RUN_ID"
    trim_specific_adapter "${r1_polya}" "${r2_polya}" "${r1_specific}" "${r2_specific}" "${sample_name}" "${output_dir}/4_specific_adapter" "$RUN_ID"
    trim_linkers "${r1_specific}" "${r2_specific}" "${r1_linker}" "${r2_linker}" "${sample_name}" "${output_dir}/5_linker_trimming" "$RUN_ID"
    trim_fastp_final "${r1_linker}" "${r2_linker}" "${r1_final}" "${r2_final}" "${sample_name}" "${output_dir}/6_fastp_final" "$RUN_ID"
    
    # Create symbolic links to final files in the main directory
    ln -sf "6_fastp_final/${sample_name}_R1.fastq.gz" "${output_dir}/${sample_name}_R1_final.fastq.gz"
    ln -sf "6_fastp_final/${sample_name}_R2.fastq.gz" "${output_dir}/${sample_name}_R2_final.fastq.gz"
    
    # Update resource summary with end time
    echo "End time: $(date)" >> "$resource_summary"
    echo "=======================================" >> "$resource_summary"
    
    log "INFO" "Processing completed for sample ${sample_name}" "$RUN_ID"
    log "INFO" "Resource usage summary available in: ${resource_summary}" "$RUN_ID"
}

# Run main function
main "$@" 