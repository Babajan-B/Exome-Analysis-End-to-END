#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import logging
import threading
import shutil
import time
import json
import gzip
import re

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Reference genome (hg19)
REFERENCE_GENOME = os.path.join(os.getcwd(), 'reference', 'hg19.fa')

def ensure_reference_indices(result_folder, analysis_id):
    """Ensure the reference genome has all required indices."""
    logger.info(f"Checking reference genome indices...")
    
    # Check for BWA index files
    bwa_index_files = [
        f"{REFERENCE_GENOME}.amb",
        f"{REFERENCE_GENOME}.ann", 
        f"{REFERENCE_GENOME}.bwt",
        f"{REFERENCE_GENOME}.pac",
        f"{REFERENCE_GENOME}.sa"
    ]
    
    # Check if any index files are missing
    missing_files = [f for f in bwa_index_files if not os.path.exists(f)]
    
    if missing_files:
        logger.info(f"Missing BWA index files: {', '.join([os.path.basename(f) for f in missing_files])}")
        logger.info(f"Creating BWA index for reference genome... (this may take a while)")
        
        # Update progress to indicate indexing
        update_progress(analysis_id, "align", "running", percent=5, 
                        note="Creating reference genome indices (this may take 15-30 minutes)")
        
        # Create BWA index
        cmd = f"bwa index {REFERENCE_GENOME}"
        logger.info(f"Running command: {cmd}")
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, universal_newlines=True)
        
        # Monitor output in real-time
        while process.poll() is None:
            stderr_line = process.stderr.readline()
            if stderr_line:
                logger.info(f"BWA index: {stderr_line.strip()}")
            stdout_line = process.stdout.readline()
            if stdout_line:
                logger.info(f"BWA index: {stdout_line.strip()}")
            time.sleep(0.1)
        
        # Get remaining output
        stdout, stderr = process.communicate()
        if stdout:
            logger.info(f"BWA index stdout: {stdout}")
        if stderr:
            logger.info(f"BWA index stderr: {stderr}")
        
        if process.returncode != 0:
            logger.error(f"BWA index creation failed with return code {process.returncode}")
            update_progress(analysis_id, "align", "error", 
                           note="Failed to create reference genome indices")
            return False
        
        logger.info(f"BWA index created successfully")
    else:
        logger.info(f"Reference genome BWA indices are already present")
    
    # Check for samtools index (.fai)
    if not os.path.exists(f"{REFERENCE_GENOME}.fai"):
        logger.info(f"Creating samtools index for reference genome...")
        
        cmd = f"samtools faidx {REFERENCE_GENOME}"
        logger.info(f"Running command: {cmd}")
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Samtools faidx failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
            update_progress(analysis_id, "align", "error", 
                           note="Failed to create samtools index for reference genome")
            return False
        
        logger.info(f"Samtools index created successfully")
    else:
        logger.info(f"Reference genome samtools index (.fai) is already present")
    
    # Check for GATK dictionary (.dict)
    dict_file = REFERENCE_GENOME.replace('.fa', '.dict')
    if not os.path.exists(dict_file):
        logger.info(f"Creating sequence dictionary for reference genome...")
        
        # Use standalone Picard JAR with native syntax
        picard_jar = os.path.join(os.getcwd(), 'tools', 'picard.jar')
        cmd = f"java -jar {picard_jar} CreateSequenceDictionary \
              R={REFERENCE_GENOME} \
              O={dict_file} \
              VALIDATION_STRINGENCY=LENIENT"
        logger.info(f"Running command: {cmd}")
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, universal_newlines=True)
        
        # Monitor output in real-time
        while process.poll() is None:
            stderr_line = process.stderr.readline()
            if stderr_line:
                logger.info(f"Picard CreateSequenceDictionary: {stderr_line.strip()}")
            stdout_line = process.stdout.readline()
            if stdout_line:
                logger.info(f"Picard CreateSequenceDictionary: {stdout_line.strip()}")
            time.sleep(0.1)
        
        # Get remaining output
        stdout, stderr = process.communicate()
        if stdout:
            logger.info(f"Picard stdout: {stdout}")
        if stderr:
            logger.info(f"Picard stderr: {stderr}")
        
        if process.returncode != 0:
            logger.error(f"Sequence dictionary creation failed with return code {process.returncode}")
            update_progress(analysis_id, "align", "error", 
                           note="Failed to create reference genome dictionary")
            return False
        
        logger.info(f"Sequence dictionary created successfully")
    else:
        logger.info(f"Reference genome sequence dictionary (.dict) is already present")
    
    return True

def verify_bam_for_gatk(bam_file, analysis_id, result_folder):
    """Verify that a BAM file is compatible with GATK before variant calling."""
    logger.info(f"Verifying BAM file compatibility with GATK: {bam_file}")
    
    # Check if BAM has read groups
    cmd = f"samtools view -H {bam_file} | grep '^@RG'"
    logger.info(f"Running command: {cmd}")
    
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    has_read_groups = False
    if process.returncode == 0 and stdout:
        has_read_groups = True
        logger.info(f"BAM file has read groups: {stdout.decode() if isinstance(stdout, bytes) else stdout}")
    else:
        logger.warning(f"BAM file does not have read groups. Adding them now.")
        
        # Add read groups to BAM file
        sample_name = f"sample_{analysis_id}"
        fixed_bam = f"{bam_file}.with_rg.bam"
        
        # Create read group information
        platform = "ILLUMINA"
        library = f"lib_{analysis_id}"
        read_group = f"ID:{analysis_id}\\tSM:{analysis_id}\\tPL:{platform}\\tLB:{library}\\tPU:unit1"
        
        # Use samtools to add read group
        cmd = f"samtools addreplacerg -r '{read_group}' -o {fixed_bam} {bam_file}"
        logger.info(f"Running command: {cmd}")
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Failed to add read groups: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
            update_progress(analysis_id, "variants", "error", note="Failed to add read groups to BAM file")
            return False
        
        # Replace the original BAM file with the fixed one
        os.rename(fixed_bam, bam_file)
        
        # Re-index the BAM
        cmd = f"samtools index {bam_file}"
        logger.info(f"Running command: {cmd}")
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Failed to index the fixed BAM: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
            update_progress(analysis_id, "variants", "error", note="Failed to index BAM file with read groups")
            return False
        
        logger.info(f"Successfully added read groups to {bam_file}")
    
    return True

# Pipeline steps
PIPELINE_STEPS = [
    {"name": "FASTQ Upload", "id": "upload", "status": "complete", "description": "Uploading raw FASTQ files", "percent": 100},
    {"name": "Quality Check", "id": "qc", "status": "pending", "description": "Quality assessment of raw reads using FastQC", "percent": 0},
    {"name": "Trimming", "id": "trim", "status": "pending", "description": "Trimming of low-quality bases and adapters using fastp", "percent": 0},
    {"name": "Alignment", "id": "align", "status": "pending", "description": "Alignment to reference genome using BWA", "percent": 0},
    {"name": "Convert SAM → BAM", "id": "sam2bam", "status": "pending", "description": "Converting SAM to BAM format", "percent": 0},
    {"name": "Sort BAM", "id": "sort", "status": "pending", "description": "Sorting BAM file by coordinates", "percent": 0},
    {"name": "Mark Duplicates", "id": "dedup", "status": "pending", "description": "Marking duplicate reads using Picard", "percent": 0},
    {"name": "Index BAM", "id": "index", "status": "pending", "description": "Indexing BAM file for fast access", "percent": 0},
    {"name": "Base Recalibration", "id": "bqsr", "status": "pending", "description": "Base quality score recalibration using GATK", "percent": 0},
    {"name": "Variant Calling", "id": "variants", "status": "pending", "description": "Calling variants using GATK HaplotypeCaller", "percent": 0},
    {"name": "Variant Filtering", "id": "filter", "status": "pending", "description": "Filtering variants by quality and depth", "percent": 0},
    {"name": "Variant Annotation", "id": "annotate", "status": "pending", "description": "Annotating variants with functional information", "percent": 0},
    {"name": "VCF Compression", "id": "compress", "status": "pending", "description": "Compressing and indexing VCF file", "percent": 0},
    {"name": "Output VCF", "id": "vcf", "status": "pending", "description": "Generating final VCF file", "percent": 0}
]

def update_progress(analysis_id, step_id, status, output_file=None, percent=None, note=None):
    """Update the progress of a pipeline step."""
    result_folder = os.path.join(os.getcwd(), 'results', analysis_id)
    progress_file = os.path.join(result_folder, 'progress.json')
    
    # Load current progress
    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            progress = json.load(f)
    else:
        # Initialize progress with all steps
        progress = {"steps": PIPELINE_STEPS.copy(), "current_step": 0}
    
    # Update the specific step
    for step in progress["steps"]:
        if step["id"] == step_id:
            step["status"] = status
            if percent is not None:
                step["percent"] = percent
            if output_file:
                step["output_file"] = output_file
            if note:
                step["note"] = note
    
    # Save progress
    with open(progress_file, 'w') as f:
        json.dump(progress, f)
    
    # Also update the log
    log_file = os.path.join(result_folder, 'pipeline.log')
    percent_str = f" ({percent}%)" if percent is not None else ""
    note_str = f" - {note}" if note else ""
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Step {step_id}: {status}{percent_str}{note_str}\n")

def run_pipeline(analysis_id, r1_path=None, r2_path=None, bam_path=None, vcf_path=None, result_folder=None, skip_qc=False, skip_trim=False):
    """Run the exome analysis pipeline."""
    threading.Thread(target=_run_pipeline_thread, args=(analysis_id, r1_path, r2_path, bam_path, vcf_path, result_folder, skip_qc, skip_trim)).start()

def _run_pipeline_thread(analysis_id, r1_path=None, r2_path=None, bam_path=None, vcf_path=None, result_folder=None, skip_qc=False, skip_trim=False):
    """Run the pipeline in a separate thread."""
    try:
        logger.info(f"Starting pipeline for analysis {analysis_id}")
        
        # Initialize result folder
        if not result_folder:
            result_folder = os.path.join(os.getcwd(), 'results', analysis_id)
        os.makedirs(result_folder, exist_ok=True)
        
        # Initialize log file
        log_file = os.path.join(result_folder, 'pipeline.log')
        with open(log_file, 'w') as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Pipeline started for analysis {analysis_id}\n")
        
        # Initialize progress
        progress_file = os.path.join(result_folder, 'progress.json')
        
        # Create a copy of pipeline steps
        steps_to_run = PIPELINE_STEPS.copy()
        
        # Debug: log actual input paths to diagnose issues
        logger.info(f"Input paths - FASTQ R1: {r1_path}, FASTQ R2: {r2_path}, BAM: {bam_path}, VCF: {vcf_path}")
        
        # If starting with a BAM file, skip QC, trimming, and alignment
        if bam_path and os.path.exists(bam_path) and bam_path.endswith('.bam'):
            logger.info(f"Starting with BAM file: {bam_path}")
            # Mark QC, trim, and align steps as skipped
            for step in steps_to_run:
                if step["id"] in ["qc", "trim", "align"]:
                    step["status"] = "skipped"
                    step["percent"] = 100
                    step["note"] = "Skipped (starting with BAM file)"
        
        # If starting with a VCF file, skip all previous steps
        elif vcf_path and os.path.exists(vcf_path) and vcf_path.endswith('.vcf'):
            logger.info(f"Starting with VCF file: {vcf_path}")
            # Mark all previous steps as skipped
            for step in steps_to_run:
                if step["id"] in ["qc", "trim", "align", "dedup"]:
                    step["status"] = "skipped"
                    step["percent"] = 100
                    step["note"] = "Skipped (starting with VCF file)"
        
        # If skipping QC or trimming, mark them as skipped
        elif r1_path and r2_path:
            if skip_qc:
                for step in steps_to_run:
                    if step["id"] == "qc":
                        step["status"] = "skipped"
                        step["percent"] = 100
                        step["note"] = "Skipped as requested by user"
            
            if skip_trim:
                for step in steps_to_run:
                    if step["id"] == "trim":
                        step["status"] = "skipped"
                        step["percent"] = 100
                        step["note"] = "Skipped as requested by user"
        
        # Initialize progress with adjusted steps
        with open(progress_file, 'w') as f:
            json.dump({"steps": steps_to_run, "current_step": 0}, f)
        
        # If starting with FASTQ files
        if r1_path and r2_path:
            # Step 1: Quality Control with FastQC
            if skip_qc:
                logger.info("Skipping quality control step as requested by user")
            else:
                update_progress(analysis_id, "qc", "running", percent=0)
                logger.info("Running FastQC for quality control")
                
                qc_output_dir = os.path.join(result_folder, 'fastqc')
                os.makedirs(qc_output_dir, exist_ok=True)
                
                cmd = f"fastqc -o {qc_output_dir} {r1_path} {r2_path}"
                logger.info(f"Running command: {cmd}")
                
                # Verify FastQC is installed and FASTQ files exist
                if not os.path.exists(r1_path) or not os.path.exists(r2_path):
                    logger.error(f"FASTQ files not found: R1={r1_path}, R2={r2_path}")
                    update_progress(analysis_id, "qc", "error", note="FASTQ files not found")
                    return
                
                try:
                    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                            text=True, bufsize=1, universal_newlines=True)
                    
                    # Monitor output in real-time with progress updates
                    start_time = time.time()
                    last_update_time = start_time
                    current_percent = 10  # Start at 10% to show immediate feedback
                    update_progress(analysis_id, "qc", "running", percent=current_percent, 
                                  note="FastQC analysis started")
                    
                    while process.poll() is None:
                        stderr_line = process.stderr.readline()
                        if stderr_line:
                            stderr_line = stderr_line.strip()
                            logger.info(f"FastQC: {stderr_line}")
                            
                            # Look for progress indicators in FastQC output
                            if "Analysis complete" in stderr_line:
                                current_percent = 95
                                update_progress(analysis_id, "qc", "running", percent=current_percent,
                                             note="Analysis complete, generating reports")
                            elif "Started analysis" in stderr_line:
                                current_percent = 20
                                update_progress(analysis_id, "qc", "running", percent=current_percent,
                                             note="Analysis started")
                            elif "Sequence format detection completed" in stderr_line:
                                current_percent = 25
                                update_progress(analysis_id, "qc", "running", percent=current_percent)
                            elif "Approx" in stderr_line and "complete" in stderr_line:
                                try:
                                    # Try to extract percentage like "Approx 50% complete"
                                    percent_match = re.search(r"Approx (\d+)% complete", stderr_line)
                                    if percent_match:
                                        extracted_percent = int(percent_match.group(1))
                                        current_percent = 25 + (extracted_percent * 0.70)  # Scale to leave room for report generation
                                        update_progress(analysis_id, "qc", "running", percent=int(current_percent))
                                except:
                                    pass
                        
                        stdout_line = process.stdout.readline()
                        if stdout_line:
                            stdout_line = stdout_line.strip()
                            logger.info(f"FastQC: {stdout_line}")
                        
                        # Time-based fallback updates
                        current_time = time.time()
                        if current_time - last_update_time > 5:  # Update every 5 seconds
                            elapsed_seconds = current_time - start_time
                            # Assume FastQC takes about 2 minutes max
                            time_based_percent = min(95, 10 + int((elapsed_seconds / 120) * 85))
                            
                            if time_based_percent > current_percent:
                                current_percent = time_based_percent
                                update_progress(analysis_id, "qc", "running", percent=current_percent, 
                                             note=f"Processing ({int(elapsed_seconds)} seconds elapsed)")
                            
                            last_update_time = current_time
                        
                        time.sleep(0.1)
                    
                    # Get remaining output
                    stdout, stderr = process.communicate()
                    
                    if process.returncode != 0:
                        logger.error(f"FastQC failed: {stderr}")
                        update_progress(analysis_id, "qc", "error", note=f"FastQC failed: {stderr}")
                        return
                except Exception as e:
                    logger.error(f"FastQC execution error: {str(e)}")
                    update_progress(analysis_id, "qc", "error", note=f"FastQC execution error: {str(e)}")
                    return
                
                # Rename output files to simpler names
                try:
                    # FastQC output filename format might vary, attempt to find the generated files
                    r1_basename = os.path.basename(r1_path)
                    r2_basename = os.path.basename(r2_path)
                    
                    # Handle various extensions
                    for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
                        if r1_basename.endswith(ext):
                            r1_fastqc_html = os.path.join(qc_output_dir, r1_basename.replace(ext, '_fastqc.html'))
                            if os.path.exists(r1_fastqc_html):
                                break
                    
                    for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
                        if r2_basename.endswith(ext):
                            r2_fastqc_html = os.path.join(qc_output_dir, r2_basename.replace(ext, '_fastqc.html'))
                            if os.path.exists(r2_fastqc_html):
                                break
                    
                    # Check if files exist before copying
                    if os.path.exists(r1_fastqc_html):
                        shutil.copy(r1_fastqc_html, os.path.join(qc_output_dir, 'r1_fastqc.html'))
                        logger.info(f"Copied {r1_fastqc_html} to r1_fastqc.html")
                    else:
                        logger.warning(f"FastQC HTML not found at expected path: {r1_fastqc_html}")
                        # Look for any HTML files in the output directory
                        html_files = [f for f in os.listdir(qc_output_dir) if f.endswith('.html')]
                        if html_files:
                            logger.info(f"Found alternative FastQC HTML files: {html_files}")
                            r1_fastqc_html = os.path.join(qc_output_dir, html_files[0])
                    
                    if os.path.exists(r2_fastqc_html):
                        shutil.copy(r2_fastqc_html, os.path.join(qc_output_dir, 'r2_fastqc.html'))
                        logger.info(f"Copied {r2_fastqc_html} to r2_fastqc.html")
                    
                    logger.info("FastQC completed successfully")
                    update_progress(analysis_id, "qc", "complete", r1_fastqc_html, percent=100)
                except Exception as e:
                    logger.warning(f"Error processing FastQC output files: {str(e)}")
                    # Still mark as complete if FastQC ran successfully
                    update_progress(analysis_id, "qc", "complete", "fastqc/r1_fastqc.html", percent=100, 
                                  note="FastQC complete but output processing had errors")

            # Step 2: Trimming with fastp
            r1_processed = r1_path
            r2_processed = r2_path
            
            if skip_trim:
                logger.info("Skipping trimming step as requested by user")
                
                # Create directory structure for consistency even when skipping
                trim_output_dir = os.path.join(result_folder, 'trimmed')
                os.makedirs(trim_output_dir, exist_ok=True)
            else:
                update_progress(analysis_id, "trim", "running", percent=0)
                logger.info("Trimming reads with fastp")
                
                trim_output_dir = os.path.join(result_folder, 'trimmed')
                os.makedirs(trim_output_dir, exist_ok=True)
                
                r1_trimmed = os.path.join(trim_output_dir, 'r1_trimmed.fastq.gz')
                r2_trimmed = os.path.join(trim_output_dir, 'r2_trimmed.fastq.gz')
                trim_report = os.path.join(trim_output_dir, 'fastp_report.html')
                
                cmd = f"fastp -i {r1_path} -I {r2_path} -o {r1_trimmed} -O {r2_trimmed} -h {trim_report}"
                logger.info(f"Running command: {cmd}")
                
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                         text=True, bufsize=1, universal_newlines=True)
                
                # Get file size in MB for progress estimation
                r1_size_mb = os.path.getsize(r1_path) / (1024 * 1024)
                
                # Monitor output in real-time with percentage updates
                start_time = time.time()
                last_update = start_time
                current_percent = 0
                
                # Regex to match fastp progress like "Read: 1000000 reads" and extract the number
                import re
                fastp_progress_pattern = re.compile(r"Read: (\d+) reads")
                total_reads_estimate = r1_size_mb * 250000  # Estimate: ~250,000 reads per MB
                
                while process.poll() is None:
                    stderr_line = process.stderr.readline()
                    if stderr_line:
                        stderr_line = stderr_line.strip()
                        logger.info(f"fastp: {stderr_line}")
                        
                        # Extract progress information if available
                        match = fastp_progress_pattern.search(stderr_line)
                        if match:
                            reads_processed = int(match.group(1))
                            percent = min(99, int((reads_processed / total_reads_estimate) * 100))
                            current_percent = percent
                            update_progress(analysis_id, "trim", "running", percent=current_percent, 
                                          note=f"Processing reads: {reads_processed:,} reads processed")
                    
                    # Check for any output from stdout
                    stdout_line = process.stdout.readline()
                    if stdout_line:
                        stdout_line = stdout_line.strip()
                        logger.info(f"fastp stdout: {stdout_line}")
                    
                    # Time-based fallback - show some progress even if we can't parse it from output
                    elapsed_seconds = time.time() - start_time
                    time_based_percent = min(99, int((elapsed_seconds / 300) * 100))  # Assume ~5 min for trimming
                    
                    if time_based_percent > current_percent:
                        current_percent = time_based_percent
                        update_progress(analysis_id, "trim", "running", percent=current_percent, 
                                      note=f"Trimming in progress (time-based estimate)")
                    
                    # Check if we should provide a progress update in the logs
                    current_time = time.time()
                    if current_time - last_update > 30:  # Update every 30 seconds
                        logger.info(f"fastp is still running... ~{current_percent}% complete (elapsed: {int(elapsed_seconds)} seconds)")
                        last_update = current_time
                    
                    # Small sleep to prevent CPU hogging
                    time.sleep(0.1)
                
                # Get remaining output
                stdout, stderr = process.communicate()
                if stdout:
                    logger.info(f"fastp stdout: {stdout}")
                if stderr:
                    logger.info(f"fastp stderr: {stderr}")
                
                if process.returncode != 0:
                    logger.error(f"Fastp failed with return code {process.returncode}")
                    update_progress(analysis_id, "trim", "error")
                    return
                
                logger.info("Trimming completed successfully")
                update_progress(analysis_id, "trim", "complete", trim_report, percent=100)
                
                # Use trimmed files for subsequent steps
                r1_processed = r1_trimmed
                r2_processed = r2_trimmed
            
            # Step 3: Alignment with BWA
            update_progress(analysis_id, "align", "running", percent=0)
            logger.info("Aligning reads to reference genome")
            logger.info(f"Using processed reads: {r1_processed} and {r2_processed}")
            
            align_output_dir = os.path.join(result_folder, 'aligned')
            os.makedirs(align_output_dir, exist_ok=True)
            
            # Ensure reference genome indices are created before alignment
            if not ensure_reference_indices(result_folder, analysis_id):
                return
            
            sample_name = f"sample_{analysis_id}"
            sam_output = os.path.join(align_output_dir, f"{sample_name}.sam")
            bam_output = os.path.join(align_output_dir, f"{sample_name}.bam")
            sorted_bam = os.path.join(align_output_dir, f"{sample_name}.sorted.bam")
            
            # BWA alignment command with read group information
            # Add proper read group with sample ID for GATK compatibility
            sample_id = f"{analysis_id}"
            platform = "ILLUMINA"
            library = f"lib_{sample_id}"
            
            # Format the read group string according to SAM specification
            read_group = f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:{platform}\\tLB:{library}\\tPU:unit1"
            
            # Updated BWA command with read group
            cmd_bwa = f"bwa mem -t 4 -R '{read_group}' {REFERENCE_GENOME} {r1_processed} {r2_processed} > {sam_output}"
            logger.info(f"Running command: {cmd_bwa}")
            
            # For BWA, we need to use shell=True and capture in a different way
            process = subprocess.Popen(cmd_bwa, shell=True, stderr=subprocess.PIPE, 
                                     text=True, bufsize=1, universal_newlines=True)
            
            # Monitor output in real-time with percentage updates
            start_time = time.time()
            last_update = start_time
            current_percent = 0
            
            # BWA progress patterns to look for in output
            bwa_progress_patterns = {
                "Started analysis": 5,
                "Approx 5% complete": 5,
                "Approx 10% complete": 10,
                "Approx 15% complete": 15,
                "Approx 20% complete": 20,
                "Approx 25% complete": 25,
                "Approx 30% complete": 30,
                "Approx 35% complete": 35,
                "Approx 40% complete": 40,
                "Approx 45% complete": 45,
                "Approx 50% complete": 50,
                "Approx 55% complete": 55,
                "Approx 60% complete": 60,
                "Approx 65% complete": 65,
                "Approx 70% complete": 70,
                "Approx 75% complete": 75,
                "Approx 80% complete": 80,
                "Approx 85% complete": 85,
                "Approx 90% complete": 90,
                "Approx 95% complete": 95,
                "Analysis complete": 100
            }
            
            while process.poll() is None:
                stderr_line = process.stderr.readline()
                if stderr_line:
                    stderr_line = stderr_line.strip()
                    logger.info(f"BWA: {stderr_line}")
                    
                    # Check for progress indicators in BWA output
                    for pattern, percent in bwa_progress_patterns.items():
                        if pattern in stderr_line:
                            current_percent = percent
                            update_progress(analysis_id, "align", "running", percent=current_percent)
                            break
                
                # Calculate time-based percentage estimate if no direct progress indicators
                elapsed_seconds = time.time() - start_time
                time_based_percent = min(99, int((elapsed_seconds / 300) * 100))  # Assume ~5 min for alignment
                
                if time_based_percent > current_percent:
                    current_percent = time_based_percent
                    update_progress(analysis_id, "align", "running", percent=current_percent)
                
                # Check if we should provide a progress update
                current_time = time.time()
                if current_time - last_update > 30:  # Update every 30 seconds
                    logger.info(f"BWA alignment is still running... ~{current_percent}% complete (elapsed: {int(elapsed_seconds)} seconds)")
                    last_update = current_time
                
                # Small sleep to prevent CPU hogging
                time.sleep(0.1)
            
            # Get remaining output
            _, stderr = process.communicate()
            if stderr:
                logger.info(f"BWA stderr: {stderr}")
            
            if process.returncode != 0:
                logger.error(f"BWA alignment failed with return code {process.returncode}")
                update_progress(analysis_id, "align", "error")
                return
            
            logger.info("Alignment completed successfully")
            update_progress(analysis_id, "align", "complete", sam_output, percent=100)
            
            # Step 4: Convert SAM to BAM
            update_progress(analysis_id, "sam2bam", "running", percent=0, 
                           note="Converting SAM to BAM format")
            
            cmd_samtools1 = f"samtools view -S -b {sam_output} > {bam_output}"
            logger.info(f"Running command: {cmd_samtools1}")
            update_progress(analysis_id, "sam2bam", "running", percent=50)
            
            process = subprocess.Popen(cmd_samtools1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                logger.error(f"Samtools view failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                update_progress(analysis_id, "sam2bam", "error")
                return
            
            logger.info("SAM to BAM conversion completed successfully")
            update_progress(analysis_id, "sam2bam", "complete", bam_output, percent=100)
            
            # Step 5: Sort BAM
            update_progress(analysis_id, "sort", "running", percent=0,
                           note="Sorting BAM file by coordinates")
            
            cmd_samtools2 = f"samtools sort {bam_output} -o {sorted_bam}"
            logger.info(f"Running command: {cmd_samtools2}")
            update_progress(analysis_id, "sort", "running", percent=50)
            
            process = subprocess.Popen(cmd_samtools2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                logger.error(f"Samtools sort failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                update_progress(analysis_id, "sort", "error")
                return
            
            logger.info("BAM sorting completed successfully")
            update_progress(analysis_id, "sort", "complete", sorted_bam, percent=100)
            
            # Step 6: Index BAM
            update_progress(analysis_id, "index", "running", percent=0,
                           note="Indexing sorted BAM file")
            
            cmd_samtools3 = f"samtools index {sorted_bam}"
            logger.info(f"Running command: {cmd_samtools3}")
            update_progress(analysis_id, "index", "running", percent=50)
            
            process = subprocess.Popen(cmd_samtools3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                update_progress(analysis_id, "index", "error")
                return
            
            logger.info("BAM indexing completed successfully")
            update_progress(analysis_id, "index", "complete", f"{sorted_bam}.bai", percent=100)
        
        # If we're starting with a BAM file, prepare the paths
        elif bam_path:
            logger.info(f"Using provided BAM file: {bam_path}")
            
            # Copy the BAM file to the alignments directory with the expected name
            align_output_dir = os.path.join(result_folder, 'aligned')
            os.makedirs(align_output_dir, exist_ok=True)
            
            sample_name = f"sample_{analysis_id}"
            sorted_bam = os.path.join(align_output_dir, f"{sample_name}.sorted.bam")
            
            # Copy or link the BAM file
            shutil.copy2(bam_path, sorted_bam)
            
            # Index BAM if needed
            update_progress(analysis_id, "index", "running", percent=0,
                           note="Indexing provided BAM file")
            
            if not os.path.exists(f"{sorted_bam}.bai"):
                logger.info("Indexing BAM file")
                cmd_samtools3 = f"samtools index {sorted_bam}"
                logger.info(f"Running command: {cmd_samtools3}")
                
                process = subprocess.Popen(cmd_samtools3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                
                if process.returncode != 0:
                    logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                    update_progress(analysis_id, "index", "error")
                    return
            
            update_progress(analysis_id, "index", "complete", f"{sorted_bam}.bai", percent=100)
            logger.info("BAM file prepared successfully")
        
        # Continue only if not starting with a VCF file
        if not vcf_path:
            # Step 7: Mark Duplicates
            update_progress(analysis_id, "dedup", "running", percent=0,
                           note="Marking duplicate reads with Picard")
            logger.info("Marking duplicate reads with Picard")
            
            # Create dedup directory
            dedup_output_dir = os.path.join(result_folder, 'dedup')
            os.makedirs(dedup_output_dir, exist_ok=True)
            
            sample_name = f"sample_{analysis_id}"
            dedup_bam = os.path.join(dedup_output_dir, f"{sample_name}.dedup.bam")
            metrics_file = os.path.join(dedup_output_dir, f"{sample_name}.dedup.metrics.txt")
            
            # Path to sorted BAM (either from alignment or directly provided)
            sorted_bam = os.path.join(result_folder, 'aligned', f"{sample_name}.sorted.bam")
            
            # Use standalone Picard JAR with native syntax
            picard_jar = os.path.join(os.getcwd(), 'tools', 'picard.jar')
            cmd_dedup = f"java -jar {picard_jar} MarkDuplicates \
                        I={sorted_bam} \
                        O={dedup_bam} \
                        M={metrics_file} \
                        VALIDATION_STRINGENCY=LENIENT \
                        CREATE_INDEX=true"
            
            logger.info(f"Running command: {cmd_dedup}")
            
            process = subprocess.Popen(cmd_dedup, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     text=True, bufsize=1, universal_newlines=True)
            
            # Monitor output in real-time and estimate progress
            start_time = time.time()
            last_update_time = start_time
            
            while process.poll() is None:
                stderr_line = process.stderr.readline()
                if stderr_line:
                    logger.info(f"Picard MarkDuplicates: {stderr_line.strip()}")
                    
                    # Update progress periodically based on log messages
                    current_time = time.time()
                    elapsed_seconds = current_time - start_time
                    
                    # Update progress every 30 seconds to avoid too many updates
                    if current_time - last_update_time > 30:
                        # Estimate progress based on elapsed time (assuming ~10 min total)
                        time_based_percent = min(95, int((elapsed_seconds / 600) * 100))
                        current_percent = max(25, time_based_percent)
                        
                        update_progress(analysis_id, "dedup", "running", percent=current_percent,
                                      note=f"Marking duplicates (running for {int(elapsed_seconds)} seconds)")
                        
                        last_update_time = current_time
                
                stdout_line = process.stdout.readline()
                if stdout_line:
                    logger.info(f"Picard MarkDuplicates: {stdout_line.strip()}")
                
                time.sleep(0.1)
            
            # Get remaining output
            stdout, stderr = process.communicate()
            if stdout:
                logger.info(f"Picard MarkDuplicates stdout: {stdout}")
            if stderr:
                logger.info(f"Picard MarkDuplicates stderr: {stderr}")
            
            if process.returncode != 0:
                logger.error(f"Picard MarkDuplicates failed with return code {process.returncode}")
                update_progress(analysis_id, "dedup", "error")
                return
            
            logger.info("Duplicate marking completed successfully")
            update_progress(analysis_id, "dedup", "complete", dedup_bam, percent=100)
            
            # Step 8: Index the BAM file
            update_progress(analysis_id, "index", "running", percent=0,
                           note="Indexing BAM file")
            
            # Check if index was already created by Picard
            bai_file = f"{dedup_bam}.bai"
            if not os.path.exists(bai_file):
                cmd_index = f"samtools index {dedup_bam}"
                logger.info(f"Running command: {cmd_index}")
                
                process = subprocess.Popen(cmd_index, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                
                if process.returncode != 0:
                    logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                    update_progress(analysis_id, "index", "error")
                    return
            else:
                logger.info("BAM index already exists from Picard")
            
            logger.info("BAM indexing completed successfully")
            update_progress(analysis_id, "index", "complete", f"{dedup_bam}.bai", percent=100)
            
            # Step 9: Base Quality Score Recalibration (BQSR)
            update_progress(analysis_id, "bqsr", "running", percent=0,
                           note="Running Base Quality Score Recalibration")
            
            bqsr_output_dir = os.path.join(result_folder, 'bqsr')
            os.makedirs(bqsr_output_dir, exist_ok=True)
            
            # First, create recalibration table
            recal_table = os.path.join(bqsr_output_dir, f"{sample_name}.recal.table")
            bqsr_bam = os.path.join(bqsr_output_dir, f"{sample_name}.bqsr.bam")
            
            # Path to known variant sites (dbSNP)
            dbsnp_vcf = os.path.join(os.getcwd(), 'reference', 'dbsnp_138.hg19.vcf')
            
            # Check if dbSNP file exists
            if not os.path.exists(dbsnp_vcf):
                logger.warning(f"dbSNP file not found at {dbsnp_vcf}. Skipping BQSR.")
                update_progress(analysis_id, "bqsr", "skipped", percent=100,
                              note="Skipped - dbSNP file not available")
            else:
                # GATK BaseRecalibrator command to create recalibration table
                gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
                
                cmd_bqsr1 = f"{gatk_path} BaseRecalibrator \
                            -R {REFERENCE_GENOME} \
                            -I {dedup_bam} \
                            --known-sites {dbsnp_vcf} \
                            -O {recal_table}"
                
                logger.info(f"Running command: {cmd_bqsr1}")
                update_progress(analysis_id, "bqsr", "running", percent=25, 
                              note="Generating recalibration table")
                
                # Start the process
                process = subprocess.Popen(cmd_bqsr1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         text=True, bufsize=1, universal_newlines=True)
                
                # Monitor output in real-time
                start_time = time.time()
                last_update_time = start_time
                
                while process.poll() is None:
                    stderr_line = process.stderr.readline()
                    if stderr_line:
                        logger.info(f"GATK BaseRecalibrator: {stderr_line.strip()}")
                        
                        # Update progress periodically
                        current_time = time.time()
                        if current_time - last_update_time > 30:
                            elapsed_seconds = current_time - start_time
                            # Estimate progress based on elapsed time (assuming ~5 min for this step)
                            time_based_percent = min(45, 25 + int((elapsed_seconds / 300) * 20))
                            
                            update_progress(analysis_id, "bqsr", "running", percent=time_based_percent)
                            last_update_time = current_time
                    
                    stdout_line = process.stdout.readline()
                    if stdout_line:
                        logger.info(f"GATK BaseRecalibrator: {stdout_line.strip()}")
                    
                    time.sleep(0.1)
                
                # Get remaining output
                stdout, stderr = process.communicate()
                if stdout:
                    logger.info(f"GATK BaseRecalibrator stdout: {stdout}")
                if stderr:
                    logger.info(f"GATK BaseRecalibrator stderr: {stderr}")
                
                if process.returncode != 0:
                    logger.error(f"GATK BaseRecalibrator failed with return code {process.returncode}")
                    update_progress(analysis_id, "bqsr", "error")
                    return
                
                logger.info("Recalibration table generated successfully")
                update_progress(analysis_id, "bqsr", "running", percent=50, 
                              note="Applying base recalibration")
                
                # Apply recalibration with ApplyBQSR
                cmd_bqsr2 = f"{gatk_path} ApplyBQSR \
                            -R {REFERENCE_GENOME} \
                            -I {dedup_bam} \
                            --bqsr-recal-file {recal_table} \
                            -O {bqsr_bam}"
                
                logger.info(f"Running command: {cmd_bqsr2}")
                
                # Start the process
                process = subprocess.Popen(cmd_bqsr2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         text=True, bufsize=1, universal_newlines=True)
                
                # Monitor output in real-time
                start_time = time.time()
                last_update_time = start_time
                
                while process.poll() is None:
                    stderr_line = process.stderr.readline()
                    if stderr_line:
                        logger.info(f"GATK ApplyBQSR: {stderr_line.strip()}")
                        
                        # Update progress periodically
                        current_time = time.time()
                        if current_time - last_update_time > 30:
                            elapsed_seconds = current_time - start_time
                            # Estimate progress based on elapsed time (assuming ~5 min for this step)
                            time_based_percent = min(95, 50 + int((elapsed_seconds / 300) * 45))
                            
                            update_progress(analysis_id, "bqsr", "running", percent=time_based_percent)
                            last_update_time = current_time
                    
                    stdout_line = process.stdout.readline()
                    if stdout_line:
                        logger.info(f"GATK ApplyBQSR: {stdout_line.strip()}")
                    
                    time.sleep(0.1)
                
                # Get remaining output
                stdout, stderr = process.communicate()
                if stdout:
                    logger.info(f"GATK ApplyBQSR stdout: {stdout}")
                if stderr:
                    logger.info(f"GATK ApplyBQSR stderr: {stderr}")
                
                if process.returncode != 0:
                    logger.error(f"GATK ApplyBQSR failed with return code {process.returncode}")
                    update_progress(analysis_id, "bqsr", "error")
                    return
                
                # Create index for BQSR BAM
                cmd_index = f"samtools index {bqsr_bam}"
                logger.info(f"Running command: {cmd_index}")
                
                process = subprocess.Popen(cmd_index, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                
                if process.returncode != 0:
                    logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                    update_progress(analysis_id, "bqsr", "error", note="Failed to index BQSR BAM file")
                    return
                
                logger.info("Base quality score recalibration completed successfully")
                update_progress(analysis_id, "bqsr", "complete", bqsr_bam, percent=100)
                
                # Use BQSR BAM for variant calling
                dedup_bam = bqsr_bam
            
            # Step 10: Variant Calling with GATK
            update_progress(analysis_id, "variants", "running", percent=0,
                           note="Calling variants with GATK HaplotypeCaller")
            
            variants_output_dir = os.path.join(result_folder, 'variants')
            os.makedirs(variants_output_dir, exist_ok=True)
            
            vcf_output = os.path.join(variants_output_dir, 'variants.vcf')
            
            # Ensure reference genome has all required indices
            if not ensure_reference_indices(result_folder, analysis_id):
                return
            
            # Verify BAM file compatibility with GATK
            if not verify_bam_for_gatk(dedup_bam, analysis_id, result_folder):
                logger.error(f"BAM file {dedup_bam} is not compatible with GATK and could not be fixed")
                update_progress(analysis_id, "variants", "error", 
                               note="BAM file is not compatible with GATK")
                return
            
            # Using HaplotypeCaller with deduplicated BAM
            gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
            cmd_gatk = f"{gatk_path} HaplotypeCaller -R {REFERENCE_GENOME} -I {dedup_bam} -O {vcf_output}"
            logger.info(f"Running command: {cmd_gatk}")
            
            # GATK prints progress to stderr
            process = subprocess.Popen(cmd_gatk, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                     text=True, bufsize=1, universal_newlines=True)
            
            # Monitor output in real-time with percentage updates
            start_time = time.time()
            last_update = start_time
            current_percent = 0
            
            # GATK progress patterns to look for in output
            gatk_progress_patterns = {
                "Started analysis": 5,
                "Approx 5% complete": 5,
                "Approx 10% complete": 10,
                "Approx 15% complete": 15,
                "Approx 20% complete": 20,
                "Approx 25% complete": 25,
                "Approx 30% complete": 30,
                "Approx 35% complete": 35,
                "Approx 40% complete": 40,
                "Approx 45% complete": 45,
                "Approx 50% complete": 50,
                "Approx 55% complete": 55,
                "Approx 60% complete": 60,
                "Approx 65% complete": 65,
                "Approx 70% complete": 70,
                "Approx 75% complete": 75,
                "Approx 80% complete": 80,
                "Approx 85% complete": 85,
                "Approx 90% complete": 90,
                "Approx 95% complete": 95,
                "Analysis complete": 100
            }
            
            while process.poll() is None:
                stderr_line = process.stderr.readline()
                if stderr_line:
                    stderr_line = stderr_line.strip()
                    logger.info(f"GATK: {stderr_line}")
                    
                    # Check for progress indicators in GATK output
                    for pattern, percent in gatk_progress_patterns.items():
                        if pattern in stderr_line:
                            current_percent = percent
                            update_progress(analysis_id, "variants", "running", percent=current_percent)
                            break
                
                # Calculate time-based percentage estimate if no direct progress indicators
                elapsed_seconds = time.time() - start_time
                time_based_percent = min(99, int((elapsed_seconds / 300) * 100))  # Assume ~5 min for variant calling
                
                if time_based_percent > current_percent:
                    current_percent = time_based_percent
                    update_progress(analysis_id, "variants", "running", percent=current_percent)
                
                # Check if we should provide a progress update
                current_time = time.time()
                if current_time - last_update > 30:  # Update every 30 seconds
                    logger.info(f"GATK variant calling is still running... ~{current_percent}% complete (elapsed: {int(elapsed_seconds)} seconds)")
                    last_update = current_time
                
                # Small sleep to prevent CPU hogging
                time.sleep(0.1)
            
            # Get remaining output
            stdout, stderr = process.communicate()
            if stdout:
                logger.info(f"GATK stdout: {stdout}")
            if stderr:
                logger.info(f"GATK stderr: {stderr}")
            
            if process.returncode != 0:
                logger.error(f"GATK variant calling failed with return code {process.returncode}")
                update_progress(analysis_id, "variants", "error")
                return
            
            # Copy final VCF to results root for easy access
            final_vcf = os.path.join(result_folder, 'variants.vcf')
            shutil.copy2(vcf_output, final_vcf)
            
            logger.info("Variant calling completed successfully")
            update_progress(analysis_id, "variants", "complete", vcf_output, percent=100)
            
            # Step 11: Variant Filtering
            update_progress(analysis_id, "filter", "running", percent=0,
                          note="Filtering variants by quality and depth")
            
            filter_output_dir = os.path.join(result_folder, 'filtered')
            os.makedirs(filter_output_dir, exist_ok=True)
            
            filtered_vcf = os.path.join(filter_output_dir, 'filtered_variants.vcf')
            
            # GATK VariantFiltration command
            cmd_filter = f"{gatk_path} VariantFiltration \
                         -R {REFERENCE_GENOME} \
                         -V {vcf_output} \
                         -O {filtered_vcf} \
                         --filter-name 'QualFilter' --filter-expression 'QUAL < 30.0' \
                         --filter-name 'DepthFilter' --filter-expression 'DP < 10' \
                         --filter-name 'QDFilter' --filter-expression 'QD < 2.0'"
            
            logger.info(f"Running command: {cmd_filter}")
            update_progress(analysis_id, "filter", "running", percent=30)
            
            process = subprocess.Popen(cmd_filter, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     text=True, bufsize=1, universal_newlines=True)
            
            # Monitor output
            start_time = time.time()
            last_update_time = start_time
            
            while process.poll() is None:
                stderr_line = process.stderr.readline()
                if stderr_line:
                    logger.info(f"GATK VariantFiltration: {stderr_line.strip()}")
                    
                    # Update progress periodically
                    current_time = time.time()
                    if current_time - last_update_time > 15:
                        elapsed_seconds = current_time - start_time
                        # Assume filtering takes about 1 minute total
                        time_based_percent = min(95, 30 + int((elapsed_seconds / 60) * 65))
                        
                        update_progress(analysis_id, "filter", "running", percent=time_based_percent)
                        last_update_time = current_time
            
                stdout_line = process.stdout.readline()
                if stdout_line:
                    logger.info(f"GATK VariantFiltration: {stdout_line.strip()}")
            
                time.sleep(0.1)
            
            # Get remaining output
            stdout, stderr = process.communicate()
            if stdout:
                logger.info(f"GATK VariantFiltration stdout: {stdout}")
            if stderr:
                logger.info(f"GATK VariantFiltration stderr: {stderr}")
            
            if process.returncode != 0:
                logger.error(f"GATK VariantFiltration failed with return code {process.returncode}")
                update_progress(analysis_id, "filter", "error")
                # Continue to next step even if this fails
            else:
                logger.info("Variant filtering completed successfully")
                update_progress(analysis_id, "filter", "complete", filtered_vcf, percent=100)
                
                # Use filtered VCF for subsequent steps
                vcf_output = filtered_vcf
            
            # Step 12: Variant Annotation
            update_progress(analysis_id, "annotate", "running", percent=0,
                          note="Annotating variants with functional information")
            
            # Check if a suitable annotation tool is available
            # For this example, we're using GATK Funcotator, but in practice might use tools like SnpEff, VEP, etc.
            annotate_output_dir = os.path.join(result_folder, 'annotated')
            os.makedirs(annotate_output_dir, exist_ok=True)
            
            annotated_vcf = os.path.join(annotate_output_dir, 'annotated_variants.vcf')
            
            # Check for data sources directory for annotation
            data_sources_dir = os.path.join(os.getcwd(), 'reference', 'funcotator_dataSources')
            
            if not os.path.exists(data_sources_dir):
                logger.warning(f"Annotation data sources not found at {data_sources_dir}. Skipping annotation.")
                update_progress(analysis_id, "annotate", "skipped", percent=100,
                              note="Skipped - annotation data sources not available")
            else:
                # GATK Funcotator command
                cmd_annotate = f"{gatk_path} Funcotator \
                               --variant {vcf_output} \
                               --reference {REFERENCE_GENOME} \
                               --data-sources-path {data_sources_dir} \
                               --output {annotated_vcf} \
                               --output-file-format VCF"
                
                logger.info(f"Running command: {cmd_annotate}")
                update_progress(analysis_id, "annotate", "running", percent=30,
                              note="Annotating variants (this may take 5-10 minutes)")
                
                process = subprocess.Popen(cmd_annotate, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         text=True, bufsize=1, universal_newlines=True)
                
                # Monitor output
                start_time = time.time()
                last_update_time = start_time
                
                while process.poll() is None:
                    stderr_line = process.stderr.readline()
                    if stderr_line:
                        logger.info(f"GATK Funcotator: {stderr_line.strip()}")
                        
                        # Update progress periodically
                        current_time = time.time()
                        if current_time - last_update_time > 30:
                            elapsed_seconds = current_time - start_time
                            # Assume annotation takes about 5-10 minutes
                            time_based_percent = min(95, 30 + int((elapsed_seconds / 600) * 65))
                            
                            update_progress(analysis_id, "annotate", "running", percent=time_based_percent)
                            last_update_time = current_time
                    
                    stdout_line = process.stdout.readline()
                    if stdout_line:
                        logger.info(f"GATK Funcotator: {stdout_line.strip()}")
                    
                    time.sleep(0.1)
                
                # Get remaining output
                stdout, stderr = process.communicate()
                if stdout:
                    logger.info(f"GATK Funcotator stdout: {stdout}")
                if stderr:
                    logger.info(f"GATK Funcotator stderr: {stderr}")
                
                if process.returncode != 0:
                    logger.error(f"GATK Funcotator failed with return code {process.returncode}")
                    update_progress(analysis_id, "annotate", "error")
                    # Continue to next step even if this fails
                else:
                    logger.info("Variant annotation completed successfully")
                    update_progress(analysis_id, "annotate", "complete", annotated_vcf, percent=100)
                    
                    # Use annotated VCF for subsequent steps
                    vcf_output = annotated_vcf
            
            # Step 13: VCF Compression and Indexing
            update_progress(analysis_id, "compress", "running", percent=0,
                          note="Compressing and indexing VCF file")
            
            compress_output_dir = os.path.join(result_folder, 'compressed')
            os.makedirs(compress_output_dir, exist_ok=True)
            
            compressed_vcf = os.path.join(compress_output_dir, 'variants.vcf.gz')
            
            # Compress VCF with bgzip
            cmd_bgzip = f"bgzip -c {vcf_output} > {compressed_vcf}"
            logger.info(f"Running command: {cmd_bgzip}")
            update_progress(analysis_id, "compress", "running", percent=30)
            
            process = subprocess.Popen(cmd_bgzip, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                logger.error(f"bgzip compression failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                update_progress(analysis_id, "compress", "error")
                # Continue to final step even if this fails
            else:
                # Index the compressed VCF with tabix
                cmd_tabix = f"tabix -p vcf {compressed_vcf}"
                logger.info(f"Running command: {cmd_tabix}")
                update_progress(analysis_id, "compress", "running", percent=70)
                
                process = subprocess.Popen(cmd_tabix, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                
                if process.returncode != 0:
                    logger.error(f"tabix indexing failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                    update_progress(analysis_id, "compress", "error")
                else:
                    logger.info("VCF compression and indexing completed successfully")
                    update_progress(analysis_id, "compress", "complete", compressed_vcf, percent=100)
                    
                    # Use compressed VCF for final step
                    vcf_output = compressed_vcf
            
            # Step 14: Output Final VCF
            update_progress(analysis_id, "vcf", "running", percent=50,
                          note="Finalizing VCF file")
            
            # Copy final VCF to the root results directory for easy access
            final_vcf = os.path.join(result_folder, 'variants.vcf')
            
            # Check if the final output is compressed
            if vcf_output.endswith('.gz'):
                # For compressed VCF, we'll create both compressed and uncompressed copies
                try:
                    # Create uncompressed copy for easier viewing
                    with gzip.open(vcf_output, 'rb') as f_in:
                        with open(final_vcf, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    
                    # Also copy the compressed version and its index
                    shutil.copy2(vcf_output, os.path.join(result_folder, 'variants.vcf.gz'))
                    if os.path.exists(f"{vcf_output}.tbi"):
                        shutil.copy2(f"{vcf_output}.tbi", os.path.join(result_folder, 'variants.vcf.gz.tbi'))
                except Exception as e:
                    logger.error(f"Error creating uncompressed VCF copy: {str(e)}")
                    # Just copy the compressed file as-is
                    shutil.copy2(vcf_output, os.path.join(result_folder, 'variants.vcf.gz'))
            else:
                # For uncompressed VCF, just copy it
                shutil.copy2(vcf_output, final_vcf)
            
            # Create a summary file with variant counts
            try:
                if os.path.exists(final_vcf):
                    # Count variants using grep
                    cmd_count = f"grep -v '^#' {final_vcf} | wc -l"
                    process = subprocess.Popen(cmd_count, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    stdout, stderr = process.communicate()
                    
                    if process.returncode == 0:
                        variant_count = stdout.decode().strip() if isinstance(stdout, bytes) else stdout.strip()
                        
                        summary_file = os.path.join(result_folder, 'variant_summary.txt')
                        with open(summary_file, 'w') as f:
                            f.write(f"Analysis ID: {analysis_id}\n")
                            f.write(f"Total variants: {variant_count}\n")
                            f.write(f"Date completed: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            except Exception as e:
                logger.error(f"Error creating variant summary: {str(e)}")
            
            logger.info("Pipeline complete. VCF file ready.")
            update_progress(analysis_id, "vcf", "complete", final_vcf, percent=100)
        
        # If starting with a VCF file, just copy it to the results directory
        elif vcf_path:
            logger.info(f"Using provided VCF file: {vcf_path}")
            
            # Create variants directory and copy VCF file
            variants_output_dir = os.path.join(result_folder, 'variants')
            os.makedirs(variants_output_dir, exist_ok=True)
            
            vcf_output = os.path.join(variants_output_dir, 'variants.vcf')
            
            # Copy the VCF file
            shutil.copy2(vcf_path, vcf_output)
            
            # Also copy to the root results directory for easy access
            final_vcf = os.path.join(result_folder, 'variants.vcf')
            shutil.copy2(vcf_output, final_vcf)
            
            logger.info("VCF file processed successfully")
            update_progress(analysis_id, "variants", "complete", vcf_output, percent=100)
            
            # Step 11: Output VCF (when starting with VCF file)
            update_progress(analysis_id, "vcf", "running", percent=50,
                          note="Finalizing VCF file")
            
            # Any VCF post-processing could be added here
            
            logger.info("Pipeline complete. VCF file ready.")
            update_progress(analysis_id, "vcf", "complete", final_vcf, percent=100)
        
        # Pipeline completed
        logger.info(f"Pipeline completed for analysis {analysis_id}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        
        # Update current step as failed
        progress_file = os.path.join(result_folder, 'progress.json')
        if os.path.exists(progress_file):
            with open(progress_file, 'r') as f:
                progress = json.load(f)
            
            current_step_idx = progress.get("current_step", 0)
            if current_step_idx < len(progress["steps"]):
                current_step = progress["steps"][current_step_idx]
                update_progress(analysis_id, current_step["id"], "error")
        
        # Ensure the log reflects the error
        log_file = os.path.join(result_folder, 'pipeline.log')
        with open(log_file, 'a') as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Pipeline error: {str(e)}\n")

def continue_pipeline(analysis_id, start_step):
    """Continue an existing analysis from a specific step."""
    threading.Thread(target=_continue_pipeline_thread, args=(analysis_id, start_step)).start()

def _continue_pipeline_thread(analysis_id, start_step):
    """Continue the pipeline from a specific step in a separate thread."""
    try:
        logger.info(f"Continuing pipeline for analysis {analysis_id} from step {start_step}")
        
        # Get result folder
        result_folder = os.path.join(os.getcwd(), 'results', analysis_id)
        if not os.path.exists(result_folder):
            logger.error(f"Result folder for analysis {analysis_id} not found")
            return
        
        # Load existing progress
        progress_file = os.path.join(result_folder, 'progress.json')
        if os.path.exists(progress_file):
            with open(progress_file, 'r') as f:
                progress = json.load(f)
        else:
            logger.error(f"Progress file for analysis {analysis_id} not found")
            return
        
        # Reset steps after the start_step
        step_ids = [step["id"] for step in PIPELINE_STEPS]
        start_index = step_ids.index(start_step) if start_step in step_ids else -1
        
        if start_index == -1:
            logger.error(f"Invalid start step: {start_step}")
            return
        
        # Reset progress for steps after the start step
        for step in progress["steps"]:
            if step_ids.index(step["id"]) >= start_index:
                if step["id"] == start_step:
                    step["status"] = "running"
                    step["percent"] = 0
                    step["note"] = "Restarting step"
                else:
                    step["status"] = "pending"
                    step["percent"] = 0
                    if "note" in step:
                        del step["note"]
        
        # Save updated progress
        with open(progress_file, 'w') as f:
            json.dump(progress, f)
        
        # Log the continuation
        log_file = os.path.join(result_folder, 'pipeline.log')
        with open(log_file, 'a') as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Continuing pipeline from step {start_step}\n")
        
        # Get required paths
        align_output_dir = os.path.join(result_folder, 'aligned')
        sample_name = f"sample_{analysis_id}"
        
        # Continue from the specific step
        if start_step == "sam2bam":
            # Convert SAM to BAM
            sam_output = os.path.join(align_output_dir, f"{sample_name}.sam")
            bam_output = os.path.join(align_output_dir, f"{sample_name}.bam")
            sorted_bam = os.path.join(align_output_dir, f"{sample_name}.sorted.bam")
            
            update_progress(analysis_id, "sam2bam", "running", percent=0, 
                           note="Converting SAM to BAM format")
            
            cmd_samtools1 = f"samtools view -S -b {sam_output} > {bam_output}"
            logger.info(f"Running command: {cmd_samtools1}")
            update_progress(analysis_id, "sam2bam", "running", percent=50)
            
            process = subprocess.Popen(cmd_samtools1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                logger.error(f"Samtools view failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                update_progress(analysis_id, "sam2bam", "error")
                return
            
            logger.info("SAM to BAM conversion completed successfully")
            update_progress(analysis_id, "sam2bam", "complete", bam_output, percent=100)
            
            # Continue to the next steps
            _continue_from_sort(analysis_id, result_folder, align_output_dir, sample_name, bam_output, sorted_bam)
            
        elif start_step == "sort":
            # Sort BAM
            bam_output = os.path.join(align_output_dir, f"{sample_name}.bam")
            sorted_bam = os.path.join(align_output_dir, f"{sample_name}.sorted.bam")
            
            _continue_from_sort(analysis_id, result_folder, align_output_dir, sample_name, bam_output, sorted_bam)
            
        elif start_step == "dedup":
            # Mark Duplicates
            sorted_bam = os.path.join(align_output_dir, f"{sample_name}.sorted.bam")
            
            _continue_from_dedup(analysis_id, result_folder, sample_name, sorted_bam)
            
        elif start_step == "bqsr":
            # Base Quality Score Recalibration
            dedup_output_dir = os.path.join(result_folder, 'dedup')
            dedup_bam = os.path.join(dedup_output_dir, f"{sample_name}.dedup.bam")
            
            if not os.path.exists(dedup_bam):
                logger.error(f"Deduplication BAM file not found: {dedup_bam}")
                update_progress(analysis_id, "bqsr", "error", 
                               note="Required deduplication BAM file not found")
                return
            
            # Start from base recalibration
            update_progress(analysis_id, "bqsr", "running", percent=0,
                          note="Running Base Quality Score Recalibration")
            
            bqsr_output_dir = os.path.join(result_folder, 'bqsr')
            os.makedirs(bqsr_output_dir, exist_ok=True)
            
            # First, create recalibration table
            recal_table = os.path.join(bqsr_output_dir, f"{sample_name}.recal.table")
            bqsr_bam = os.path.join(bqsr_output_dir, f"{sample_name}.bqsr.bam")
            
            # Path to known variant sites (dbSNP)
            dbsnp_vcf = os.path.join(os.getcwd(), 'reference', 'dbsnp_138.hg19.vcf')
            
            # Check if dbSNP file exists
            if not os.path.exists(dbsnp_vcf):
                logger.warning(f"dbSNP file not found at {dbsnp_vcf}. Skipping BQSR.")
                update_progress(analysis_id, "bqsr", "skipped", percent=100,
                              note="Skipped - dbSNP file not available")
                
                # Continue with variant calling
                _continue_from_variants(analysis_id, result_folder, dedup_bam)
            else:
                # GATK BaseRecalibrator command to create recalibration table
                gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
                
                cmd_bqsr1 = f"{gatk_path} BaseRecalibrator \
                            -R {REFERENCE_GENOME} \
                            -I {dedup_bam} \
                            --known-sites {dbsnp_vcf} \
                            -O {recal_table}"
                
                logger.info(f"Running command: {cmd_bqsr1}")
                update_progress(analysis_id, "bqsr", "running", percent=25, 
                              note="Generating recalibration table")
                
                # Start the process
                process = subprocess.Popen(cmd_bqsr1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         text=True, bufsize=1, universal_newlines=True)
                
                # Monitor output in real-time
                start_time = time.time()
                last_update_time = start_time
                
                while process.poll() is None:
                    stderr_line = process.stderr.readline()
                    if stderr_line:
                        logger.info(f"GATK BaseRecalibrator: {stderr_line.strip()}")
                        
                        # Update progress periodically
                        current_time = time.time()
                        if current_time - last_update_time > 30:
                            elapsed_seconds = current_time - start_time
                            # Estimate progress based on elapsed time (assuming ~5 min for this step)
                            time_based_percent = min(45, 25 + int((elapsed_seconds / 300) * 20))
                            
                            update_progress(analysis_id, "bqsr", "running", percent=time_based_percent)
                            last_update_time = current_time
                    
                    stdout_line = process.stdout.readline()
                    if stdout_line:
                        logger.info(f"GATK BaseRecalibrator: {stdout_line.strip()}")
                    
                    time.sleep(0.1)
                
                # Get remaining output
                stdout, stderr = process.communicate()
                if stdout:
                    logger.info(f"GATK BaseRecalibrator stdout: {stdout}")
                if stderr:
                    logger.info(f"GATK BaseRecalibrator stderr: {stderr}")
                
                if process.returncode != 0:
                    logger.error(f"GATK BaseRecalibrator failed with return code {process.returncode}")
                    update_progress(analysis_id, "bqsr", "error")
                    return
                
                logger.info("Recalibration table generated successfully")
                update_progress(analysis_id, "bqsr", "running", percent=50, 
                              note="Applying base recalibration")
                
                # Apply recalibration with ApplyBQSR
                cmd_bqsr2 = f"{gatk_path} ApplyBQSR \
                            -R {REFERENCE_GENOME} \
                            -I {dedup_bam} \
                            --bqsr-recal-file {recal_table} \
                            -O {bqsr_bam}"
                
                logger.info(f"Running command: {cmd_bqsr2}")
                
                # Start the process
                process = subprocess.Popen(cmd_bqsr2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         text=True, bufsize=1, universal_newlines=True)
                
                # Monitor output in real-time
                start_time = time.time()
                last_update_time = start_time
                
                while process.poll() is None:
                    stderr_line = process.stderr.readline()
                    if stderr_line:
                        logger.info(f"GATK ApplyBQSR: {stderr_line.strip()}")
                        
                        # Update progress periodically
                        current_time = time.time()
                        if current_time - last_update_time > 30:
                            elapsed_seconds = current_time - start_time
                            # Estimate progress based on elapsed time (assuming ~5 min for this step)
                            time_based_percent = min(95, 50 + int((elapsed_seconds / 300) * 45))
                            
                            update_progress(analysis_id, "bqsr", "running", percent=time_based_percent)
                            last_update_time = current_time
                    
                    stdout_line = process.stdout.readline()
                    if stdout_line:
                        logger.info(f"GATK ApplyBQSR: {stdout_line.strip()}")
                    
                    time.sleep(0.1)
                
                # Get remaining output
                stdout, stderr = process.communicate()
                if stdout:
                    logger.info(f"GATK ApplyBQSR stdout: {stdout}")
                if stderr:
                    logger.info(f"GATK ApplyBQSR stderr: {stderr}")
                
                if process.returncode != 0:
                    logger.error(f"GATK ApplyBQSR failed with return code {process.returncode}")
                    update_progress(analysis_id, "bqsr", "error")
                    return
                
                # Create index for BQSR BAM
                cmd_index = f"samtools index {bqsr_bam}"
                logger.info(f"Running command: {cmd_index}")
                
                process = subprocess.Popen(cmd_index, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                
                if process.returncode != 0:
                    logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                    update_progress(analysis_id, "bqsr", "error", note="Failed to index BQSR BAM file")
                    return
                
                logger.info("Base quality score recalibration completed successfully")
                update_progress(analysis_id, "bqsr", "complete", bqsr_bam, percent=100)
                
                # Use BQSR BAM for variant calling
                _continue_from_variants(analysis_id, result_folder, bqsr_bam)
            
        elif start_step == "variants":
            # Variant Calling
            # First check for BQSR BAM
            bqsr_output_dir = os.path.join(result_folder, 'bqsr')
            bqsr_bam = os.path.join(bqsr_output_dir, f"{sample_name}.bqsr.bam")
            
            # If BQSR BAM exists, use it; otherwise use dedup BAM
            if os.path.exists(bqsr_bam):
                input_bam = bqsr_bam
            else:
                dedup_output_dir = os.path.join(result_folder, 'dedup')
                input_bam = os.path.join(dedup_output_dir, f"{sample_name}.dedup.bam")
            
            if not os.path.exists(input_bam):
                logger.error(f"Required BAM file not found: {input_bam}")
                update_progress(analysis_id, "variants", "error", 
                               note="Required BAM file not found")
                return
            
            _continue_from_variants(analysis_id, result_folder, input_bam)
        
        elif start_step == "filter":
            # Variant Filtering
            variants_output_dir = os.path.join(result_folder, 'variants')
            vcf_output = os.path.join(variants_output_dir, 'variants.vcf')
            
            if not os.path.exists(vcf_output):
                logger.error(f"Required VCF file not found: {vcf_output}")
                update_progress(analysis_id, "filter", "error", 
                               note="Required VCF file not found")
                return
            
            update_progress(analysis_id, "filter", "running", percent=0,
                          note="Filtering variants by quality and depth")
            
            filter_output_dir = os.path.join(result_folder, 'filtered')
            os.makedirs(filter_output_dir, exist_ok=True)
            
            filtered_vcf = os.path.join(filter_output_dir, 'filtered_variants.vcf')
            
            # GATK VariantFiltration command
            gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
            cmd_filter = f"{gatk_path} VariantFiltration \
                         -R {REFERENCE_GENOME} \
                         -V {vcf_output} \
                         -O {filtered_vcf} \
                         --filter-name 'QualFilter' --filter-expression 'QUAL < 30.0' \
                         --filter-name 'DepthFilter' --filter-expression 'DP < 10' \
                         --filter-name 'QDFilter' --filter-expression 'QD < 2.0'"
            
            logger.info(f"Running command: {cmd_filter}")
            update_progress(analysis_id, "filter", "running", percent=30)
            
            process = subprocess.Popen(cmd_filter, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                logger.error(f"GATK VariantFiltration failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
                update_progress(analysis_id, "filter", "error")
                return
            
            logger.info("Variant filtering completed successfully")
            update_progress(analysis_id, "filter", "complete", filtered_vcf, percent=100)
            
            # Continue with annotation
            # Call a continuation function that handles annotation, compression, and final output
            _continue_from_filter(analysis_id, result_folder, filtered_vcf)
            
        elif start_step == "annotate":
            # Variant Annotation
            # First check for filtered VCF
            filter_output_dir = os.path.join(result_folder, 'filtered')
            filtered_vcf = os.path.join(filter_output_dir, 'filtered_variants.vcf')
            
            # If filtered VCF exists, use it; otherwise use raw VCF
            if os.path.exists(filtered_vcf):
                input_vcf = filtered_vcf
            else:
                variants_output_dir = os.path.join(result_folder, 'variants')
                input_vcf = os.path.join(variants_output_dir, 'variants.vcf')
            
            if not os.path.exists(input_vcf):
                logger.error(f"Required VCF file not found: {input_vcf}")
                update_progress(analysis_id, "annotate", "error", 
                               note="Required VCF file not found")
                return
            
            # Continue with annotation
            _continue_from_annotation(analysis_id, result_folder, input_vcf)
            
        elif start_step == "compress":
            # VCF Compression
            # First check for annotated VCF
            annotate_output_dir = os.path.join(result_folder, 'annotated')
            annotated_vcf = os.path.join(annotate_output_dir, 'annotated_variants.vcf')
            
            # If annotated VCF exists, use it; otherwise check for filtered VCF
            if os.path.exists(annotated_vcf):
                input_vcf = annotated_vcf
            else:
                filter_output_dir = os.path.join(result_folder, 'filtered')
                filtered_vcf = os.path.join(filter_output_dir, 'filtered_variants.vcf')
                
                # If filtered VCF exists, use it; otherwise use raw VCF
                if os.path.exists(filtered_vcf):
                    input_vcf = filtered_vcf
                else:
                    variants_output_dir = os.path.join(result_folder, 'variants')
                    input_vcf = os.path.join(variants_output_dir, 'variants.vcf')
            
            if not os.path.exists(input_vcf):
                logger.error(f"Required VCF file not found: {input_vcf}")
                update_progress(analysis_id, "compress", "error", 
                               note="Required VCF file not found")
                return
            
            # Continue with compression
            _continue_from_compression(analysis_id, result_folder, input_vcf)
            
        else:
            logger.error(f"Unsupported start step: {start_step}")
            return
            
    except Exception as e:
        logger.error(f"Continuation failed: {str(e)}")
        
        # Update current step as failed
        progress_file = os.path.join(result_folder, 'progress.json')
        if os.path.exists(progress_file):
            with open(progress_file, 'r') as f:
                progress = json.load(f)
            
            # Find the current running step and mark it as error
            for step in progress["steps"]:
                if step["status"] == "running":
                    update_progress(analysis_id, step["id"], "error", note=f"Error: {str(e)}")
                    break
        
        # Ensure the log reflects the error
        log_file = os.path.join(result_folder, 'pipeline.log')
        with open(log_file, 'a') as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Pipeline continuation error: {str(e)}\n")

def _continue_from_sort(analysis_id, result_folder, align_output_dir, sample_name, bam_output, sorted_bam):
    """Continue pipeline from the sort step."""
    update_progress(analysis_id, "sort", "running", percent=0,
                   note="Sorting BAM file by coordinates")
    
    cmd_samtools2 = f"samtools sort {bam_output} -o {sorted_bam}"
    logger.info(f"Running command: {cmd_samtools2}")
    update_progress(analysis_id, "sort", "running", percent=50)
    
    process = subprocess.Popen(cmd_samtools2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        logger.error(f"Samtools sort failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
        update_progress(analysis_id, "sort", "error")
        return
    
    logger.info("BAM sorting completed successfully")
    update_progress(analysis_id, "sort", "complete", sorted_bam, percent=100)
    
    # Index BAM
    update_progress(analysis_id, "index", "running", percent=0,
                   note="Indexing sorted BAM file")
    
    cmd_samtools3 = f"samtools index {sorted_bam}"
    logger.info(f"Running command: {cmd_samtools3}")
    update_progress(analysis_id, "index", "running", percent=50)
    
    process = subprocess.Popen(cmd_samtools3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
        update_progress(analysis_id, "index", "error")
        return
    
    logger.info("BAM indexing completed successfully")
    update_progress(analysis_id, "index", "complete", f"{sorted_bam}.bai", percent=100)
    
    # Continue with deduplication
    _continue_from_dedup(analysis_id, result_folder, sample_name, sorted_bam)

def _continue_from_dedup(analysis_id, result_folder, sample_name, sorted_bam):
    """Continue pipeline from the duplicate marking step."""
    # Mark Duplicates
    update_progress(analysis_id, "dedup", "running", percent=0,
                  note="Marking duplicate reads using Picard")
    
    dedup_output_dir = os.path.join(result_folder, 'dedup')
    os.makedirs(dedup_output_dir, exist_ok=True)
    
    # Output files
    dedup_bam = os.path.join(dedup_output_dir, f"{sample_name}.dedup.bam")
    metrics_file = os.path.join(dedup_output_dir, f"{sample_name}.dedup.metrics.txt")
    
    # Use standalone Picard tool instead of GATK
    picard_jar = os.path.join(os.getcwd(), 'tools', 'picard.jar')
    
    # Picard MarkDuplicates command
    cmd_dedup = f"java -jar {picard_jar} MarkDuplicates \
                I={sorted_bam} \
                O={dedup_bam} \
                M={metrics_file} \
                VALIDATION_STRINGENCY=LENIENT \
                CREATE_INDEX=true"
    
    logger.info(f"Running command: {cmd_dedup}")
    update_progress(analysis_id, "dedup", "running", percent=25, 
                  note="This process may take 5-15 minutes depending on file size")
    
    # Start the process
    process = subprocess.Popen(cmd_dedup, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True, bufsize=1, universal_newlines=True)
    
    # Monitor output in real-time and estimate progress
    start_time = time.time()
    last_update_time = start_time
    
    while process.poll() is None:
        stderr_line = process.stderr.readline()
        if stderr_line:
            logger.info(f"Picard MarkDuplicates: {stderr_line.strip()}")
            
            # Update progress periodically based on log messages
            current_time = time.time()
            elapsed_seconds = current_time - start_time
            
            # Update progress every 30 seconds to avoid too many updates
            if current_time - last_update_time > 30:
                # Estimate progress based on elapsed time (assuming ~10 min total)
                time_based_percent = min(95, int((elapsed_seconds / 600) * 100))
                current_percent = max(25, time_based_percent)
                
                update_progress(analysis_id, "dedup", "running", percent=current_percent,
                              note=f"Marking duplicates (running for {int(elapsed_seconds)} seconds)")
                
                last_update_time = current_time
        
        stdout_line = process.stdout.readline()
        if stdout_line:
            logger.info(f"Picard MarkDuplicates: {stdout_line.strip()}")
        
        time.sleep(0.1)
    
    # Get remaining output
    stdout, stderr = process.communicate()
    if stdout:
        logger.info(f"Picard MarkDuplicates stdout: {stdout}")
    if stderr:
        logger.info(f"Picard MarkDuplicates stderr: {stderr}")
    
    if process.returncode != 0:
        logger.error(f"Picard MarkDuplicates failed with return code {process.returncode}")
        update_progress(analysis_id, "dedup", "error")
        return
    
    logger.info("Duplicate marking completed successfully")
    update_progress(analysis_id, "dedup", "complete", dedup_bam, percent=100)
    
    # Step 8: Index the BAM file
    update_progress(analysis_id, "index", "running", percent=0,
                  note="Indexing BAM file")
    
    # Check if index was already created by Picard
    bai_file = f"{dedup_bam}.bai"
    if not os.path.exists(bai_file):
        cmd_index = f"samtools index {dedup_bam}"
        logger.info(f"Running command: {cmd_index}")
        
        process = subprocess.Popen(cmd_index, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
            update_progress(analysis_id, "index", "error")
            return
    else:
        logger.info("BAM index already exists from Picard")
    
    logger.info("BAM indexing completed successfully")
    update_progress(analysis_id, "index", "complete", f"{dedup_bam}.bai", percent=100)
    
    # Step 9: Base Quality Score Recalibration (BQSR)
    update_progress(analysis_id, "bqsr", "running", percent=0,
                  note="Running Base Quality Score Recalibration")
    
    bqsr_output_dir = os.path.join(result_folder, 'bqsr')
    os.makedirs(bqsr_output_dir, exist_ok=True)
    
    # First, create recalibration table
    recal_table = os.path.join(bqsr_output_dir, f"{sample_name}.recal.table")
    bqsr_bam = os.path.join(bqsr_output_dir, f"{sample_name}.bqsr.bam")
    
    # Path to known variant sites (dbSNP)
    dbsnp_vcf = os.path.join(os.getcwd(), 'reference', 'dbsnp_138.hg19.vcf')
    
    # Check if dbSNP file exists
    if not os.path.exists(dbsnp_vcf):
        logger.warning(f"dbSNP file not found at {dbsnp_vcf}. Skipping BQSR.")
        update_progress(analysis_id, "bqsr", "skipped", percent=100,
                      note="Skipped - dbSNP file not available")
    else:
        # GATK BaseRecalibrator command to create recalibration table
        gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
        
        cmd_bqsr1 = f"{gatk_path} BaseRecalibrator \
                    -R {REFERENCE_GENOME} \
                    -I {dedup_bam} \
                    --known-sites {dbsnp_vcf} \
                    -O {recal_table}"
        
        logger.info(f"Running command: {cmd_bqsr1}")
        update_progress(analysis_id, "bqsr", "running", percent=25, 
                      note="Generating recalibration table")
        
        # Start the process
        process = subprocess.Popen(cmd_bqsr1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, universal_newlines=True)
        
        # Monitor output in real-time
        start_time = time.time()
        last_update_time = start_time
        
        while process.poll() is None:
            stderr_line = process.stderr.readline()
            if stderr_line:
                logger.info(f"GATK BaseRecalibrator: {stderr_line.strip()}")
                
                # Update progress periodically
                current_time = time.time()
                if current_time - last_update_time > 30:
                    elapsed_seconds = current_time - start_time
                    # Estimate progress based on elapsed time (assuming ~5 min for this step)
                    time_based_percent = min(45, 25 + int((elapsed_seconds / 300) * 20))
                    
                    update_progress(analysis_id, "bqsr", "running", percent=time_based_percent)
                    last_update_time = current_time
            
            stdout_line = process.stdout.readline()
            if stdout_line:
                logger.info(f"GATK BaseRecalibrator: {stdout_line.strip()}")
            
            time.sleep(0.1)
        
        # Get remaining output
        stdout, stderr = process.communicate()
        if stdout:
            logger.info(f"GATK BaseRecalibrator stdout: {stdout}")
        if stderr:
            logger.info(f"GATK BaseRecalibrator stderr: {stderr}")
        
        if process.returncode != 0:
            logger.error(f"GATK BaseRecalibrator failed with return code {process.returncode}")
            update_progress(analysis_id, "bqsr", "error")
            return
        
        logger.info("Recalibration table generated successfully")
        update_progress(analysis_id, "bqsr", "running", percent=50, 
                      note="Applying base recalibration")
        
        # Apply recalibration with ApplyBQSR
        cmd_bqsr2 = f"{gatk_path} ApplyBQSR \
                    -R {REFERENCE_GENOME} \
                    -I {dedup_bam} \
                    --bqsr-recal-file {recal_table} \
                    -O {bqsr_bam}"
        
        logger.info(f"Running command: {cmd_bqsr2}")
        
        # Start the process
        process = subprocess.Popen(cmd_bqsr2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, universal_newlines=True)
        
        # Monitor output in real-time
        start_time = time.time()
        last_update_time = start_time
        
        while process.poll() is None:
            stderr_line = process.stderr.readline()
            if stderr_line:
                logger.info(f"GATK ApplyBQSR: {stderr_line.strip()}")
                
                # Update progress periodically
                current_time = time.time()
                if current_time - last_update_time > 30:
                    elapsed_seconds = current_time - start_time
                    # Estimate progress based on elapsed time (assuming ~5 min for this step)
                    time_based_percent = min(95, 50 + int((elapsed_seconds / 300) * 45))
                    
                    update_progress(analysis_id, "bqsr", "running", percent=time_based_percent)
                    last_update_time = current_time
            
            stdout_line = process.stdout.readline()
            if stdout_line:
                logger.info(f"GATK ApplyBQSR: {stdout_line.strip()}")
            
            time.sleep(0.1)
        
        # Get remaining output
        stdout, stderr = process.communicate()
        if stdout:
            logger.info(f"GATK ApplyBQSR stdout: {stdout}")
        if stderr:
            logger.info(f"GATK ApplyBQSR stderr: {stderr}")
        
        if process.returncode != 0:
            logger.error(f"GATK ApplyBQSR failed with return code {process.returncode}")
            update_progress(analysis_id, "bqsr", "error")
            return
        
        # Create index for BQSR BAM
        cmd_index = f"samtools index {bqsr_bam}"
        logger.info(f"Running command: {cmd_index}")
        
        process = subprocess.Popen(cmd_index, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Samtools index failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
            update_progress(analysis_id, "bqsr", "error", note="Failed to index BQSR BAM file")
            return
        
        logger.info("Base quality score recalibration completed successfully")
        update_progress(analysis_id, "bqsr", "complete", bqsr_bam, percent=100)
        
        # Use BQSR BAM for variant calling
        dedup_bam = bqsr_bam
    
    # Continue with variant calling
    _continue_from_variants(analysis_id, result_folder, dedup_bam)

def _continue_from_variants(analysis_id, result_folder, dedup_bam):
    """Continue pipeline from the variant calling step."""
    # Variant Calling with GATK
    update_progress(analysis_id, "variants", "running", percent=0,
                   note="Calling variants with GATK HaplotypeCaller")
    logger.info("Calling variants with GATK")
    
    variants_output_dir = os.path.join(result_folder, 'variants')
    os.makedirs(variants_output_dir, exist_ok=True)
    
    vcf_output = os.path.join(variants_output_dir, 'variants.vcf')
    
    # Ensure reference genome has all required indices
    if not ensure_reference_indices(result_folder, analysis_id):
        return
    
    # Verify BAM file compatibility with GATK
    if not verify_bam_for_gatk(dedup_bam, analysis_id, result_folder):
        logger.error(f"BAM file {dedup_bam} is not compatible with GATK and could not be fixed")
        update_progress(analysis_id, "variants", "error", 
                       note="BAM file is not compatible with GATK")
        return
    
    # Using HaplotypeCaller with deduplicated BAM
    gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
    cmd_gatk = f"{gatk_path} HaplotypeCaller \
                  -R {REFERENCE_GENOME} \
                  -I {dedup_bam} \
                  -O {vcf_output} \
                  --standard-min-confidence-threshold-for-calling 30 \
                  --minimum-mapping-quality 20"
    
    logger.info(f"Running command: {cmd_gatk}")
    
    # Start the process
    process = subprocess.Popen(cmd_gatk, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True, bufsize=1, universal_newlines=True)
    
    # Monitor output in real-time and estimate progress
    start_time = time.time()
    last_update_time = start_time
    
    while process.poll() is None:
        stderr_line = process.stderr.readline()
        if stderr_line:
            logger.info(f"GATK HaplotypeCaller: {stderr_line.strip()}")
            
            # Update progress periodically
            current_time = time.time()
            if current_time - last_update_time > 30:
                elapsed_seconds = current_time - start_time
                
                # For HaplotypeCaller, try to detect progress from the output, but fall back to time-based
                progress_match = None
                
                if progress_match:
                    current_percent = int(float(progress_match.group(1)))
                else:
                    # Estimate progress based on elapsed time (assume ~5 min for variant calling)
                    time_based_percent = min(99, int((elapsed_seconds / 300) * 100))  # Assume ~5 min for variant calling
                    current_percent = time_based_percent
                
                update_progress(analysis_id, "variants", "running", percent=current_percent)
                last_update_time = current_time
                
                logger.info(f"GATK variant calling is still running... ~{current_percent}% complete (elapsed: {int(elapsed_seconds)} seconds)")
        
        stdout_line = process.stdout.readline()
        if stdout_line:
            logger.info(f"GATK HaplotypeCaller: {stdout_line.strip()}")
        
        time.sleep(0.1)
    
    # Get remaining output
    stdout, stderr = process.communicate()
    if stdout:
        logger.info(f"GATK HaplotypeCaller stdout: {stdout}")
    if stderr:
        logger.info(f"GATK HaplotypeCaller stderr: {stderr}")
    
    if process.returncode != 0:
        logger.error(f"GATK variant calling failed with return code {process.returncode}")
        update_progress(analysis_id, "variants", "error")
        return
    
    # Copy output to root results directory for easy access
    final_vcf = os.path.join(result_folder, 'variants.vcf')
    shutil.copy2(vcf_output, final_vcf)
    
    logger.info("Variant calling completed successfully")
    update_progress(analysis_id, "variants", "complete", vcf_output, percent=100)
    
    # Step 11: Variant Filtering
    update_progress(analysis_id, "filter", "running", percent=0,
                  note="Filtering variants by quality and depth")
    
    filter_output_dir = os.path.join(result_folder, 'filtered')
    os.makedirs(filter_output_dir, exist_ok=True)
    
    filtered_vcf = os.path.join(filter_output_dir, 'filtered_variants.vcf')
    
    # GATK VariantFiltration command
    cmd_filter = f"{gatk_path} VariantFiltration \
                 -R {REFERENCE_GENOME} \
                 -V {vcf_output} \
                 -O {filtered_vcf} \
                 --filter-name 'QualFilter' --filter-expression 'QUAL < 30.0' \
                 --filter-name 'DepthFilter' --filter-expression 'DP < 10' \
                 --filter-name 'QDFilter' --filter-expression 'QD < 2.0'"
    
    logger.info(f"Running command: {cmd_filter}")
    update_progress(analysis_id, "filter", "running", percent=30)
    
    process = subprocess.Popen(cmd_filter, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True, bufsize=1, universal_newlines=True)
    
    # Monitor output
    start_time = time.time()
    last_update_time = start_time
    
    while process.poll() is None:
        stderr_line = process.stderr.readline()
        if stderr_line:
            logger.info(f"GATK VariantFiltration: {stderr_line.strip()}")
            
            # Update progress periodically
            current_time = time.time()
            if current_time - last_update_time > 15:
                elapsed_seconds = current_time - start_time
                # Assume filtering takes about 1 minute total
                time_based_percent = min(95, 30 + int((elapsed_seconds / 60) * 65))
                
                update_progress(analysis_id, "filter", "running", percent=time_based_percent)
                last_update_time = current_time
        
        stdout_line = process.stdout.readline()
        if stdout_line:
            logger.info(f"GATK VariantFiltration: {stdout_line.strip()}")
        
        time.sleep(0.1)
    
    # Get remaining output
    stdout, stderr = process.communicate()
    if stdout:
        logger.info(f"GATK VariantFiltration stdout: {stdout}")
    if stderr:
        logger.info(f"GATK VariantFiltration stderr: {stderr}")
    
    if process.returncode != 0:
        logger.error(f"GATK VariantFiltration failed with return code {process.returncode}")
        update_progress(analysis_id, "filter", "error")
        # Continue to next step even if this fails
    else:
        logger.info("Variant filtering completed successfully")
        update_progress(analysis_id, "filter", "complete", filtered_vcf, percent=100)
        
        # Use filtered VCF for subsequent steps
        vcf_output = filtered_vcf
    
    # Step 12: Variant Annotation
    update_progress(analysis_id, "annotate", "running", percent=0,
                  note="Annotating variants with functional information")
    
    # Check if a suitable annotation tool is available
    # For this example, we're using GATK Funcotator, but in practice might use tools like SnpEff, VEP, etc.
    annotate_output_dir = os.path.join(result_folder, 'annotated')
    os.makedirs(annotate_output_dir, exist_ok=True)
    
    annotated_vcf = os.path.join(annotate_output_dir, 'annotated_variants.vcf')
    
    # Check for data sources directory for annotation
    data_sources_dir = os.path.join(os.getcwd(), 'reference', 'funcotator_dataSources')
    
    if not os.path.exists(data_sources_dir):
        logger.warning(f"Annotation data sources not found at {data_sources_dir}. Skipping annotation.")
        update_progress(analysis_id, "annotate", "skipped", percent=100,
                      note="Skipped - annotation data sources not available")
    else:
        # GATK Funcotator command
        cmd_annotate = f"{gatk_path} Funcotator \
                       --variant {vcf_output} \
                       --reference {REFERENCE_GENOME} \
                       --data-sources-path {data_sources_dir} \
                       --output {annotated_vcf} \
                       --output-file-format VCF"
        
        logger.info(f"Running command: {cmd_annotate}")
        update_progress(analysis_id, "annotate", "running", percent=30,
                      note="Annotating variants (this may take 5-10 minutes)")
        
        process = subprocess.Popen(cmd_annotate, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, universal_newlines=True)
        
        # Monitor output
        start_time = time.time()
        last_update_time = start_time
        
        while process.poll() is None:
            stderr_line = process.stderr.readline()
            if stderr_line:
                logger.info(f"GATK Funcotator: {stderr_line.strip()}")
                
                # Update progress periodically
                current_time = time.time()
                if current_time - last_update_time > 30:
                    elapsed_seconds = current_time - start_time
                    # Assume annotation takes about 5-10 minutes
                    time_based_percent = min(95, 30 + int((elapsed_seconds / 600) * 65))
                    
                    update_progress(analysis_id, "annotate", "running", percent=time_based_percent)
                    last_update_time = current_time
            
            stdout_line = process.stdout.readline()
            if stdout_line:
                logger.info(f"GATK Funcotator: {stdout_line.strip()}")
            
            time.sleep(0.1)
        
        # Get remaining output
        stdout, stderr = process.communicate()
        if stdout:
            logger.info(f"GATK Funcotator stdout: {stdout}")
        if stderr:
            logger.info(f"GATK Funcotator stderr: {stderr}")
        
        if process.returncode != 0:
            logger.error(f"GATK Funcotator failed with return code {process.returncode}")
            update_progress(analysis_id, "annotate", "error")
            # Continue to next step even if this fails
        else:
            logger.info("Variant annotation completed successfully")
            update_progress(analysis_id, "annotate", "complete", annotated_vcf, percent=100)
            
            # Use annotated VCF for subsequent steps
            vcf_output = annotated_vcf
    
    # Step 13: VCF Compression and Indexing
    update_progress(analysis_id, "compress", "running", percent=0,
                  note="Compressing and indexing VCF file")
    
    compress_output_dir = os.path.join(result_folder, 'compressed')
    os.makedirs(compress_output_dir, exist_ok=True)
    
    compressed_vcf = os.path.join(compress_output_dir, 'variants.vcf.gz')
    
    # Compress VCF with bgzip
    cmd_bgzip = f"bgzip -c {vcf_output} > {compressed_vcf}"
    logger.info(f"Running command: {cmd_bgzip}")
    update_progress(analysis_id, "compress", "running", percent=30)
    
    process = subprocess.Popen(cmd_bgzip, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        logger.error(f"bgzip compression failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
        update_progress(analysis_id, "compress", "error")
        # Continue to final step even if this fails
    else:
        # Index the compressed VCF with tabix
        cmd_tabix = f"tabix -p vcf {compressed_vcf}"
        logger.info(f"Running command: {cmd_tabix}")
        update_progress(analysis_id, "compress", "running", percent=70)
        
        process = subprocess.Popen(cmd_tabix, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"tabix indexing failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
            update_progress(analysis_id, "compress", "error")
        else:
            logger.info("VCF compression and indexing completed successfully")
            update_progress(analysis_id, "compress", "complete", compressed_vcf, percent=100)
            
            # Use compressed VCF for final step
            vcf_output = compressed_vcf
    
    # Step 14: Output Final VCF
    update_progress(analysis_id, "vcf", "running", percent=50,
                  note="Finalizing VCF file")
    
    # Copy final VCF to the root results directory for easy access
    final_vcf = os.path.join(result_folder, 'variants.vcf')
    
    # Check if the final output is compressed
    if vcf_output.endswith('.gz'):
        # For compressed VCF, we'll create both compressed and uncompressed copies
        try:
            # Create uncompressed copy for easier viewing
            with gzip.open(vcf_output, 'rb') as f_in:
                with open(final_vcf, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Also copy the compressed version and its index
            shutil.copy2(vcf_output, os.path.join(result_folder, 'variants.vcf.gz'))
            if os.path.exists(f"{vcf_output}.tbi"):
                shutil.copy2(f"{vcf_output}.tbi", os.path.join(result_folder, 'variants.vcf.gz.tbi'))
        except Exception as e:
            logger.error(f"Error creating uncompressed VCF copy: {str(e)}")
            # Just copy the compressed file as-is
            shutil.copy2(vcf_output, os.path.join(result_folder, 'variants.vcf.gz'))
    else:
        # For uncompressed VCF, just copy it
        shutil.copy2(vcf_output, final_vcf)
    
    # Create a summary file with variant counts
    try:
        if os.path.exists(final_vcf):
            # Count variants using grep
            cmd_count = f"grep -v '^#' {final_vcf} | wc -l"
            process = subprocess.Popen(cmd_count, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode == 0:
                variant_count = stdout.decode().strip() if isinstance(stdout, bytes) else stdout.strip()
                
                summary_file = os.path.join(result_folder, 'variant_summary.txt')
                with open(summary_file, 'w') as f:
                    f.write(f"Analysis ID: {analysis_id}\n")
                    f.write(f"Total variants: {variant_count}\n")
                    f.write(f"Date completed: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    except Exception as e:
        logger.error(f"Error creating variant summary: {str(e)}")
    
    logger.info("Pipeline complete. VCF file ready.")
    update_progress(analysis_id, "vcf", "complete", final_vcf, percent=100)

def _continue_from_filter(analysis_id, result_folder, filtered_vcf):
    """Continue pipeline from the variant filtering step."""
    sample_name = f"sample_{analysis_id}"
    
    # Step 12: Variant Annotation
    update_progress(analysis_id, "annotate", "running", percent=0,
                  note="Annotating variants with functional information")
    
    # Check if a suitable annotation tool is available
    annotate_output_dir = os.path.join(result_folder, 'annotated')
    os.makedirs(annotate_output_dir, exist_ok=True)
    
    annotated_vcf = os.path.join(annotate_output_dir, 'annotated_variants.vcf')
    
    # Check for data sources directory for annotation
    data_sources_dir = os.path.join(os.getcwd(), 'reference', 'funcotator_dataSources')
    
    if not os.path.exists(data_sources_dir):
        logger.warning(f"Annotation data sources not found at {data_sources_dir}. Skipping annotation.")
        update_progress(analysis_id, "annotate", "skipped", percent=100,
                      note="Skipped - annotation data sources not available")
        
        # Continue with compression
        _continue_from_compression(analysis_id, result_folder, filtered_vcf)
    else:
        # GATK Funcotator command
        gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
        cmd_annotate = f"{gatk_path} Funcotator \
                       --variant {filtered_vcf} \
                       --reference {REFERENCE_GENOME} \
                       --data-sources-path {data_sources_dir} \
                       --output {annotated_vcf} \
                       --output-file-format VCF"
        
        logger.info(f"Running command: {cmd_annotate}")
        update_progress(analysis_id, "annotate", "running", percent=30,
                      note="Annotating variants (this may take 5-10 minutes)")
        
        process = subprocess.Popen(cmd_annotate, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, universal_newlines=True)
        
        # Monitor output
        start_time = time.time()
        last_update_time = start_time
        
        while process.poll() is None:
            stderr_line = process.stderr.readline()
            if stderr_line:
                logger.info(f"GATK Funcotator: {stderr_line.strip()}")
                
                # Update progress periodically
                current_time = time.time()
                if current_time - last_update_time > 30:
                    elapsed_seconds = current_time - start_time
                    # Assume annotation takes about 5-10 minutes
                    time_based_percent = min(95, 30 + int((elapsed_seconds / 600) * 65))
                    
                    update_progress(analysis_id, "annotate", "running", percent=time_based_percent)
                    last_update_time = current_time
            
            stdout_line = process.stdout.readline()
            if stdout_line:
                logger.info(f"GATK Funcotator: {stdout_line.strip()}")
            
            time.sleep(0.1)
        
        # Get remaining output
        stdout, stderr = process.communicate()
        if stdout:
            logger.info(f"GATK Funcotator stdout: {stdout}")
        if stderr:
            logger.info(f"GATK Funcotator stderr: {stderr}")
        
        if process.returncode != 0:
            logger.error(f"GATK Funcotator failed with return code {process.returncode}")
            update_progress(analysis_id, "annotate", "error")
            
            # Continue with compression using the filtered VCF
            _continue_from_compression(analysis_id, result_folder, filtered_vcf)
        else:
            logger.info("Variant annotation completed successfully")
            update_progress(analysis_id, "annotate", "complete", annotated_vcf, percent=100)
            
            # Continue with compression using the annotated VCF
            _continue_from_compression(analysis_id, result_folder, annotated_vcf)

def _continue_from_annotation(analysis_id, result_folder, input_vcf):
    """Continue pipeline from the annotation step."""
    sample_name = f"sample_{analysis_id}"
    
    # Check for data sources directory for annotation
    data_sources_dir = os.path.join(os.getcwd(), 'reference', 'funcotator_dataSources')
    
    if not os.path.exists(data_sources_dir):
        logger.warning(f"Annotation data sources not found at {data_sources_dir}. Skipping annotation.")
        update_progress(analysis_id, "annotate", "skipped", percent=100,
                      note="Skipped - annotation data sources not available")
        
        # Continue with compression
        _continue_from_compression(analysis_id, result_folder, input_vcf)
    else:
        # Step 12: Variant Annotation
        update_progress(analysis_id, "annotate", "running", percent=0,
                      note="Annotating variants with functional information")
        
        annotate_output_dir = os.path.join(result_folder, 'annotated')
        os.makedirs(annotate_output_dir, exist_ok=True)
        
        annotated_vcf = os.path.join(annotate_output_dir, 'annotated_variants.vcf')
        
        # GATK Funcotator command
        gatk_path = os.path.join(os.getcwd(), 'gatk-4.6.2.0', 'gatk')
        cmd_annotate = f"{gatk_path} Funcotator \
                       --variant {input_vcf} \
                       --reference {REFERENCE_GENOME} \
                       --data-sources-path {data_sources_dir} \
                       --output {annotated_vcf} \
                       --output-file-format VCF"
        
        logger.info(f"Running command: {cmd_annotate}")
        update_progress(analysis_id, "annotate", "running", percent=30,
                      note="Annotating variants (this may take 5-10 minutes)")
        
        process = subprocess.Popen(cmd_annotate, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, universal_newlines=True)
        
        # Monitor output
        start_time = time.time()
        last_update_time = start_time
        
        while process.poll() is None:
            stderr_line = process.stderr.readline()
            if stderr_line:
                logger.info(f"GATK Funcotator: {stderr_line.strip()}")
                
                # Update progress periodically
                current_time = time.time()
                if current_time - last_update_time > 30:
                    elapsed_seconds = current_time - start_time
                    # Assume annotation takes about 5-10 minutes
                    time_based_percent = min(95, 30 + int((elapsed_seconds / 600) * 65))
                    
                    update_progress(analysis_id, "annotate", "running", percent=time_based_percent)
                    last_update_time = current_time
            
            stdout_line = process.stdout.readline()
            if stdout_line:
                logger.info(f"GATK Funcotator: {stdout_line.strip()}")
            
            time.sleep(0.1)
        
        # Get remaining output
        stdout, stderr = process.communicate()
        if stdout:
            logger.info(f"GATK Funcotator stdout: {stdout}")
        if stderr:
            logger.info(f"GATK Funcotator stderr: {stderr}")
        
        if process.returncode != 0:
            logger.error(f"GATK Funcotator failed with return code {process.returncode}")
            update_progress(analysis_id, "annotate", "error")
            
            # Continue with compression using the input VCF
            _continue_from_compression(analysis_id, result_folder, input_vcf)
        else:
            logger.info("Variant annotation completed successfully")
            update_progress(analysis_id, "annotate", "complete", annotated_vcf, percent=100)
            
            # Continue with compression using the annotated VCF
            _continue_from_compression(analysis_id, result_folder, annotated_vcf)

def _continue_from_compression(analysis_id, result_folder, input_vcf):
    """Continue pipeline from the VCF compression step."""
    sample_name = f"sample_{analysis_id}"
    
    # Step 13: VCF Compression and Indexing
    update_progress(analysis_id, "compress", "running", percent=0,
                  note="Compressing and indexing VCF file")
    
    compress_output_dir = os.path.join(result_folder, 'compressed')
    os.makedirs(compress_output_dir, exist_ok=True)
    
    compressed_vcf = os.path.join(compress_output_dir, 'variants.vcf.gz')
    
    # Compress VCF with bgzip
    cmd_bgzip = f"bgzip -c {input_vcf} > {compressed_vcf}"
    logger.info(f"Running command: {cmd_bgzip}")
    update_progress(analysis_id, "compress", "running", percent=30)
    
    process = subprocess.Popen(cmd_bgzip, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        logger.error(f"bgzip compression failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
        update_progress(analysis_id, "compress", "error")
        
        # Continue to final step with the uncompressed VCF
        _finalize_vcf(analysis_id, result_folder, input_vcf)
    else:
        # Index the compressed VCF with tabix
        cmd_tabix = f"tabix -p vcf {compressed_vcf}"
        logger.info(f"Running command: {cmd_tabix}")
        update_progress(analysis_id, "compress", "running", percent=70)
        
        process = subprocess.Popen(cmd_tabix, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"tabix indexing failed: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
            update_progress(analysis_id, "compress", "error")
            
            # Continue to final step with the compressed but unindexed VCF
            _finalize_vcf(analysis_id, result_folder, compressed_vcf)
        else:
            logger.info("VCF compression and indexing completed successfully")
            update_progress(analysis_id, "compress", "complete", compressed_vcf, percent=100)
            
            # Continue to final step with the compressed and indexed VCF
            _finalize_vcf(analysis_id, result_folder, compressed_vcf)

def _finalize_vcf(analysis_id, result_folder, vcf_output):
    """Final step to prepare the VCF file for download."""
    # Step 14: Output Final VCF
    update_progress(analysis_id, "vcf", "running", percent=50,
                  note="Finalizing VCF file")
    
    # Copy final VCF to the root results directory for easy access
    final_vcf = os.path.join(result_folder, 'variants.vcf')
    
    # Check if the final output is compressed
    if vcf_output.endswith('.gz'):
        # For compressed VCF, we'll create both compressed and uncompressed copies
        try:
            import gzip
            # Create uncompressed copy for easier viewing
            with gzip.open(vcf_output, 'rb') as f_in:
                with open(final_vcf, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Also copy the compressed version and its index
            shutil.copy2(vcf_output, os.path.join(result_folder, 'variants.vcf.gz'))
            if os.path.exists(f"{vcf_output}.tbi"):
                shutil.copy2(f"{vcf_output}.tbi", os.path.join(result_folder, 'variants.vcf.gz.tbi'))
        except Exception as e:
            logger.error(f"Error creating uncompressed VCF copy: {str(e)}")
            # Just copy the compressed file as-is
            shutil.copy2(vcf_output, os.path.join(result_folder, 'variants.vcf.gz'))
    else:
        # For uncompressed VCF, just copy it
        shutil.copy2(vcf_output, final_vcf)
    
    # Create a summary file with variant counts
    try:
        if os.path.exists(final_vcf):
            # Count variants using grep
            cmd_count = f"grep -v '^#' {final_vcf} | wc -l"
            process = subprocess.Popen(cmd_count, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode == 0:
                variant_count = stdout.decode().strip() if isinstance(stdout, bytes) else stdout.strip()
                
                summary_file = os.path.join(result_folder, 'variant_summary.txt')
                with open(summary_file, 'w') as f:
                    f.write(f"Analysis ID: {analysis_id}\n")
                    f.write(f"Total variants: {variant_count}\n")
                    f.write(f"Date completed: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    except Exception as e:
        logger.error(f"Error creating variant summary: {str(e)}")
    
    logger.info("Pipeline complete. VCF file ready.")
    update_progress(analysis_id, "vcf", "complete", final_vcf, percent=100)