#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import logging
import argparse
import glob

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def fix_read_groups(analysis_id):
    """Fix a BAM file by adding read groups to make it compatible with GATK."""
    logger.info(f"Processing analysis {analysis_id}")
    
    # Set up paths
    results_folder = os.path.join(os.getcwd(), 'results', analysis_id)
    dedup_dir = os.path.join(results_folder, 'dedup')
    sample_name = f"sample_{analysis_id}"
    dedup_bam = os.path.join(dedup_dir, f"{sample_name}.dedup.bam")
    fixed_bam = os.path.join(dedup_dir, f"{sample_name}.fixed.bam")
    
    # Check if the BAM file exists
    if not os.path.exists(dedup_bam):
        logger.error(f"BAM file not found: {dedup_bam}")
        return False
    
    # Create read group information
    platform = "ILLUMINA"
    library = f"lib_{analysis_id}"
    read_group = f"ID:{analysis_id}\tSM:{analysis_id}\tPL:{platform}\tLB:{library}\tPU:unit1"
    
    # Use samtools to add read group
    cmd = f"samtools addreplacerg -r '{read_group}' -o {fixed_bam} {dedup_bam}"
    logger.info(f"Running command: {cmd}")
    
    # Execute the command
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        logger.error(f"Failed to add read groups: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
        return False
    
    # Replace the original BAM file with the fixed one
    cmd = f"mv {fixed_bam} {dedup_bam}"
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    
    # Index the fixed BAM
    cmd = f"samtools index {dedup_bam}"
    logger.info(f"Running command: {cmd}")
    
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        logger.error(f"Failed to index the fixed BAM: {stderr.decode() if isinstance(stderr, bytes) else stderr}")
        return False
    
    logger.info(f"Successfully added read groups to {dedup_bam}")
    return True

def process_all_analyses():
    """Process all analyses in the results directory."""
    results_dir = os.path.join(os.getcwd(), 'results')
    analyses = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
    
    fixed_count = 0
    failed_count = 0
    
    for analysis_id in analyses:
        if fix_read_groups(analysis_id):
            fixed_count += 1
        else:
            failed_count += 1
    
    logger.info(f"Fixed {fixed_count} BAM files, {failed_count} failures")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fix BAM files by adding read groups for GATK compatibility')
    parser.add_argument('-a', '--analysis', help='Specific analysis ID to fix')
    parser.add_argument('-all', '--all', action='store_true', help='Process all analyses')
    
    args = parser.parse_args()
    
    if args.all:
        process_all_analyses()
    elif args.analysis:
        fix_read_groups(args.analysis)
    else:
        parser.print_help() 