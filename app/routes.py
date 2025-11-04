#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import uuid
import time
import shutil
import json
import subprocess
from flask import Blueprint, render_template, request, redirect, url_for, flash, jsonify, send_from_directory, current_app

from app.pipeline.pipeline import run_pipeline

main_bp = Blueprint('main', __name__)

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in current_app.config['ALLOWED_EXTENSIONS'] or \
           filename.endswith('.fastq.gz') or filename.endswith('.fq.gz')

@main_bp.route('/')
def index():
    return render_template('index.html')

@main_bp.route('/upload', methods=['POST'])
def upload_file():
    # Check if files are in the request
    if 'r1_file' not in request.files or 'r2_file' not in request.files:
        flash('No file part', 'error')
        return redirect(request.url)
    
    r1_file = request.files['r1_file']
    r2_file = request.files['r2_file']
    
    # If user does not select file, browser also submits an empty part without filename
    if r1_file.filename == '' or r2_file.filename == '':
        flash('No selected file', 'error')
        return redirect(request.url)
    
    # Get pipeline options
    skip_qc = 'skip_qc' in request.form
    skip_trim = 'skip_trim' in request.form
    
    if r1_file and r2_file and allowed_file(r1_file.filename) and allowed_file(r2_file.filename):
        # Create a unique ID for this analysis
        analysis_id = str(uuid.uuid4())
        
        # Create folder for this analysis
        analysis_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], analysis_id)
        os.makedirs(analysis_folder, exist_ok=True)
        
        # Save the files
        r1_path = os.path.join(analysis_folder, 'r1.fastq.gz')
        r2_path = os.path.join(analysis_folder, 'r2.fastq.gz')
        r1_file.save(r1_path)
        r2_file.save(r2_path)
        
        # Start the pipeline
        result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
        os.makedirs(result_folder, exist_ok=True)
        
        # Run the pipeline in a background thread
        run_pipeline(analysis_id, r1_path, r2_path, result_folder, skip_qc, skip_trim)
        
        return redirect(url_for('main.analysis_status', analysis_id=analysis_id))
    
    flash('File type not allowed', 'error')
    return redirect(url_for('main.index'))

@main_bp.route('/upload_direct', methods=['POST'])
def upload_direct():
    """Handle direct file path uploads"""
    r1_path = request.form.get('r1_path')
    r2_path = request.form.get('r2_path')
    
    if not r1_path or not r2_path:
        flash('Both file paths are required', 'error')
        return redirect(url_for('main.index'))
    
    # Check if files exist
    if not os.path.exists(r1_path):
        flash(f'R1 file not found at: {r1_path}', 'error')
        return redirect(url_for('main.index'))
    
    if not os.path.exists(r2_path):
        flash(f'R2 file not found at: {r2_path}', 'error')
        return redirect(url_for('main.index'))
    
    # Check file types
    if not (allowed_file(r1_path) and allowed_file(r2_path)):
        flash('File type not allowed. Must be FASTQ format (.fastq, .fq, .fastq.gz, .fq.gz)', 'error')
        return redirect(url_for('main.index'))
    
    # Create a unique ID for this analysis
    analysis_id = str(uuid.uuid4())
    
    # Create folder for this analysis
    analysis_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], analysis_id)
    os.makedirs(analysis_folder, exist_ok=True)
    
    # Copy the files to the upload folder
    dest_r1_path = os.path.join(analysis_folder, 'r1.fastq.gz')
    dest_r2_path = os.path.join(analysis_folder, 'r2.fastq.gz')
    
    shutil.copy2(r1_path, dest_r1_path)
    shutil.copy2(r2_path, dest_r2_path)
    
    # Start the pipeline
    result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
    os.makedirs(result_folder, exist_ok=True)
    
    # Run the pipeline in a background thread
    run_pipeline(analysis_id, dest_r1_path, dest_r2_path, result_folder)
    
    return redirect(url_for('main.analysis_status', analysis_id=analysis_id))

@main_bp.route('/upload_bam', methods=['POST'])
def upload_bam():
    """Handle BAM file upload."""
    if 'bam_file' not in request.files:
        flash('No file part', 'error')
        return redirect(url_for('main.index'))
    
    bam_file = request.files['bam_file']
    
    if bam_file.filename == '':
        flash('No selected file', 'error')
        return redirect(url_for('main.index'))
    
    if bam_file and bam_file.filename.endswith('.bam'):
        # Create a unique ID for this analysis
        analysis_id = str(uuid.uuid4())
        
        # Create folder for this analysis
        analysis_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], analysis_id)
        os.makedirs(analysis_folder, exist_ok=True)
        
        # Save the file
        bam_path = os.path.join(analysis_folder, 'sample.bam')
        bam_file.save(bam_path)
        
        # Start the pipeline
        result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
        os.makedirs(result_folder, exist_ok=True)
        
        # Run the pipeline in a background thread, starting from the BAM file
        run_pipeline(analysis_id, bam_path=bam_path, result_folder=result_folder)
        
        return redirect(url_for('main.analysis_status', analysis_id=analysis_id))
    
    flash('File type not allowed. Must be BAM format (.bam)', 'error')
    return redirect(url_for('main.index'))

@main_bp.route('/upload_vcf', methods=['POST'])
def upload_vcf():
    """Handle VCF file upload."""
    if 'vcf_file' not in request.files:
        flash('No file part', 'error')
        return redirect(url_for('main.index'))
    
    vcf_file = request.files['vcf_file']
    
    if vcf_file.filename == '':
        flash('No selected file', 'error')
        return redirect(url_for('main.index'))
    
    if vcf_file and vcf_file.filename.endswith('.vcf'):
        # Create a unique ID for this analysis
        analysis_id = str(uuid.uuid4())
        
        # Create folder for this analysis
        analysis_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], analysis_id)
        os.makedirs(analysis_folder, exist_ok=True)
        
        # Save the file
        vcf_path = os.path.join(analysis_folder, 'variants.vcf')
        vcf_file.save(vcf_path)
        
        # Start the pipeline
        result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
        os.makedirs(result_folder, exist_ok=True)
        
        # Run the pipeline in a background thread, starting with the VCF file
        run_pipeline(analysis_id, vcf_path=vcf_path, result_folder=result_folder)
        
        return redirect(url_for('main.analysis_status', analysis_id=analysis_id))
    
    flash('File type not allowed. Must be VCF format (.vcf)', 'error')
    return redirect(url_for('main.index'))

@main_bp.route('/status/<analysis_id>')
def analysis_status(analysis_id):
    # Check if results are ready
    result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
    vcf_file = os.path.join(result_folder, 'variants.vcf')
    
    if os.path.exists(vcf_file):
        status = 'complete'
    else:
        status = 'processing'
    
    return render_template('status.html', analysis_id=analysis_id, status=status)

@main_bp.route('/api/status/<analysis_id>')
def api_status(analysis_id):
    # Check if results are ready
    result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
    vcf_file = os.path.join(result_folder, 'variants.vcf')
    log_file = os.path.join(result_folder, 'pipeline.log')
    progress_file = os.path.join(result_folder, 'progress.json')
    
    if os.path.exists(vcf_file):
        status = 'complete'
    else:
        status = 'processing'
        
    # Get log contents if available
    log_content = ""
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            log_content = f.read()
    
    # Get progress information if available
    progress = {}
    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            progress = json.load(f)
    
    return jsonify({
        'status': status,
        'log': log_content,
        'progress': progress
    })

@main_bp.route('/download/<analysis_id>/<path:filename>')
def download_file(analysis_id, filename):
    """Download a specific file from the results directory."""
    result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
    
    # Special case for FastQC and fastp reports to ensure they exist
    if filename.startswith('fastqc/'):
        if not os.path.exists(os.path.join(result_folder, filename)):
            # Check if the directory exists
            fastqc_dir = os.path.join(result_folder, 'fastqc')
            if os.path.isdir(fastqc_dir):
                # List all HTML files in the fastqc directory
                html_files = [f for f in os.listdir(fastqc_dir) if f.endswith('.html')]
                if html_files:
                    if 'r1' in filename.lower():
                        # Find a suitable R1 report
                        r1_files = [f for f in html_files if 'r1' in f.lower() or '_1_' in f.lower() or '_1.' in f.lower()]
                        if r1_files:
                            return send_from_directory(fastqc_dir, r1_files[0])
                    elif 'r2' in filename.lower():
                        # Find a suitable R2 report
                        r2_files = [f for f in html_files if 'r2' in f.lower() or '_2_' in f.lower() or '_2.' in f.lower()]
                        if r2_files:
                            return send_from_directory(fastqc_dir, r2_files[0])
                    # Default to the first HTML file if no match
                    return send_from_directory(fastqc_dir, html_files[0])
                else:
                    # No HTML files found, show a message
                    return f"<html><body><h1>FastQC reports not found</h1><p>No FastQC HTML reports were found in the results folder.</p></body></html>"
    
    elif filename.startswith('trimmed/'):
        if not os.path.exists(os.path.join(result_folder, filename)):
            # Check if the directory exists
            trimmed_dir = os.path.join(result_folder, 'trimmed')
            if os.path.isdir(trimmed_dir):
                # List all HTML files in the trimmed directory
                html_files = [f for f in os.listdir(trimmed_dir) if f.endswith('.html')]
                if html_files:
                    # Return the first HTML file (likely the fastp report)
                    return send_from_directory(trimmed_dir, html_files[0])
                else:
                    # No HTML files found, show a message
                    return f"<html><body><h1>Trimming report not found</h1><p>No fastp HTML report was found in the results folder.</p></body></html>"
    
    # Handle directory requests (like /fastqc)
    if os.path.isdir(os.path.join(result_folder, filename)):
        files = os.listdir(os.path.join(result_folder, filename))
        html_files = [f for f in files if f.endswith('.html')]
        
        if html_files:
            # Redirect to the first HTML file
            return redirect(url_for('main.download_file', analysis_id=analysis_id, filename=os.path.join(filename, html_files[0])))
        else:
            # Display a list of files
            html = '<html><head><title>Directory Listing</title></head><body>'
            html += f'<h1>Files in {filename}</h1><ul>'
            for f in files:
                html += f'<li><a href="{url_for("main.download_file", analysis_id=analysis_id, filename=os.path.join(filename, f))}" target="_blank">{f}</a></li>'
            html += '</ul></body></html>'
            return html
    
    # Check if file exists
    if not os.path.exists(os.path.join(result_folder, filename)):
        return f"<html><body><h1>File not found</h1><p>The requested file {filename} was not found in the results folder.</p></body></html>"
        
    # Send the file
    return send_from_directory(result_folder, filename)

@main_bp.route('/visualize_bam', methods=['POST'])
def visualize_bam():
    """Handle BAM file visualization with IGV."""
    if 'bam_file' not in request.files:
        flash('No BAM file part', 'error')
        return redirect(url_for('main.index'))
    
    bam_file = request.files['bam_file']
    
    if bam_file.filename == '':
        flash('No selected BAM file', 'error')
        return redirect(url_for('main.index'))
    
    if not bam_file.filename.endswith('.bam'):
        flash('File type not allowed. Must be BAM format (.bam)', 'error')
        return redirect(url_for('main.index'))
    
    # Create a unique ID for this visualization
    visualization_id = str(uuid.uuid4())
    
    # Create folder for this visualization
    visualization_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], visualization_id)
    os.makedirs(visualization_folder, exist_ok=True)
    
    # Save the BAM file
    bam_path = os.path.join(visualization_folder, 'sample.bam')
    bam_file.save(bam_path)
    
    # Check if the BAM file is valid
    try:
        # First check if the file is a valid BAM using samtools quickcheck
        check_cmd = f"samtools quickcheck -v {bam_path} 2>&1"
        result = subprocess.run(check_cmd, shell=True, check=False, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_message = result.stdout or result.stderr or "Unknown error"
            if "EOF" in error_message:
                flash(f'The BAM file appears to be truncated or corrupted (missing EOF block). This often happens when file uploads are interrupted. Please try the following solutions: 1) Upload a fresh copy of the BAM file, 2) Try a smaller BAM file, or 3) Fix the BAM file using "samtools reheader" locally.', 'error')
            else:
                flash(f'Invalid BAM file: {error_message}', 'error')
            # Clean up the invalid file
            shutil.rmtree(visualization_folder)
            return redirect(url_for('main.index'))
            
        # Check if the BAM file actually contains reads
        count_cmd = f"samtools view -c {bam_path}"
        count_result = subprocess.run(count_cmd, shell=True, check=True, capture_output=True, text=True)
        read_count = int(count_result.stdout.strip())
        
        if read_count == 0:
            flash('The BAM file does not contain any reads. Please upload a BAM file with aligned reads.', 'error')
            shutil.rmtree(visualization_folder)
            return redirect(url_for('main.index'))
            
    except Exception as e:
        flash(f'Error checking BAM file: {str(e)}', 'error')
        # Clean up in case of error
        if os.path.exists(visualization_folder):
            shutil.rmtree(visualization_folder)
        return redirect(url_for('main.index'))
    
    # Handle BAI file (index)
    bai_path = os.path.join(visualization_folder, 'sample.bam.bai')
    bai_created = False
    
    # Check if user uploaded a BAI file
    if 'bai_file' in request.files and request.files['bai_file'].filename != '':
        bai_file = request.files['bai_file']
        if bai_file.filename.endswith('.bai'):
            bai_file.save(bai_path)
            
            # Verify the BAI file is valid
            try:
                # A simple test to see if the BAI file can be used to fetch a region
                test_cmd = f"samtools view {bam_path} chr1:1-1000 -c"
                subprocess.run(test_cmd, shell=True, check=True, capture_output=True)
                bai_created = True
            except:
                flash('The provided BAI index file appears to be invalid or incompatible with the BAM file.', 'error')
                # Remove the invalid BAI
                if os.path.exists(bai_path):
                    os.remove(bai_path)
    
    # If BAI wasn't uploaded or was invalid, try to create one
    if not bai_created:
        try:
            # First check if the BAM is sorted
            sort_check_cmd = f"samtools view -H {bam_path} | grep 'SO:coordinate'"
            sort_result = subprocess.run(sort_check_cmd, shell=True, check=False, capture_output=True)
            
            if sort_result.returncode != 0 or not sort_result.stdout:
                # BAM is not sorted
                flash('BAM file is not sorted by coordinate. Sorting now (this may take a while for large files)...', 'info')
                
                # Sort the BAM file
                sorted_bam = os.path.join(visualization_folder, 'sample.sorted.bam')
                sort_cmd = f"samtools sort {bam_path} -o {sorted_bam}"
                subprocess.run(sort_cmd, shell=True, check=True)
                
                # Replace original with sorted version
                os.remove(bam_path)
                os.rename(sorted_bam, bam_path)
                
                flash('BAM file has been sorted successfully.', 'success')
            
            # Create BAI index
            flash('Creating BAM index file...', 'info')
            index_cmd = f"samtools index {bam_path}"
            subprocess.run(index_cmd, shell=True, check=True)
            
            # Verify the index was created
            if os.path.exists(bai_path) and os.path.getsize(bai_path) > 0:
                bai_created = True
                flash('BAM index created successfully.', 'success')
            else:
                flash('Failed to create BAM index file.', 'error')
                
        except subprocess.CalledProcessError as e:
            flash(f'Error during BAM processing: {str(e)}. The file may be corrupted or in an unsupported format.', 'error')
            if os.path.exists(visualization_folder):
                shutil.rmtree(visualization_folder)
            return redirect(url_for('main.index'))
    
    # Final check to ensure we have a valid index
    if not bai_created or not os.path.exists(bai_path) or os.path.getsize(bai_path) == 0:
        flash('Could not create or find a valid BAM index file (.bai). The BAM file might be corrupted.', 'error')
        if os.path.exists(visualization_folder):
            shutil.rmtree(visualization_folder)
        return redirect(url_for('main.index'))
    
    # Prepare IGV.js launch URL
    igv_url = url_for('main.view_igv', visualization_id=visualization_id)
    
    return redirect(igv_url)

@main_bp.route('/igv/<visualization_id>')
def view_igv(visualization_id):
    """Display the IGV.js browser for a BAM file."""
    bam_url = url_for('main.get_bam_file', visualization_id=visualization_id, filename='sample.bam')
    bai_url = url_for('main.get_bam_file', visualization_id=visualization_id, filename='sample.bam.bai')
    
    return render_template('igv.html', 
                          visualization_id=visualization_id,
                          bam_url=bam_url,
                          bai_url=bai_url)

@main_bp.route('/igv/data/<visualization_id>/<path:filename>')
def get_bam_file(visualization_id, filename):
    """Serve BAM and BAI files for IGV."""
    visualization_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], visualization_id)
    return send_from_directory(visualization_folder, filename)

@main_bp.route('/continue_analysis', methods=['POST'])
def continue_analysis():
    """Handle continuing an existing analysis."""
    analysis_id = request.form.get('analysis_id')
    continue_step = request.form.get('continue_step', '')  # Default to auto-detect if not specified
    
    if not analysis_id:
        flash('Analysis ID is required', 'error')
        return redirect(url_for('main.index'))
    
    # Check if this analysis exists
    result_folder = os.path.join(current_app.config['RESULTS_FOLDER'], analysis_id)
    if not os.path.exists(result_folder):
        flash(f'Analysis {analysis_id} not found', 'error')
        return redirect(url_for('main.index'))
    
    # Valid steps for continuation
    valid_steps = ["sam2bam", "sort", "dedup", "bqsr", "variants", "filter", "annotate", "compress"]
    
    if continue_step and continue_step not in valid_steps:
        flash(f'Invalid step: {continue_step}', 'error')
        return redirect(url_for('main.index'))
    
    # Auto-detect step if not specified
    if not continue_step:
        # Load progress file
        progress_file = os.path.join(result_folder, 'progress.json')
        if os.path.exists(progress_file):
            with open(progress_file, 'r') as f:
                progress = json.load(f)
            
            # Find the first incomplete step
            for step in progress.get('steps', []):
                if step['status'] in ['pending', 'error']:
                    continue_step = step['id']
                    break
        
        if not continue_step:
            flash('Could not auto-detect next step. Please select a specific step.', 'error')
            return redirect(url_for('main.index'))
    
    # Start continuing the pipeline
    from app.pipeline.pipeline import continue_pipeline
    continue_pipeline(analysis_id, continue_step)
    
    return redirect(url_for('main.analysis_status', analysis_id=analysis_id)) 