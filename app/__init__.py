#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from flask import Flask
from werkzeug.utils import secure_filename

def create_app(test_config=None):
    # Create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    
    app.config.from_mapping(
        SECRET_KEY='dev',
        UPLOAD_FOLDER=os.path.join(os.getcwd(), 'uploads'),
        RESULTS_FOLDER=os.path.join(os.getcwd(), 'results'),
        MAX_CONTENT_LENGTH=16 * 1024 * 1024 * 1024,  # 16GB max upload size
        ALLOWED_EXTENSIONS={'fastq', 'fq', 'fastq.gz', 'fq.gz', 'bam', 'vcf'}
    )

    if test_config is None:
        # Load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # Load the test config if passed in
        app.config.from_mapping(test_config)

    # Ensure the upload and results folders exist
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

    # Register blueprints
    from app.routes import main_bp
    app.register_blueprint(main_bp)

    return app 