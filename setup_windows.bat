@echo off
REM NGS Exome Analysis Pipeline - Windows Setup Script
REM This script sets up the Python virtual environment and installs dependencies

echo ========================================
echo NGS Exome Analysis Pipeline Setup
echo ========================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is not installed or not in PATH
    echo Please install Python 3.7+ from https://www.python.org/downloads/
    echo Make sure to check "Add Python to PATH" during installation
    pause
    exit /b 1
)

echo [1/5] Python found
python --version

REM Check if Java is installed
java -version >nul 2>&1
if errorlevel 1 (
    echo WARNING: Java is not installed or not in PATH
    echo Java is required for GATK. Please install from:
    echo https://www.oracle.com/java/technologies/downloads/
    echo.
    echo You can continue setup and install Java later.
    pause
) else (
    echo [2/5] Java found
    java -version
)

REM Create virtual environment
echo.
echo [3/5] Creating Python virtual environment...
if exist venv (
    echo Virtual environment already exists. Skipping creation.
) else (
    python -m venv venv
    if errorlevel 1 (
        echo ERROR: Failed to create virtual environment
        pause
        exit /b 1
    )
    echo Virtual environment created successfully
)

REM Activate virtual environment
echo.
echo [4/5] Activating virtual environment...
call venv\Scripts\activate.bat
if errorlevel 1 (
    echo ERROR: Failed to activate virtual environment
    pause
    exit /b 1
)

REM Upgrade pip
echo.
echo [5/5] Installing Python dependencies...
python -m pip install --upgrade pip

REM Install requirements
pip install -r requirements.txt
if errorlevel 1 (
    echo ERROR: Failed to install Python dependencies
    pause
    exit /b 1
)

REM Create necessary directories
echo.
echo Creating project directories...
if not exist uploads mkdir uploads
if not exist results mkdir results
if not exist reference mkdir reference
if not exist test_data mkdir test_data

echo.
echo ========================================
echo Setup completed successfully!
echo ========================================
echo.
echo IMPORTANT NOTES FOR WINDOWS USERS:
echo.
echo 1. We recommend using WSL2 (Windows Subsystem for Linux) for better
echo    compatibility with bioinformatics tools.
echo.
echo 2. To install WSL2, run this command in PowerShell as Administrator:
echo    wsl --install
echo.
echo 3. After WSL2 is installed, you can run the Linux version of this pipeline
echo    which provides better performance and compatibility.
echo.
echo 4. To start the application:
echo    - Open Command Prompt or PowerShell
echo    - Navigate to this directory
echo    - Run: venv\Scripts\activate
echo    - Run: python run.py
echo    - Open browser to: http://localhost:5008
echo.
echo For detailed instructions, see: docs\INSTALL_WINDOWS.md
echo.
pause

