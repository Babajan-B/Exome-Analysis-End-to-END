# Installation Guide for Linux

This guide covers installation on Ubuntu/Debian-based systems. For other distributions, adjust package manager commands accordingly.

## 📋 Prerequisites

- Ubuntu 20.04+ or Debian 10+ (or equivalent Linux distribution)
- Sudo/root access
- At least 50GB free disk space
- Internet connection

## 🚀 Step-by-Step Installation

### Step 1: Update System Packages

```bash
# Update package lists
sudo apt update

# Upgrade existing packages
sudo apt upgrade -y
```

### Step 2: Install System Dependencies

```bash
# Install build essentials and tools
sudo apt install -y build-essential wget curl git unzip

# Install Python 3 and pip
sudo apt install -y python3 python3-pip python3-venv

# Install Java (required for GATK)
sudo apt install -y default-jre default-jdk

# Verify installations
python3 --version
java -version
git --version
```

### Step 3: Install Bioinformatics Tools

#### FastQC
```bash
sudo apt install -y fastqc
fastqc --version
```

#### fastp
```bash
# Download and install fastp
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
sudo mv ./fastp /usr/local/bin/
fastp --version
```

#### BWA (Burrows-Wheeler Aligner)
```bash
sudo apt install -y bwa
bwa
```

#### Samtools
```bash
sudo apt install -y samtools
samtools --version
```

#### BCFtools (optional but useful)
```bash
sudo apt install -y bcftools
bcftools --version
```

#### Tabix
```bash
sudo apt install -y tabix
tabix --version
```

### Step 4: Install GATK

```bash
# Download GATK
cd /tmp
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip

# Extract
unzip gatk-4.6.2.0.zip

# Move to /opt
sudo mv gatk-4.6.2.0 /opt/

# Create symbolic link
sudo ln -s /opt/gatk-4.6.2.0/gatk /usr/local/bin/gatk

# Test GATK
gatk --version
```

### Step 5: Clone the Repository

```bash
# Navigate to desired installation directory
cd ~

# Clone the repository
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git

# Navigate to project directory
cd Exome-Analysis-End-to-END
```

### Step 6: Set Up Python Virtual Environment

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install Python dependencies
pip install -r requirements.txt
```

### Step 7: Download Reference Genome

```bash
# Create reference directory
mkdir -p reference
cd reference

# Download hg19 reference genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# Extract
gunzip hg19.fa.gz

# Index for BWA
bwa index hg19.fa

# Create FASTA index
samtools faidx hg19.fa

# Create sequence dictionary for GATK
gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict

cd ..
```

**For hg38 (newer reference):**
```bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bwa index hg38.fa
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
```

### Step 8: Download Known Variants (Optional but Recommended)

```bash
# Create known sites directory
mkdir -p reference/known_sites
cd reference/known_sites

# Download dbSNP for hg19
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi

# Rename
mv 00-All.vcf.gz dbsnp_138.hg19.vcf.gz
mv 00-All.vcf.gz.tbi dbsnp_138.hg19.vcf.gz.tbi

cd ../..
```

### Step 9: Run Automated Setup (Optional)

```bash
# Make setup script executable
chmod +x setup.sh

# Run setup script
./setup.sh
```

### Step 10: Start the Application

```bash
# Activate virtual environment (if not already activated)
source venv/bin/activate

# Start Flask server
python run.py
```

### Step 11: Access Web Interface

Open your browser and navigate to:
```
http://localhost:5008
```

Or if running on a remote server:
```
http://your-server-ip:5008
```

## 🐧 Distribution-Specific Instructions

### CentOS/RHEL/Fedora

```bash
# Update system
sudo yum update -y  # or dnf update -y

# Install dependencies
sudo yum install -y gcc gcc-c++ make wget curl git unzip python3 python3-pip java-11-openjdk

# Install bioinformatics tools
sudo yum install -y fastqc bwa samtools

# For fastp and GATK, follow manual installation steps above
```

### Arch Linux

```bash
# Update system
sudo pacman -Syu

# Install dependencies
sudo pacman -S base-devel wget curl git unzip python python-pip jdk-openjdk

# Install bioinformatics tools
sudo pacman -S fastqc bwa samtools

# Install from AUR
yay -S fastp gatk
```

## 🐳 Docker Installation (Alternative)

If you prefer using Docker:

```bash
# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker

# Build Docker image
cd Exome-Analysis-End-to-END
docker build -t ngs-exome-analysis .

# Run container
docker run -d -p 5008:5008 \
    -v $(pwd)/uploads:/app/uploads \
    -v $(pwd)/results:/app/results \
    ngs-exome-analysis
```

## 🖥️ Running as a Service (Production)

### Using systemd

Create a service file:

```bash
sudo nano /etc/systemd/system/ngs-exome.service
```

Add the following content:

```ini
[Unit]
Description=NGS Exome Analysis Pipeline
After=network.target

[Service]
Type=simple
User=your-username
WorkingDirectory=/home/your-username/Exome-Analysis-End-to-END
Environment="PATH=/home/your-username/Exome-Analysis-End-to-END/venv/bin"
ExecStart=/home/your-username/Exome-Analysis-End-to-END/venv/bin/gunicorn -w 4 -b 0.0.0.0:5008 run:app
Restart=always

[Install]
WantedBy=multi-user.target
```

Enable and start the service:

```bash
# Reload systemd
sudo systemctl daemon-reload

# Enable service
sudo systemctl enable ngs-exome.service

# Start service
sudo systemctl start ngs-exome.service

# Check status
sudo systemctl status ngs-exome.service
```

### Using Nginx as Reverse Proxy

```bash
# Install Nginx
sudo apt install -y nginx

# Create Nginx configuration
sudo nano /etc/nginx/sites-available/ngs-exome
```

Add configuration:

```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://127.0.0.1:5008;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        client_max_body_size 20G;
    }
}
```

Enable and restart:

```bash
sudo ln -s /etc/nginx/sites-available/ngs-exome /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

## ⚡ Performance Optimization

### Increase File Descriptor Limits

```bash
# Edit limits
sudo nano /etc/security/limits.conf

# Add lines:
* soft nofile 65536
* hard nofile 65536

# Apply changes
sudo sysctl -p
```

### Optimize for Large Files

```bash
# Increase inotify limits
echo fs.inotify.max_user_watches=524288 | sudo tee -a /etc/sysctl.conf
sudo sysctl -p
```

### Use tmpfs for Temporary Files (if you have enough RAM)

```bash
# Edit /etc/fstab
sudo nano /etc/fstab

# Add line:
tmpfs /tmp tmpfs defaults,size=8G 0 0

# Mount
sudo mount -a
```

## ❓ Troubleshooting

### Issue: Permission denied

```bash
# Fix permissions on project directory
chmod -R u+rwx ~/Exome-Analysis-End-to-END

# Fix virtual environment
chmod +x venv/bin/activate
```

### Issue: Port 5008 already in use

```bash
# Find process using port
sudo lsof -i :5008

# Kill process
sudo kill -9 <PID>

# Or change port in run.py
```

### Issue: Out of memory

```bash
# Check memory usage
free -h

# Monitor during analysis
watch -n 1 free -h

# Reduce threads in pipeline or add swap
sudo fallocate -l 8G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

### Issue: Tool not found

```bash
# Check if tool is installed
which fastqc
which bwa
which samtools
which gatk

# Check PATH
echo $PATH

# Add to PATH if needed
export PATH="/usr/local/bin:$PATH"
echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.bashrc
```

### Issue: Java heap space error (GATK)

```bash
# Increase Java heap size
export _JAVA_OPTIONS="-Xmx8g"

# Or edit GATK command in pipeline.py
```

## 🔒 Security Considerations

### Firewall Configuration

```bash
# Allow port 5008
sudo ufw allow 5008/tcp

# Enable firewall
sudo ufw enable

# Check status
sudo ufw status
```

### SSL/TLS with Let's Encrypt

```bash
# Install certbot
sudo apt install -y certbot python3-certbot-nginx

# Get certificate
sudo certbot --nginx -d your-domain.com

# Auto-renewal
sudo certbot renew --dry-run
```

## 📊 Monitoring

### View Logs

```bash
# Application logs
tail -f results/*/pipeline.log

# System service logs
sudo journalctl -u ngs-exome.service -f

# Nginx logs
sudo tail -f /var/log/nginx/access.log
```

### Resource Monitoring

```bash
# Install htop
sudo apt install -y htop

# Monitor resources
htop

# Or use top
top
```

## 🧪 Test Installation

```bash
# Quick test
curl http://localhost:5008

# Should return HTML of the homepage
```

## 📚 Additional Resources

- [Ubuntu Documentation](https://help.ubuntu.com/)
- [GATK Documentation](https://gatk.broadinstitute.org/)
- [BWA Manual](http://bio-bwa.sourceforge.net/)
- [Samtools Documentation](http://www.htslib.org/doc/)

## 🆘 Getting Help

If you encounter issues:
1. Check logs: `results/[analysis_id]/pipeline.log`
2. Verify all tools are installed: `which fastqc bwa samtools gatk`
3. Check system resources: `free -h`, `df -h`
4. Open an issue on [GitHub](https://github.com/Babajan-B/Exome-Analysis-End-to-END/issues)

---

**Linux Installation Complete!** 🎉

Access the application at: http://localhost:5008

