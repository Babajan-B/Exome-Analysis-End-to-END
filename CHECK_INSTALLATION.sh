#!/bin/bash
# Quick Installation Check Script
# Verifies all tools are installed and in PATH

echo "╔════════════════════════════════════════════════════════════╗"
echo "║           CHECKING INSTALLATION STATUS                     ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Function to check command
check_cmd() {
    local cmd=$1
    local path=$2
    
    if command -v $cmd &> /dev/null; then
        echo "✅ $cmd: $(which $cmd)"
        return 0
    elif [ -f "$path" ]; then
        echo "⚠️  $cmd: Found at $path but NOT in PATH"
        return 1
    else
        echo "❌ $cmd: NOT FOUND"
        return 2
    fi
}

# Check tools
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Checking Bioinformatics Tools:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

check_cmd "fastqc" "~/NGS/tools/FastQC/fastqc"
check_cmd "fastp" "~/NGS/tools/fastp"
check_cmd "bwa" "/usr/bin/bwa"
check_cmd "samtools" "/usr/bin/samtools"
check_cmd "bgzip" "/usr/bin/bgzip"
check_cmd "tabix" "/usr/bin/tabix"

# Check GATK
echo ""
if [ -f "/opt/gatk-4.6.2.0/gatk" ]; then
    echo "✅ GATK: /opt/gatk-4.6.2.0/gatk"
elif [ -f "~/NGS/tools/gatk-4.6.2.0/gatk" ]; then
    echo "⚠️  GATK: Found at ~/NGS/tools/gatk-4.6.2.0/gatk but NOT in /opt/"
else
    echo "❌ GATK: NOT FOUND"
fi

# Check ANNOVAR
echo ""
if [ -d "~/NGS/tools/annovar" ]; then
    echo "✅ ANNOVAR: ~/NGS/tools/annovar"
else
    echo "❌ ANNOVAR: NOT FOUND"
fi

# Check reference genome
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Checking Reference Genome:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "~/NGS/reference/hg19.fa" ]; then
    echo "✅ hg19.fa: ~/NGS/reference/hg19.fa"
    
    if [ -f "~/NGS/reference/hg19.fa.bwt" ]; then
        echo "✅ BWA index: Present"
    else
        echo "❌ BWA index: Missing"
    fi
    
    if [ -f "~/NGS/reference/hg19.fa.fai" ]; then
        echo "✅ FASTA index: Present"
    else
        echo "❌ FASTA index: Missing"
    fi
    
    if [ -f "~/NGS/reference/hg19.dict" ]; then
        echo "✅ GATK dict: Present"
    else
        echo "❌ GATK dict: Missing"
    fi
else
    echo "❌ hg19.fa: NOT FOUND"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Storage Usage:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
du -sh ~/NGS 2>/dev/null || echo "~/NGS directory not found"
df -h ~ | tail -1

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Recommendation:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "If tools are missing or not in PATH, run:"
echo "  bash install_all.sh"
echo ""
echo "Or if tools exist but not in PATH, add them:"
echo "  export PATH=\$PATH:~/NGS/tools/FastQC:~/NGS/tools"
echo ""

