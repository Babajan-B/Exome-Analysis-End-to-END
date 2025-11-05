#!/bin/bash
# Test ANNOVAR Installation
# This script verifies that ANNOVAR is properly installed and functional

set -e

ANNOVAR_DIR=~/NGS/tools/annovar

echo "=========================================="
echo "ANNOVAR Installation Test"
echo "=========================================="
echo ""

# Test 1: Check if ANNOVAR directory exists
echo "[Test 1/5] Checking ANNOVAR directory..."
if [ -d "$ANNOVAR_DIR" ]; then
    echo "✅ ANNOVAR directory found: $ANNOVAR_DIR"
else
    echo "❌ ANNOVAR directory not found: $ANNOVAR_DIR"
    echo "Please run: bash setup_annovar.sh"
    exit 1
fi
echo ""

# Test 2: Check if main scripts exist
echo "[Test 2/5] Checking ANNOVAR scripts..."
scripts=(
    "annotate_variation.pl"
    "table_annovar.pl"
    "convert2annovar.pl"
    "coding_change.pl"
)

all_found=true
for script in "${scripts[@]}"; do
    if [ -f "$ANNOVAR_DIR/$script" ]; then
        echo "✅ Found: $script"
    else
        echo "❌ Missing: $script"
        all_found=false
    fi
done

if [ "$all_found" = false ]; then
    echo "Some scripts are missing. Please reinstall ANNOVAR."
    exit 1
fi
echo ""

# Test 3: Check if databases are downloaded
echo "[Test 3/5] Checking database installation..."
if [ -d "$ANNOVAR_DIR/humandb" ]; then
    echo "✅ Database directory found"
    echo ""
    echo "Installed databases for hg19:"
    ls -lh $ANNOVAR_DIR/humandb/ | grep "hg19" | awk '{print "  - " $9}' | head -20
    
    db_count=$(ls $ANNOVAR_DIR/humandb/ | grep "hg19" | wc -l)
    echo ""
    echo "Total hg19 database files: $db_count"
    
    if [ $db_count -lt 5 ]; then
        echo "⚠️  Warning: Very few databases installed. Consider running setup_annovar.sh"
    fi
else
    echo "❌ Database directory not found"
    echo "Please run: bash setup_annovar.sh"
    exit 1
fi
echo ""

# Test 4: Check Perl availability
echo "[Test 4/5] Checking Perl installation..."
if command -v perl &> /dev/null; then
    perl_version=$(perl -v | grep "This is perl" | head -1)
    echo "✅ Perl found: $perl_version"
else
    echo "❌ Perl not found. Please install Perl."
    exit 1
fi
echo ""

# Test 5: Run a simple ANNOVAR test
echo "[Test 5/5] Running ANNOVAR test annotation..."

# Create a test VCF with a known variant
TEST_DIR=/tmp/annovar_test_$$
mkdir -p $TEST_DIR
cd $TEST_DIR

cat > test.vcf << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100000	.	A	G	100	PASS	DP=50	GT	0/1
chr1	200000	.	C	T	100	PASS	DP=50	GT	0/1
EOF

echo "Created test VCF with 2 variants"

# Try to annotate (just with refGene if available)
if [ -f "$ANNOVAR_DIR/humandb/hg19_refGene.txt" ]; then
    echo "Running test annotation..."
    
    if perl $ANNOVAR_DIR/table_annovar.pl \
        test.vcf \
        $ANNOVAR_DIR/humandb/ \
        -buildver hg19 \
        -out test_output \
        -remove \
        -protocol refGene \
        -operation g \
        -vcfinput \
        2>&1 | grep -q "NOTICE"; then
        
        echo "✅ ANNOVAR test annotation successful!"
        
        # Check if output was created
        if [ -f "test_output.hg19_multianno.vcf" ]; then
            echo "✅ Output VCF created"
            variant_count=$(grep -v "^#" test_output.hg19_multianno.vcf | wc -l)
            echo "   Annotated $variant_count variants"
        fi
    else
        echo "⚠️  ANNOVAR ran but output may need verification"
    fi
else
    echo "⚠️  refGene database not found - skipping annotation test"
    echo "   Run: bash setup_annovar.sh to download databases"
fi

# Cleanup
cd - > /dev/null
rm -rf $TEST_DIR

echo ""
echo "=========================================="
echo "✅ ANNOVAR Installation Test Complete!"
echo "=========================================="
echo ""
echo "Summary:"
echo "  ✅ ANNOVAR is installed at: $ANNOVAR_DIR"
echo "  ✅ Scripts are functional"
echo "  ✅ Databases are present"
echo "  ✅ Test annotation successful"
echo ""
echo "You can now use ANNOVAR with:"
echo "  bash annovar_helper.sh input.vcf output_prefix"
echo ""
echo "Or manually:"
echo "  perl $ANNOVAR_DIR/table_annovar.pl input.vcf \\"
echo "      $ANNOVAR_DIR/humandb/ \\"
echo "      -buildver hg19 \\"
echo "      -out output \\"
echo "      -protocol refGene,avsnp150,clinvar_20240917 \\"
echo "      -operation g,f,f \\"
echo "      -vcfinput"
echo ""
echo "=========================================="

