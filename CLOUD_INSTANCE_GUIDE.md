# Cloud Instance Setup Guide

Step-by-step instructions to run the NGS Exome Analysis Pipeline on any shared cloud server (AWS, GCP, Azure, on-prem, etc.).

---

## 🎯 Quick Reference

1. Launch a Linux VM (Ubuntu 20.04+ recommended).
2. Minimum configuration: 16 vCPU, 32 GB RAM, 100 GB SSD.
3. SSH into the instance and run the commands below.

### Command 1 – install everything (run once)
```bash
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git ~/NGS
cd ~/NGS
bash install_all.sh
```

### Command 2 – upload data
```bash
mkdir -p ~/NGS/data
# copy your *R1.fastq.gz and *R2.fastq.gz files into ~/NGS/data/
```

### Command 3 – run the ultimate pipeline
```bash
bash ULTIMATE_MASTER_PIPELINE.sh ~/NGS/data 16
```

### Download results
- Archive generated at: `~/NGS/NGS_Results_Complete_<timestamp>.zip`
- Includes annotated TXT files, compressed VCFs, HTML reports, and summary.

---

## 🚀 Instance Selection Tips

| Use case | Recommended instance |
|----------|----------------------|
| Test / tutorial | 8 vCPU, 16 GB RAM |
| Typical exome | 16 vCPU, 32 GB RAM |
| Multiple samples | 32 vCPU, 64 GB RAM |

- Make sure the disk has at least 100 GB free space.
- Enable outbound internet access to download tools and reference data.
- Close the VM when finished to avoid charges.

---

## 🔐 Security & Access

- Use SSH key authentication when possible.
- Restrict inbound firewall rules to your IP.
- Store sensitive FASTQ/VCF files in encrypted storage if required by policy.

---

## 🛠 Common Tasks

**Check disk usage**
```bash
df -h ~
```

**Monitor pipeline log**
```bash
tail -f ~/NGS/results/<sample>/pipeline.log
```

**Resume after interruption**
- Re-run `ULTIMATE_MASTER_PIPELINE.sh` – it skips completed steps automatically.

---

## ✅ Tear-down checklist

- Download the ZIP archive with results.
- Delete intermediate files if needed (`rm -rf ~/NGS/results`).
- Shut down or delete the cloud instance to stop billing.

---

Need a laptop-specific guide instead? See `INSTALL_MAC.md`, `INSTALL_LINUX.md`, or `INSTALL_WINDOWS.md` in the `docs/` folder.
