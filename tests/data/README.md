# Test data for the metagenomics pipeline

This directory holds tiny synthetic FASTQ files used for pipeline integration
testing (dry-run validation).

## Generating minimal test reads

```bash
# Create tiny paired-end FASTQ files (10 reads each)
python - <<'EOF'
import gzip, random, string, pathlib

pathlib.Path("data/reads").mkdir(parents=True, exist_ok=True)

def random_read(length=150):
    bases = "ACGT"
    seq  = "".join(random.choices(bases, k=length))
    qual = "I" * length          # phred 40 for all bases
    return seq, qual

for sample in ["sample1", "sample2"]:
    for mate in ["R1", "R2"]:
        fname = f"data/reads/{sample}_{mate}.fastq.gz"
        with gzip.open(fname, "wt") as fh:
            for i in range(10):
                seq, qual = random_read()
                fh.write(f"@{sample}_{mate}_read{i}\n{seq}\n+\n{qual}\n")

print("Test reads generated.")
EOF
```
