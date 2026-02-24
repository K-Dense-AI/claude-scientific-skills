# Genomic NRM (0-2 Scale)

## Input -> Process -> Output

### Input
- Simulated marker matrix (0/1/2 dosage)
- Marker allele frequencies

### Process
1. Center genotypes by allele frequency
2. Compute genomic relationship matrix on 0-2 scale
3. Export matrix and heatmap

### Output
- `output/genomic_nrm_0_2.csv`
- `output/genomic_nrm_heatmap.png`

## Run
```bash
python run_genomic_nrm.py
```
