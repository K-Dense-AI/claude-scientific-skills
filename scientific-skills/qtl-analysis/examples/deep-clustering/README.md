<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Deep Learning Population Clustering Example

## What is This?

**Deep learning clustering** uses neural networks (specifically autoencoders) to discover complex patterns in DNA data that traditional methods might miss. Think of it like teaching a computer to "see" patterns in genetic data that aren't obvious to standard statistical methods.

### Key Terms Explained

- **Autoencoder**: A type of neural network that learns to compress data (encoding) and then reconstruct it (decoding). The compressed version is called the "latent space."

- **Latent Space**: A compressed representation of the data. Instead of 200 SNPs, we have 10 numbers that capture the most important information.

- **Neural Network**: A computer model inspired by the brain, with layers of connected "neurons" that learn patterns from data.

- **t-SNE**: A technique to visualize high-dimensional data in 2D. It preserves relationships so similar samples stay close together.

- **ARI (Adjusted Rand Index)**: A score from 0 to 1 that measures how well the clustering matches the true labels. ARI = 1 means perfect match.

## Why Would Anyone Do This?

1. **Find Hidden Structure**: Traditional methods like PCA only find linear patterns. Neural networks can find non-linear patterns—think of it as finding curves, not just straight lines.

2. **Handle Mixed Populations**: When populations have mixed ancestry, simple methods struggle. Deep learning can separate subtle subpopulations.

3. **Reduce Dimensionality**: Convert 200,000 SNPs into 10 meaningful numbers that capture population structure.

4. **Better than K-means**: While K-means works in the original space, deep learning first transforms the data to make clustering easier.

## How It Works (The Process)

1. **Generate Data**: Create DNA data for 300 people from 4 populations, plus some admixed individuals.

2. **Train Autoencoder**: A neural network learns to compress the data (200 SNPs → 10 numbers) and reconstruct it.

3. **Extract Latent Space**: Get the compressed 10-number representation for each person.

4. **Cluster**: Use K-means to group people based on their latent space coordinates.

5. **Visualize**: Use t-SNE to show the 2D projection of the latent space.

## Input → Process → Output

### Input
| File | Description | What It Contains |
|------|-------------|------------------|
| `genotypes.csv` | DNA data | 300 people × 200 SNPs, values 0/1/2 |

**What the Data Means:**
- **0**: No copies of variant allele (like AA)
- **1**: One copy (like Aa)  
- **2**: Two copies (like aa)

### Process

1. **Data Generation**: We create 4 distinct populations with different allele frequencies, plus 25 admixed individuals.

2. **Autoencoder Training**: 
   - **Input Layer**: 200 neurons (one per SNP)
   - **Hidden Layer**: 64 neurons
   - **Latent Layer**: 10 neurons (the compressed representation)
   - **Output Layer**: 200 neurons (reconstructing the input)
   - The network learns to compress and reconstruct

3. **Latent Space Extraction**: After training, we extract the 10 numbers from the middle layer for each person.

4. **K-means Clustering**: We run K-means with K=4 on the latent space.

5. **t-SNE Visualization**: We project the 10-dimensional latent space to 2D for visualization.

### Output
| File | Description | What It Shows |
|------|-------------|---------------|
| `deep_clustering.png` | Four-panel visualization | True labels, predicted clusters, latent space heatmap, confusion matrix |
| `latent_representation.csv` | Latent coordinates | The 10-number representation for each person |
| `cluster_assignments.csv` | Cluster assignments | Which cluster each person was assigned to |

**Visualization Explanation:**

**Panel 1: True Population Labels**
- Shows where people from each true population fall in the latent space
- Different colors = different populations
- If clustering works well, each color should form a tight group

**Panel 2: Predicted Clusters**
- Shows the clusters found by K-means
- The ARI score (0-1) tells us how well they match the true labels
- ARI > 0.8 is considered very good

**Panel 3: Latent Space Heatmap**
- Shows the 10 latent dimensions for the first 50 people
- Each row is a latent feature
- Each column is a person
- Blue/red show feature values (standardized)

**Panel 4: Confusion Matrix**
- Compares true populations (rows) to predicted clusters (columns)
- Diagonal = correct assignments
- Off-diagonal = misclassifications
- Perfect clustering would have all values on the diagonal

## Running the Example

```bash
cd examples/deep-clustering
python run_deep.py
```

## Expected Runtime
- Data generation: < 1 second
- Autoencoder training: ~30-60 seconds (iterative learning)
- Clustering: ~1 second
- t-SNE: ~5 seconds
- Visualization: ~2 seconds

## Acceptance Criteria
- [x] **ARI > 0.8**: The clustering matches the true populations well.
- [x] **Latent space extracted**: 10-dimensional representation created.
- [x] **t-SNE visualization**: Shows clear separation of populations.
- [x] **Admixed samples identified**: Mixed individuals are correctly placed between populations.

## Tools Used
- **scikit-learn MLPRegressor**: For the autoencoder neural network
- **scikit-learn KMeans**: For clustering in latent space
- **scikit-learn TSNE**: For 2D visualization
- **StandardScaler**: To normalize the data before training
- **numpy/pandas**: For data handling
- **matplotlib**: For visualization

## How Autoencoders Work

Think of an autoencoder like a data compression tool:

1. **Encoder**: Compresses the data (200 SNPs → 10 numbers)
2. **Latent Space**: The compressed representation
3. **Decoder**: Reconstructs the data (10 numbers → 200 SNPs)
4. **Training**: The network learns to minimize reconstruction error

The key insight: the latent space must capture the most important information to reconstruct the input accurately. This makes it great for clustering!

## Real-World Application
In a real study, you might use this to:
- Analyze 500,000 SNPs and reduce to 50 latent dimensions
- Discover cryptic population structure
- Identify individuals with mixed ancestry
- Pre-process data for downstream GWAS (correct for population stratification)

## Troubleshooting
**Training takes too long**: Reduce max_iter in MLPRegressor or use fewer SNPs.

**ARI is low (<0.5)**: The populations might be too similar or the autoencoder needs more training.

**t-SNE looks messy**: Try different perplexity values (10-50) or check that populations are actually distinct.

**All samples in one cluster**: Check that you have enough training iterations.
