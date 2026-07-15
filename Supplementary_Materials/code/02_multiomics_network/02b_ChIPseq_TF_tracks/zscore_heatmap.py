#!/usr/bin/env python3
"""
Z-score normalize each sample in a deepTools computeMatrix output,
write a new .gz matrix file, then run plotHeatmap.

Usage:
    python zscore_heatmap.py <input_matrix.gz> <output_prefix>
"""
import gzip, json, sys, os, subprocess
import numpy as np

matrix_file = sys.argv[1] if len(sys.argv) > 1 else "analysis/TF_binding_matrix_v2.gz"
out_prefix = sys.argv[2] if len(sys.argv) > 2 else "analysis/TF_binding_matrix_v3_zscore"

# Read the matrix (JSON format)
with gzip.open(matrix_file, 'rt') as f:
    data = json.load(f)

# Extract the matrix data (list of lists: [n_samples][n_regions][n_bins])
# deepTools stores it as a list: for each sample, a list of region arrays
scores = data['matrix']  # list of (n_regions, n_bins) arrays per sample
sample_labels = data['sample_labels']
n_samples = len(scores)
n_regions = len(scores[0])

print(f"Matrix: {n_samples} samples, {n_regions} regions, {len(scores[0][0])} bins each")

# Z-score normalize per sample: for each sample, compute mean & std across all bins
# of all regions, then z-score normalize
for s in range(n_samples):
    # Convert to numpy array: shape (n_regions, n_bins)
    mat = np.array(scores[s], dtype=np.float64)

    # Compute global mean and std for this sample (across all bins and regions)
    global_mean = np.mean(mat)
    global_std = np.std(mat)

    # Z-score normalize
    mat_z = (mat - global_mean) / global_std

    # Store back
    scores[s] = mat_z.tolist()

    print(f"  {sample_labels[s]}: mean={global_mean:.4f}, std={global_std:.4f}")
    print(f"    After z-score: min={np.min(mat_z):.2f}, max={np.max(mat_z):.2f}")

# Update the matrix in the data structure
data['matrix'] = scores

# Write modified matrix to new file
out_file = out_prefix + '.gz'
with gzip.open(out_file, 'wt') as f:
    json.dump(data, f)

print(f"\nWritten: {out_file} ({os.path.getsize(out_file)} bytes)")

# Now run plotHeatmap on the z-score normalized matrix
heatmap_file = out_prefix + '_heatmap.pdf'
cmd = [
    'plotHeatmap',
    '-m', out_file,
    '-out', heatmap_file,
    '--colorMap', 'RdYlBu_r',
    '--heatmapHeight', '12',
    '--heatmapWidth', '3',
    '--yAxisLabel', 'Signal (z-score)',
    '--regionsLabel', 'Prognostic eRNAs',
    '--whatToShow', 'heatmap and colorbar',
    '--dpi', '300'
]
for label in sample_labels:
    cmd.extend(['--samplesLabel', label])

print(f"Running: {' '.join(cmd[:6])} ...")
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode == 0:
    print(f"✅ Heatmap: {heatmap_file}")
else:
    print(f"❌ Error: {result.stderr}")

print("\nDone!")
