# Before pathway clustering, plotting the pathways*cancers matrix and cancer similarity matrix,
# creating heatmaps

from plot_boot import *
boot_plot_folder()      # coloring scheme using cycler
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from definitions import *

# --- Setup ---
output_folder = PLOTS_P
os.makedirs(output_folder, exist_ok=True)
Q_SIG = 0.05   # Significance threshold

# --- 1. Load Data ---
csv_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
delta_frames = []
qval_frames  = []

for csv_path in csv_files:
    cancer_type = os.path.splitext(os.path.basename(csv_path))[0]
    df = pd.read_csv(csv_path)
    # Ensure no duplicate pathways per file
    df = df.drop_duplicates(subset='pathway', keep='first').set_index('pathway')
    delta_frames.append(df['delta_means'].rename(cancer_type))
    qval_frames.append(df['q_value'].rename(cancer_type))

# Construct master matrices (Pathways as columns, Cancers as rows)
delta_mat = pd.concat(delta_frames, axis=1).T  
qval_mat  = pd.concat(qval_frames,  axis=1).T

# --- 2. Robust Bi-directional Filtering ---
# Step A: Keep pathways significant in AT LEAST one cancer
sig_pathway_mask = (qval_mat <= Q_SIG).any(axis=0)
delta_filt = delta_mat.loc[:, sig_pathway_mask]
qval_filt  = qval_mat.loc[:, sig_pathway_mask]

# Step B: Keep cancers significant in AT LEAST one pathway (removes "empty" cancers)
sig_cancer_mask = (qval_filt <= Q_SIG).any(axis=1)
delta_filt = delta_filt.loc[sig_cancer_mask, :]
qval_filt  = qval_filt.loc[sig_cancer_mask, :]

# Step C: Final Cleanup (Remove any remaining all-NaN rows/cols that might break clustering)
delta_filt = delta_filt.dropna(how='all', axis=0).dropna(how='all', axis=1)
qval_filt  = qval_filt.loc[delta_filt.index, delta_filt.columns]

print(f"Filtered Matrix: {delta_filt.shape[0]} Cancers x {delta_filt.shape[1]} Pathways")

# --- 3. Clustering Helper ---
def hierarchical_order(df, axis=0, method='ward', metric='euclidean'):
    """Return reordered index/columns using Ward hierarchical clustering."""
    # Fill NaNs with 0 for the purpose of distance calculation
    data = df.fillna(0).values if axis == 0 else df.fillna(0).T.values
    if data.shape[0] < 2:
        return df.index if axis == 0 else df.columns
    Z = linkage(data, method=method, metric=metric)
    order = leaves_list(Z)
    return (df.index if axis == 0 else df.columns)[order]

# --- 4. Plotting Stage A: Clustered Heatmap ---
print("Plotting A: Clustered Heatmap...")
row_order = hierarchical_order(delta_filt, axis=0)
col_order = hierarchical_order(delta_filt, axis=1)

delta_ordered = delta_filt.loc[row_order, col_order]
qval_ordered  = qval_filt.loc[row_order, col_order]
delta_display = delta_ordered.fillna(0)

# Set color scale based on 95th percentile of absolute values
vmax = np.nanpercentile(np.abs(delta_filt.values), 95)
vmin = -vmax

fig, ax = plt.subplots(figsize=(min(2 + delta_display.shape[1] * 0.15, 45),
                                 min(2 + delta_display.shape[0] * 0.30, 25)))

# Use interpolation='nearest' and edgecolor='none' for clean high-density rendering
im = ax.imshow(delta_display.values, aspect='auto', cmap='RdBu_r',
               vmin=vmin, vmax=vmax, interpolation='nearest')

# Overlay significance dots
sig_mask = qval_ordered.values <= Q_SIG
rows_sig, cols_sig = np.where(sig_mask)
ax.scatter(cols_sig, rows_sig, s=7, c='black', marker='o', zorder=3, alpha=0.8)

# Formatting
ax.set_xticks(range(delta_display.shape[1]))
ax.set_xticklabels(delta_display.columns, rotation=90, fontsize=5)
ax.set_yticks(range(delta_display.shape[0]))
ax.set_yticklabels(delta_display.index, fontsize=9)

plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label='Delta Means')
ax.legend(handles=[mpatches.Patch(color='black', label=f'Q <= {Q_SIG}')], 
          loc='upper right', bbox_to_anchor=(1.05, 1))

ax.set_title(f'Clustered Heatmap: Delta Means\n({delta_display.shape[0]} Cancers, {delta_display.shape[1]} Sig Pathways)', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'p11A_clustered_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()

# --- 5. Plotting Stage B: Similarity Matrix ---
print("Plotting B: Cancer-Cancer Similarity Matrix...")
# Fill missing values with column mean for similarity calculation
delta_filled = delta_filt.fillna(delta_filt.mean(axis=0))
# Cosine similarity
normed = delta_filled.values / (np.linalg.norm(delta_filled.values, axis=1, keepdims=True) + 1e-12)
sim_mat = pd.DataFrame(normed @ normed.T, index=delta_filled.index, columns=delta_filled.index)

cancer_order = hierarchical_order(sim_mat, axis=0)
sim_ordered  = sim_mat.loc[cancer_order, cancer_order]

fig, ax = plt.subplots(figsize=(max(10, len(cancer_order)*0.4), max(9, len(cancer_order)*0.4)))
im = ax.imshow(sim_ordered.values, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1)

ax.set_xticks(range(len(cancer_order)))
ax.set_xticklabels(sim_ordered.columns, rotation=90, fontsize=8)
ax.set_yticks(range(len(cancer_order)))
ax.set_yticklabels(sim_ordered.index, fontsize=8)

plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Cosine Similarity')
ax.set_title('Cancer-Cancer Pathway Profile Similarity', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'p11B_cancer_similarity_matrix.png'), dpi=150)
plt.close()

# --- 6. Plotting Stage D: PCA ---
def plot_pca(df_input, pca_name, scree_name, suffix=''):
    print(f"Plotting D: PCA Scatter {suffix}...")
    # Pre-processing: Fill NaNs with column means, then scale
    X = df_input.fillna(df_input.mean(axis=0))
    X_scaled = StandardScaler().fit_transform(X.values)

    pca = PCA(n_components=min(10, X_scaled.shape[0] - 1))
    coords = pca.fit_transform(X_scaled)
    
    var_exp = pca.explained_variance_ratio_ * 100

    # PCA Scatter
    plt.figure(figsize=(12, 9))
    sns.set_style("whitegrid")
    
    # Use a diverse color palette
    palette = sns.color_cycle if hasattr(sns, 'color_cycle') else plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    for i, cancer in enumerate(df_input.index):
        plt.scatter(coords[i, 0], coords[i, 1], s=100, edgecolors='white', alpha=0.8, color=palette[i % len(palette)])
        plt.annotate(cancer, (coords[i, 0], coords[i, 1]), fontsize=8, xytext=(0,5), textcoords='offset points', ha='center')

    plt.axhline(0, color='grey', lw=1, ls='--')
    plt.axvline(0, color='grey', lw=1, ls='--')
    plt.xlabel(f'PC1 ({var_exp[0]:.1f}%)', fontsize=12)
    plt.ylabel(f'PC2 ({var_exp[1]:.1f}%)', fontsize=12)
    plt.title(f'PCA: Cancer Distribution in Pathway Space{suffix}', fontsize=14)
    sns.despine()
    plt.savefig(os.path.join(output_folder, pca_name), dpi=150, bbox_inches='tight')
    plt.close()

    # Scree Plot
    plt.figure(figsize=(6, 4))
    plt.bar(range(1, len(var_exp)+1), var_exp, color='steelblue')
    plt.ylabel('Explained Variance (%)')
    plt.xlabel('Principal Component')
    plt.title(f'Scree Plot{suffix}')
    sns.despine()
    plt.savefig(os.path.join(output_folder, scree_name), dpi=150)
    plt.close()

# Run PCA on full filtered set
plot_pca(delta_filt, 'p11D_pca_cancers.png', 'p11D_pca_scree.png')

# Run PCA excluding UM (outlier check)
if 'um' in delta_filt.index:
    plot_pca(delta_filt.drop('um'), 'p11D_pca_cancers_no_UM.png', 'p11D_pca_scree_no_UM.png', suffix=' (excl. UM)')

print(f"\nDone. All plots saved to: {output_folder}")