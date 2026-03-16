import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import powerlaw

# =========================
# 0. 路径设置
# =========================
DATA_FILE = r'G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\0_data\primary_simple_gene_cancer_edges.xlsx'
OUT_DIR = r'G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M1'

os.makedirs(OUT_DIR, exist_ok=True)

DEGREE_TABLE_FILE = os.path.join(OUT_DIR, 'gene_degree_table.xlsx')
FREQ_TABLE_FILE = os.path.join(OUT_DIR, 'degree_frequency_table.xlsx')
SUMMARY_FILE = os.path.join(OUT_DIR, 'powerlaw_fit_summary.txt')
COMPARE_FILE = os.path.join(OUT_DIR, 'powerlaw_distribution_comparison.xlsx')
FIG_HIST = os.path.join(OUT_DIR, 'figure_M1_degree_histogram.png')
FIG_CCDF = os.path.join(OUT_DIR, 'figure_M1_ccdf_loglog.png')

# =========================
# 1. 读取数据
# =========================
df = pd.read_excel(DATA_FILE)

print("=" * 70)
print("Input data preview")
print("=" * 70)
print(df.head())
print("\nShape:", df.shape)
print("Columns:", list(df.columns))

required_cols = {'Gene', 'Cancer'}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    raise ValueError(f"输入文件缺少必要列: {missing_cols}")

# =========================
# 2. 计算 gene degree
# =========================
# 虽然 D2 已经是 unique (gene, cancer)，这里仍用 nunique() 更稳妥
degree_df = (
    df.groupby('Gene', as_index=False)['Cancer']
    .nunique()
    .rename(columns={'Cancer': 'Degree'})
    .sort_values(by=['Degree', 'Gene'], ascending=[False, True])
    .reset_index(drop=True)
)

degree_df.to_excel(DEGREE_TABLE_FILE, index=False)

degrees = degree_df['Degree'].values

# degree frequency table
freq_df = (
    degree_df.groupby('Degree', as_index=False)
    .size()
    .rename(columns={'size': 'Gene_count'})
    .sort_values(by='Degree')
    .reset_index(drop=True)
)
freq_df.to_excel(FREQ_TABLE_FILE, index=False)

# =========================
# 3. 基础统计
# =========================
n_genes = len(degree_df)
n_edges = len(df)
n_cancers = df['Cancer'].nunique()

degree_min = int(np.min(degrees))
degree_max = int(np.max(degrees))
degree_mean = float(np.mean(degrees))
degree_median = float(np.median(degrees))

print("\n" + "=" * 70)
print("Degree summary")
print("=" * 70)
print(f"Genes: {n_genes}")
print(f"Cancers: {n_cancers}")
print(f"Edges: {n_edges}")
print(f"Degree min: {degree_min}")
print(f"Degree median: {degree_median}")
print(f"Degree mean: {degree_mean:.4f}")
print(f"Degree max: {degree_max}")

# =========================
# 4. power-law 拟合
# =========================
# discrete=True 因为 degree 是离散值
fit = powerlaw.Fit(degrees, discrete=True, verbose=False)

alpha = fit.power_law.alpha
xmin = fit.power_law.xmin
D_ks = fit.power_law.KS()

# powerlaw 包里 sigma 常作为 alpha 的标准误
sigma = fit.power_law.sigma

# 近似 95% CI（正态近似）
alpha_ci_low = alpha - 1.96 * sigma
alpha_ci_high = alpha + 1.96 * sigma

# goodness-of-fit p-value（bootstrap）
# 注意：small samples 时这个值解释也要谨慎
try:
    R_gof, p_gof = fit.distribution_compare('power_law', 'lognormal', normalized_ratio=True)
except Exception:
    R_gof, p_gof = np.nan, np.nan

# =========================
# 5. 分布比较：Clauset-style likelihood ratio tests
# =========================
comparisons = []

def safe_compare(dist1, dist2):
    try:
        R, p = fit.distribution_compare(dist1, dist2, normalized_ratio=True)
        if np.isnan(R):
            favored = 'undetermined'
        elif R > 0:
            favored = dist1
        elif R < 0:
            favored = dist2
        else:
            favored = 'tie'
        return R, p, favored
    except Exception as e:
        return np.nan, np.nan, f'error: {str(e)}'

for d1, d2 in [
    ('power_law', 'lognormal'),
    ('power_law', 'exponential'),
    ('power_law', 'truncated_power_law'),
]:
    R, p, favored = safe_compare(d1, d2)
    comparisons.append({
        'distribution_1': d1,
        'distribution_2': d2,
        'loglikelihood_ratio_R': R,
        'p_value': p,
        'favored_model': favored
    })

compare_df = pd.DataFrame(comparisons)
compare_df.to_excel(COMPARE_FILE, index=False)

print("\n" + "=" * 70)
print("Distribution comparison")
print("=" * 70)
print(compare_df)

# =========================
# 6. 作图 1：degree histogram
# =========================
plt.figure(figsize=(7, 5))
bins = np.arange(degree_min, degree_max + 2) - 0.5
plt.hist(degrees, bins=bins, edgecolor='black', alpha=0.85)
plt.xlabel('Gene degree (number of distinct cancer types)')
plt.ylabel('Number of genes')
plt.title('Degree distribution of the primary simple gene–cancer network')
plt.xticks(range(degree_min, degree_max + 1))
plt.tight_layout()
plt.savefig(FIG_HIST, dpi=300)
plt.close()

# =========================
# 7. 作图 2：CCDF log-log
# =========================
plt.figure(figsize=(7, 5))
fit.plot_ccdf(color='black', linewidth=2, label='Empirical data')
fit.power_law.plot_ccdf(color='red', linestyle='--', linewidth=2, label='Power-law fit')
plt.xlabel('Degree')
plt.ylabel('CCDF')
plt.title('Log-log CCDF of gene degree distribution')
plt.legend()
plt.tight_layout()
plt.savefig(FIG_CCDF, dpi=300)
plt.close()

# =========================
# 8. 输出 summary txt
# =========================
hub_k2 = int((degree_df['Degree'] >= 2).sum())
hub_k3 = int((degree_df['Degree'] >= 3).sum())
hub_k4 = int((degree_df['Degree'] >= 4).sum())
hub_k5 = int((degree_df['Degree'] >= 5).sum())

top_genes = degree_df.head(20)

summary_lines = []
summary_lines.append("M1: Primary network degree distribution and power-law analysis")
summary_lines.append("=" * 70)
summary_lines.append(f"Input file: {DATA_FILE}")
summary_lines.append(f"Number of genes: {n_genes}")
summary_lines.append(f"Number of cancers: {n_cancers}")
summary_lines.append(f"Number of edges: {n_edges}")
summary_lines.append("")
summary_lines.append("Degree summary")
summary_lines.append("-" * 70)
summary_lines.append(f"Minimum degree: {degree_min}")
summary_lines.append(f"Median degree: {degree_median}")
summary_lines.append(f"Mean degree: {degree_mean:.4f}")
summary_lines.append(f"Maximum degree: {degree_max}")
summary_lines.append("")
summary_lines.append("Hub counts under different thresholds")
summary_lines.append("-" * 70)
summary_lines.append(f"Degree >= 2: {hub_k2}")
summary_lines.append(f"Degree >= 3: {hub_k3}")
summary_lines.append(f"Degree >= 4: {hub_k4}")
summary_lines.append(f"Degree >= 5: {hub_k5}")
summary_lines.append("")
summary_lines.append("Power-law fit")
summary_lines.append("-" * 70)
summary_lines.append(f"alpha (gamma): {alpha:.4f}")
summary_lines.append(f"sigma (SE): {sigma:.4f}")
summary_lines.append(f"Approx. 95% CI for alpha: [{alpha_ci_low:.4f}, {alpha_ci_high:.4f}]")
summary_lines.append(f"xmin: {xmin}")
summary_lines.append(f"KS statistic: {D_ks:.4f}")
summary_lines.append("")
summary_lines.append("Likelihood ratio comparisons")
summary_lines.append("-" * 70)
for _, row in compare_df.iterrows():
    summary_lines.append(
        f"{row['distribution_1']} vs {row['distribution_2']}: "
        f"R={row['loglikelihood_ratio_R']}, p={row['p_value']}, favored={row['favored_model']}"
    )
summary_lines.append("")
summary_lines.append("Top genes by degree")
summary_lines.append("-" * 70)
for _, row in top_genes.iterrows():
    summary_lines.append(f"{row['Gene']}: {row['Degree']}")

with open(SUMMARY_FILE, 'w', encoding='utf-8') as f:
    f.write('\n'.join(summary_lines))

# =========================
# 9. 控制台输出
# =========================
print("\n" + "=" * 70)
print("Saved files")
print("=" * 70)
print("Degree table:", DEGREE_TABLE_FILE)
print("Degree frequency table:", FREQ_TABLE_FILE)
print("Summary:", SUMMARY_FILE)
print("Comparison table:", COMPARE_FILE)
print("Histogram figure:", FIG_HIST)
print("CCDF figure:", FIG_CCDF)

print("\n✅ M1 analysis completed.")
