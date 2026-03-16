import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import powerlaw

# =========================
# 路径
# =========================
DATA_FILE = r'G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\0_data\primary_simple_gene_cancer_edges.xlsx'
OUT_DIR = r'G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M1'
os.makedirs(OUT_DIR, exist_ok=True)

FIG_BAR = os.path.join(OUT_DIR, 'figure_M1_degree_frequency_barplot.png')
FIG_CCDF2 = os.path.join(OUT_DIR, 'figure_M1_manual_ccdf_loglog.png')
FREQ_TABLE_FILE = os.path.join(OUT_DIR, 'degree_frequency_table.xlsx')

# =========================
# 读取数据
# =========================
df = pd.read_excel(DATA_FILE)

degree_df = (
    df.groupby('Gene', as_index=False)['Cancer']
    .nunique()
    .rename(columns={'Cancer': 'Degree'})
)

degrees = np.sort(degree_df['Degree'].values)
n = len(degrees)

# =========================
# 1. degree frequency bar plot
# =========================
freq_df = (
    degree_df.groupby('Degree', as_index=False)
    .size()
    .rename(columns={'size': 'Gene_count'})
    .sort_values('Degree')
)

freq_df.to_excel(FREQ_TABLE_FILE, index=False)

plt.figure(figsize=(7, 5))
plt.bar(freq_df['Degree'], freq_df['Gene_count'], width=0.7, edgecolor='black')
plt.xlabel('Gene degree (number of distinct cancer types)')
plt.ylabel('Number of genes')
plt.title('Degree frequency in the primary simple gene–cancer network')
plt.xticks(freq_df['Degree'])
plt.tight_layout()
plt.savefig(FIG_BAR, dpi=300)
plt.close()

# =========================
# 2. 手工经验 CCDF
#    P(K >= k)
# =========================
unique_k = np.sort(np.unique(degrees))
ccdf = np.array([(degrees >= k).sum() / n for k in unique_k])

# power-law fit
fit = powerlaw.Fit(degrees, discrete=True, verbose=False)
alpha = fit.power_law.alpha
xmin = fit.power_law.xmin

# 只画 tail 的拟合线（k >= xmin）
k_tail = unique_k[unique_k >= xmin]
if len(k_tail) > 0:
    # 连续比例常数只用于可视化，不做严格密度解释
    y_fit = (k_tail / xmin) ** (-(alpha - 1))
    # 让拟合线在 xmin 点与经验 CCDF 对齐
    ccdf_xmin = ccdf[np.where(unique_k == xmin)[0][0]]
    y_fit = y_fit * ccdf_xmin
else:
    k_tail = np.array([])
    y_fit = np.array([])

plt.figure(figsize=(7, 5))
plt.loglog(unique_k, ccdf, 'o-', color='black', label='Empirical CCDF')
if len(k_tail) > 0:
    plt.loglog(k_tail, y_fit, '--', color='red', linewidth=2, label='Power-law tail fit')
plt.xlabel('Degree')
plt.ylabel('P(K ≥ k)')
plt.title('Empirical CCDF of gene degree distribution')
plt.legend()
plt.tight_layout()
plt.savefig(FIG_CCDF2, dpi=300)
plt.close()

print("Saved:")
print(FIG_BAR)
print(FIG_CCDF2)
print("alpha =", alpha)
print("xmin =", xmin)
