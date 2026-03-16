import pandas as pd
import os

# =========================
# 1. 路径设置
# =========================
BASE_DIR = r'G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\0_data'

RAW_FILE = os.path.join(BASE_DIR, 'raw_association_revision.xlsx')
MAP_FILE = os.path.join(BASE_DIR, 'cancer_harmonization_map.xlsx')

OUT_STUDY = os.path.join(BASE_DIR, 'study_level_unique_gene_cancer_pmid.xlsx')
OUT_PRIMARY = os.path.join(BASE_DIR, 'primary_simple_gene_cancer_edges.xlsx')
OUT_HARMONIZED_ALL = os.path.join(BASE_DIR, 'harmonized_association_revision.xlsx')

# =========================
# 2. 读取文件
# =========================
raw_df = pd.read_excel(RAW_FILE)
map_df = pd.read_excel(MAP_FILE)

print("=" * 60)
print("Raw data preview")
print("=" * 60)
print(raw_df.head())
print("\nRaw data shape:", raw_df.shape)

print("\n" + "=" * 60)
print("Cancer harmonization map preview")
print("=" * 60)
print(map_df.head())
print("\nMap shape:", map_df.shape)

# =========================
# 3. 基础检查
# =========================
required_raw_cols = {'Source', 'Target', 'PMID'}
required_map_cols = {'original_cancer_name', 'harmonized_cancer_name', 'include_in_primary_analysis'}

missing_raw = required_raw_cols - set(raw_df.columns)
missing_map = required_map_cols - set(map_df.columns)

if missing_raw:
    raise ValueError(f"raw_association_revision.xlsx 缺少列: {missing_raw}")

if missing_map:
    raise ValueError(f"cancer_harmonization_map.xlsx 缺少列: {missing_map}")

# 去掉列名首尾空格
raw_df.columns = raw_df.columns.str.strip()
map_df.columns = map_df.columns.str.strip()

# 统一字符串首尾空格
raw_df['Target'] = raw_df['Target'].astype(str).str.strip()
raw_df['Source'] = raw_df['Source'].astype(str).str.strip()
raw_df['PMID'] = raw_df['PMID'].astype(str).str.strip()

map_df['original_cancer_name'] = map_df['original_cancer_name'].astype(str).str.strip()
map_df['harmonized_cancer_name'] = map_df['harmonized_cancer_name'].astype(str).str.strip()
map_df['include_in_primary_analysis'] = map_df['include_in_primary_analysis'].astype(str).str.strip()

# =========================
# 4. 检查 raw 里的癌种标签是否都在 map 中
# =========================
raw_cancers = set(raw_df['Target'].unique())
mapped_cancers = set(map_df['original_cancer_name'].unique())

unmapped = raw_cancers - mapped_cancers

print("\n" + "=" * 60)
print("Cancer label coverage check")
print("=" * 60)
print(f"Raw unique cancer labels: {len(raw_cancers)}")
print(f"Mapped original cancer labels: {len(mapped_cancers)}")

if unmapped:
    print("\n以下癌种标签在 raw 文件中出现，但不在 harmonization map 中：")
    for x in sorted(unmapped):
        print(" -", x)
    raise ValueError("请先补全 cancer_harmonization_map.xlsx 后再运行。")
else:
    print("\n所有 raw 中的癌种标签都已在 harmonization map 中找到。")

# =========================
# 5. 合并 harmonization map
# =========================
harmonized_df = raw_df.merge(
    map_df,
    how='left',
    left_on='Target',
    right_on='original_cancer_name'
)

# 检查 merge 后是否有缺失
if harmonized_df['harmonized_cancer_name'].isna().any():
    missing_rows = harmonized_df[harmonized_df['harmonized_cancer_name'].isna()]
    print("\n以下记录未成功匹配 harmonized cancer name：")
    print(missing_rows[['Source', 'Target', 'PMID']].drop_duplicates().head(20))
    raise ValueError("存在未匹配的 cancer label，请检查 harmonization map。")

# =========================
# 6. 输出 harmonized 全量文件（可选但非常推荐）
# =========================
harmonized_df.to_excel(OUT_HARMONIZED_ALL, index=False)
print("\n已输出 harmonized 全量文件：")
print(OUT_HARMONIZED_ALL)

# =========================
# 7. 生成 study-level 文件
#    规则：unique (gene, harmonized_cancer, PMID)
# =========================
study_df = harmonized_df.copy()

study_df = study_df.drop_duplicates(subset=['Source', 'harmonized_cancer_name', 'PMID'])

# 可选：保留核心列在前面
study_cols_front = ['Source', 'Target', 'harmonized_cancer_name', 'PMID', 'include_in_primary_analysis']
other_cols = [c for c in study_df.columns if c not in study_cols_front]
study_df = study_df[study_cols_front + other_cols]

study_df.to_excel(OUT_STUDY, index=False)

print("\n" + "=" * 60)
print("Study-level dataset generated")
print("=" * 60)
print("Shape:", study_df.shape)
print("Unique genes:", study_df['Source'].nunique())
print("Unique harmonized cancer labels:", study_df['harmonized_cancer_name'].nunique())
print("Saved to:", OUT_STUDY)

# =========================
# 8. 生成 primary simple network 文件
#    规则：
#    - include_in_primary_analysis == Yes
#    - unique (gene, harmonized_cancer)
# =========================
primary_df = harmonized_df.copy()

primary_df = primary_df[
    primary_df['include_in_primary_analysis'].str.lower() == 'yes'
].copy()

primary_df = primary_df.drop_duplicates(subset=['Source', 'harmonized_cancer_name'])

# 只保留主网络需要的核心列
primary_df = primary_df[['Source', 'harmonized_cancer_name']].copy()
primary_df.columns = ['Gene', 'Cancer']

# 排序（可选）
primary_df = primary_df.sort_values(by=['Gene', 'Cancer']).reset_index(drop=True)

primary_df.to_excel(OUT_PRIMARY, index=False)

print("\n" + "=" * 60)
print("Primary simple network dataset generated")
print("=" * 60)
print("Shape:", primary_df.shape)
print("Unique genes:", primary_df['Gene'].nunique())
print("Unique cancers:", primary_df['Cancer'].nunique())
print("Saved to:", OUT_PRIMARY)

# =========================
# 9. 简单统计输出
# =========================
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Raw rows: {len(raw_df)}")
print(f"Harmonized rows: {len(harmonized_df)}")
print(f"Study-level unique (gene, cancer, PMID): {len(study_df)}")
print(f"Primary simple unique (gene, cancer): {len(primary_df)}")

print("\nPrimary cancer labels:")
for c in sorted(primary_df['Cancer'].unique()):
    print(" -", c)

print("\n✅ 全部完成。")
