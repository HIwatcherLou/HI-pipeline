import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky 
import warnings

warnings.filterwarnings('ignore')

# 1. 配置路径和参数
sofia_cat_path = '/home/jclou/loujc_work/sofia/GAMA/1300-1350/Dec-0011_09_05_arcdrift/Dec-0011_09_05_arcdrift_cat.txt'
gama_local_path = '/home/jclou/loujc_work/GAMA_data/GAMA_DR4_SpecObj.csv' 

top_n = 217 
search_radius = 3 * u.arcmin 
f_rest = 1420.40575 

# 2. 加载数据
print("正在加载 SoFiA 和 GAMA 数据...")
columns = ["name", "id", "x", "y", "z", "x_min", "x_max", "y_min", "y_max", "z_min", "z_max", 
           "n_pix", "f_min", "f_max", "f_sum", "rel", "flag", "rms", "w20", "w50", "wm50", 
           "z_w20", "z_w50", "z_wm50", "ell_maj", "ell_min", "ell_pa", "ell3s_maj", "ell3s_min", 
           "ell3s_pa", "kin_pa", "err_x", "err_y", "err_z", "err_f_sum", "snr", "snr_max", 
           "ra", "dec", "freq", "x_peak", "y_peak", "z_peak", "ra_peak", "dec_peak", "freq_peak"]
df_sofia = pd.read_csv(sofia_cat_path, sep='\s+', comment='#', names=columns)
df_top = df_sofia.sort_values('snr', ascending=False).head(top_n).reset_index(drop=True)

# 统一计算 SoFiA 源的 z_HI
df_top['f_obs_MHz'] = df_top['freq'] / 1e6
df_top['z_hi'] = (f_rest / df_top['f_obs_MHz']) - 1

# 加载本地 GAMA 数据库
df_gama = pd.read_csv(gama_local_path)

# 动态识别 GAMA 中的核心列
id_col = 'CATAID' if 'CATAID' in df_gama.columns else df_gama.columns[df_gama.columns.str.lower() == 'cataid'][0]
ra_col = 'RAJ2000' if 'RAJ2000' in df_gama.columns else df_gama.columns[df_gama.columns.str.lower().str.contains('ra')][0]
dec_col = 'DEJ2000' if 'DEJ2000' in df_gama.columns else df_gama.columns[df_gama.columns.str.lower().str.contains('de')][0]
z_col = 'z' if 'z' in df_gama.columns else df_gama.columns[df_gama.columns.str.lower() == 'z'][0]

# 3. 构建天球坐标 (SkyCoord)
coords_sofia = SkyCoord(ra=df_top['ra'].values*u.deg, dec=df_top['dec'].values*u.deg)
coords_gama = SkyCoord(ra=df_gama[ra_col].values*u.deg, dec=df_gama[dec_col].values*u.deg)

# 4. 寻找半径范围内的【所有】候选源
print("正在执行空间索引矩阵交叉匹配...")
idx_sofia, idx_gama, d2d, d3d = search_around_sky(coords_sofia, coords_gama, search_radius)

# 5. 红移筛选与最佳距离优选
best_matches = {}

for i in range(len(idx_sofia)):
    s_idx = idx_sofia[i]  
    g_idx = idx_gama[i]  
    dist = d2d[i]         
    
    z_hi = df_top.iloc[s_idx]['z_hi']
    gama_z = df_gama.iloc[g_idx][z_col]
    gama_id = int(df_gama.iloc[g_idx][id_col])
    
    # 必须满足红移误差要求 delta_z < 0.002
    if abs(gama_z - z_hi) < 0.002:
        # 寻找满足红移条件下，空间距离最近的宿主星系
        if s_idx not in best_matches or dist < best_matches[s_idx]['dist']:
            best_matches[s_idx] = {
                'gama_id': gama_id,
                'z_gama': gama_z,
                'dist': dist
            }

# 6. 生成最终的统计表格
results = []
for i in range(len(df_top)):
    s_id = int(df_top.iloc[i]['id'])
    z_hi = df_top.iloc[i]['z_hi']
    snr = df_top.iloc[i]['snr']
    
    if i in best_matches:
        match = best_matches[i]
        results.append({
            'ID': s_id, 'SNR': snr, 'z_HI': z_hi, 
            'GAMA_Match': f"GAMA_{match['gama_id']}", 
            'z_GAMA': match['z_gama'], 
            'Dist(arcmin)': round(match['dist'].arcmin, 3), 
            'Status': '✅ 成功配对'
        })
    else:
        results.append({
            'ID': s_id, 'SNR': snr, 'z_HI': z_hi, 
            'GAMA_Match': 'N/A', 'z_GAMA': None, 'Dist(arcmin)': None, 
            'Status': '❌ 范围内无对应红移源'
        })

# 7. 展示结果
print("\n" + "="*50)
print("检索完成！结果汇总如下：")
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
summary_df = pd.DataFrame(results)
display(summary_df)

# 统计成功率
success_count = len([r for r in results if '✅' in r['Status']])
print(f"\n 总结: 在 {top_n} 个源中，成功找到了 {success_count} 个具有光学对应体及红移一致的星系！")
