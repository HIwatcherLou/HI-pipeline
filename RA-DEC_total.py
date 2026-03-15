#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 尝试导入 hifast
try:
    from hifast.core.radec import get_radec
except ImportError:
    print("错误: 无法导入 hifast。请确保你处于 hifast 环境中。")
    sys.exit(1)

# --- 定义一组深色、高对比度的颜色库 (18种) ---
# 这些颜色经过挑选，没有白色，没有浅色，在白底上对比度极高
DARK_COLORS = [
    '#000000', # 黑色
    '#FF0000', # 纯红
    '#0000FF', # 纯蓝
    '#008000', # 纯绿
    '#800080', # 紫色
    '#A52A2A', # 棕色
    '#FF8C00', # 深橙色
    '#008080', # 蓝绿色 (Teal)
    '#000080', # 海军蓝 (Navy)
    '#800000', # 栗色 (Maroon)
    '#556B2F', # 深橄榄绿
    '#4B0082', # 靛青 (Indigo)
    '#DC143C', # 猩红 (Crimson)
    '#8B4513', # 鞍褐色
    '#2F4F4F', # 深岩灰
    '#C71585', # 中紫红
    '#191970', # 午夜蓝
    '#8B0000', # 深红
]

# --- 辅助函数 ---
def _tight_ra(ra):
    ra_s = np.sort(ra)
    diff_s = np.diff(ra_s)
    if len(diff_s) == 0: return ra
    ind_max = np.argmax(diff_s)
    if diff_s[ind_max] < ra_s[0] + 360 - ra_s[-1]:
        return ra
    else:
        ra = np.copy(ra)
        is_c = ra >= ra_s[ind_max + 1]
        ra[is_c] = ra[is_c] - 360
        return ra

def get_continuous_blocks(mask):
    mask = mask.astype(bool)
    d = np.diff(np.concatenate(([0], mask.view(np.int8), [0])))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    return list(zip(starts, ends))

def process_data_memory_only(fname, args):
    """
    仅计算数据，返回 (Success, Data, Mask)
    """
    basename = os.path.basename(fname)
    print(f"读取: {basename:<30} ... ", end='', flush=True)

    try:
        radec_data = get_radec(fname, nproc=args.nproc, backend=args.backend)
    except Exception as e:
        print(f"[Error] {e}")
        return False, None, None

    if radec_data is None:
        print("[Skip] 数据为空")
        return False, None, None

    # 排序
    ind_sort = np.argsort(radec_data['mjd'])
    for k in radec_data:
        if isinstance(radec_data[k], np.ndarray) and len(radec_data[k]) == len(ind_sort):
            radec_data[k] = radec_data[k][ind_sort]

    dec_ref = radec_data['dec1']
    ra_ref = _tight_ra(radec_data['ra1'])
    
    # 筛选不稳定数据
    window_size = 50
    std_dec = pd.Series(dec_ref).rolling(window=window_size, center=True).std().fillna(np.inf)
    mask_dec_stable = (std_dec < 0.005).values
    
    # 筛选短扫描
    blocks = get_continuous_blocks(mask_dec_stable)
    final_mask = np.zeros_like(mask_dec_stable, dtype=bool)
    
    for start, end in blocks:
        block_ra = ra_ref[start:end]
        if len(block_ra) < 50: continue
        ra_span = np.max(block_ra) - np.min(block_ra)
        if ra_span > 0.1:
            final_mask[start:end] = True

    valid_count = np.sum(final_mask)
    if valid_count == 0:
        print("[Warning] 无有效数据 (扫描过短或不稳定)")
        return False, None, None
    
    print("OK")
    return True, radec_data, final_mask

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot combined RA-DEC tracks (Dark Colors, No White).')
    parser.add_argument('files', nargs='+', help='Input xlsx files')
    parser.add_argument('--outdir', help='Output directory')
    parser.add_argument('--backend', choices=['erfa', 'astropy'], default='astropy', help='backend')
    parser.add_argument('-n', '--nproc', type=int, default=1, help='parallel process number')
    
    args = parser.parse_args()

    # 1. 获取文件列表
    file_list = []
    for f in args.files:
        expanded = glob.glob(f, recursive=True)
        if expanded: file_list.extend(expanded)
        else: file_list.append(f)
    
    file_list = sorted(list(set(file_list)))
    total_files = len(file_list)
    
    if total_files == 0:
        print("未找到文件。")
        sys.exit(0)

    print(f"找到 {total_files} 个文件。")

    # 2. 初始化画布
    fig, ax = plt.subplots(1, 1, figsize=(15, 12), dpi=150)
    
    all_ra_points = []
    all_dec_points = []
    
    plotted_files_count = 0
    
    # 3. 循环处理
    for idx, fname in enumerate(file_list):
        if not (os.path.exists(fname) and fname.endswith('.xlsx')):
            continue
            
        success, data, mask = process_data_memory_only(fname, args)
        
        if success:
            # === 颜色选择逻辑 ===
            # 使用取模运算循环使用深色列表
            # 这样即使有 20 个文件，第 19 个会复用第 1 个颜色，但都保证是深色
            color_idx = plotted_files_count % len(DARK_COLORS)
            current_color = DARK_COLORS[color_idx]
            
            plotted_something = False
            for i in range(1, 20):
                key_ra, key_dec = f'ra{i}', f'dec{i}'
                if key_ra in data:
                    ra = data[key_ra]
                    dec = data[key_dec]
                    ra = _tight_ra(ra)
                    
                    ra_plot = ra.copy()
                    dec_plot = dec.copy()
                    ra_plot[~mask] = np.nan
                    dec_plot[~mask] = np.nan
                    
                    # === 绘图样式调整 ===
                    # linewidth=0.6: 稍微加粗
                    # alpha=0.9: 几乎不透明，颜色更深
                    ax.plot(ra_plot, dec_plot, linewidth=0.6, alpha=0.9, color=current_color)
                    
                    valid_ra = ra[mask]
                    valid_dec = dec[mask]
                    if len(valid_ra) > 0:
                        all_ra_points.append(valid_ra)
                        all_dec_points.append(valid_dec)
                        plotted_something = True
            
            if plotted_something:
                plotted_files_count += 1
    
    print(f"\n处理完毕: 成功绘制了 {plotted_files_count} / {total_files} 个文件的数据。")
    if plotted_files_count < total_files:
        print("提示: 如果有文件未被绘制，请检查上方日志中是否有 [Warning] 或 [Error]。")

    # 4. 设置坐标轴范围
    if len(all_ra_points) > 0:
        full_ra = np.concatenate(all_ra_points)
        full_dec = np.concatenate(all_dec_points)
        
        min_ra, max_ra = np.min(full_ra), np.max(full_ra)
        min_dec, max_dec = np.min(full_dec), np.max(full_dec)
        
        span_ra = max_ra - min_ra
        span_dec = max_dec - min_dec
        
        pad_ra = span_ra * 0.02 if span_ra > 0 else 0.05
        pad_dec = span_dec * 0.02 if span_dec > 0 else 0.05
        
        ax.set_xlim(min_ra - pad_ra, max_ra + pad_ra)
        ax.set_ylim(min_dec - pad_dec, max_dec + pad_dec)
    
    # 5. 样式美化
    ax.set_xlabel('RA (degree)', fontsize=14, fontweight='bold')
    ax.set_ylabel('DEC (degree)', fontsize=14, fontweight='bold')
    ax.set_title('Drift Tracking(all)', fontsize=18, fontweight='bold', pad=15)
    
    ax.minorticks_on()
    # 主网格深一点
    ax.grid(True, which='major', linestyle='-', linewidth=0.7, alpha=0.5, color='gray')
    # 次网格浅一点
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.3, color='gray')
    
    # 6. 保存
    save_dir = args.outdir if args.outdir else '.' 
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    out_pdf = os.path.join(save_dir, "Combined_Tracks_Dark.pdf")
    
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    
    print(f"图片已保存至: {out_pdf}")
