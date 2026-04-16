import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import os
import argparse
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import warnings
import matplotlib.lines as mlines 

def main():
    parser = argparse.ArgumentParser(description="Calculate and plot RMS spatial distribution with dynamic core-slicing.")
    parser.add_argument('--indir', type=str, required=True, help="Input FITS file path")
    parser.add_argument('--outdir', type=str, required=True, help="Output directory path")
    parser.add_argument('--num_cuts', type=int, default=6, help="Number of horizontal/vertical cuts")
    parser.add_argument('--threshold_mult', type=float, default=2.0, help="RMS threshold multiplier based on Mean RMS (默认: 2.0)")
    
    args = parser.parse_args()

    filename = args.indir
    save_path = args.outdir
    num_cuts = args.num_cuts
    threshold_mult = args.threshold_mult
    
    base_name = os.path.splitext(os.path.basename(filename))[0]

    if not os.path.exists(save_path): 
        os.makedirs(save_path)

    print(f"正在分析数据: {filename}")
    with fits.open(filename) as hdul:
        header = hdul[0].header
        data = hdul[0].data  
        if data.ndim == 4: data = data[0]
        wcs = WCS(header).celestial

        print("计算 RMS 空间分布中...")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            rms_map_raw = np.nanstd(data, axis=0)

    # 动态数据分离
    mean_rms = np.nanmean(rms_map_raw)
    threshold_val = mean_rms * threshold_mult
    
    print(f"统计信息:")
    print(f" -> 平均 RMS: {mean_rms:.2e} Jy/beam")
    print(f" -> 设定阈值 ({threshold_mult}x 平均值): {threshold_val:.2e} Jy/beam")
    
    rms_map_sliced = rms_map_raw.copy()
    rms_map_sliced[rms_map_sliced > threshold_val] = np.nan

    print("绘制 RMS 图  ...")
    
    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(2, 2, width_ratios=[4, 1.2], height_ratios=[4, 1.5], hspace=0.05, wspace=0.05)
    
    ax_main = fig.add_subplot(gs[0, 0], projection=wcs) 
    ax_bottom = fig.add_subplot(gs[1, 0])               
    ax_right = fig.add_subplot(gs[0, 1])                

    vmin = np.nanpercentile(rms_map_raw, 2)
    vmax = np.nanpercentile(rms_map_raw, 98)

    # 边界脱钩计算
    # 1. 主图的全景边界 (用于绘图显示全貌，和坐标轴对齐)
    non_nan_mask_raw = ~np.isnan(rms_map_raw)
    rows_raw = np.any(non_nan_mask_raw, axis=1)
    cols_raw = np.any(non_nan_mask_raw, axis=0)
    ymin_raw, ymax_raw = np.where(rows_raw)[0][[0, -1]]
    xmin_raw, xmax_raw = np.where(cols_raw)[0][[0, -1]]

    # 2. 切片的核心边界 (只在青色圈包围的有效范围内取线)
    non_nan_mask_core = ~np.isnan(rms_map_sliced)
    rows_core = np.any(non_nan_mask_core, axis=1)
    cols_core = np.any(non_nan_mask_core, axis=0)
    ymin_core, ymax_core = np.where(rows_core)[0][[0, -1]]
    xmin_core, xmax_core = np.where(cols_core)[0][[0, -1]]
    
    # 使用 core 边界来生成抽样线，保证刀刀切在肉上
    y_cuts = np.linspace(ymin_core + (ymax_core-ymin_core)*0.05, ymax_core - (ymax_core-ymin_core)*0.05, num_cuts).astype(int)
    x_cuts = np.linspace(xmin_core + (xmax_core-xmin_core)*0.05, xmax_core - (xmax_core-xmin_core)*0.05, num_cuts).astype(int)
    cut_colors = plt.cm.tab10(np.linspace(0, 1, num_cuts))
    # ------------------------------------------------

    # 主图 (ax_main)
    im = ax_main.imshow(rms_map_raw, origin='lower', cmap='inferno', 
                        aspect='auto', vmin=vmin, vmax=vmax)

    try:
        ax_main.contour(rms_map_raw, levels=[threshold_val], colors='cyan', 
                        linewidths=2.5, linestyles='--')
        contour_line = mlines.Line2D([], [], color='cyan', linestyle='--', linewidth=2.5, 
                                     label=f'Slice Valid Area (< {threshold_mult}x Mean)')
        ax_main.legend(handles=[contour_line], loc='upper right', facecolor='black', labelcolor='white')
    except Exception as e:
        print(f"提示: 等高线绘制跳过 ({e})")

    # 画参考线 (这些线现在全都在有效区域内了)
    for i, y_cut in enumerate(y_cuts):
        ax_main.axhline(y_cut, color=cut_colors[i], ls='-', alpha=0.8, lw=1.5)
    for i, x_cut in enumerate(x_cuts):
        ax_main.axvline(x_cut, color=cut_colors[i], ls='--', alpha=0.8, lw=1.5)

    lon = ax_main.coords['ra']
    lon.set_ticklabel_visible(False)
    lon.set_axislabel('')
    # 主图依然使用全景边界来显示
    ax_main.set_xlim(xmin_raw, xmax_raw)
    ax_main.set_ylim(ymin_raw, ymax_raw) 

    lat = ax_main.coords['dec']
    lat.set_major_formatter('dd:mm')
    lat.set_ticks(number=8)
    lat.set_axislabel('Declination (J2000)', fontsize=12)

    fig.suptitle(f'RMS Spatial Distribution - {base_name}', fontsize=16, y=0.96)

    # 底部 RA 切片 (ax_bottom)
    # 底部和侧边坐标轴也要使用全景边界 (raw)，这是为了和主图完美对齐
    coord_bl = wcs.pixel_to_world(xmin_raw, ymin_raw)
    coord_tr = wcs.pixel_to_world(xmax_raw, ymax_raw)
    left_ra = coord_bl.ra.degree
    right_ra = coord_tr.ra.degree
    
    x_coords_deg = np.zeros(rms_map_raw.shape[1])
    for x in range(rms_map_raw.shape[1]):
        coord = wcs.pixel_to_world(x, rms_map_raw.shape[0]//2)
        x_coords_deg[x] = coord.ra.degree

    for i, y_cut in enumerate(y_cuts):
        profile = rms_map_sliced[y_cut, :]
        valid = ~np.isnan(profile)
        ax_bottom.plot(x_coords_deg[valid], profile[valid], color=cut_colors[i], ls='-', alpha=0.8, lw=1)

    def ra_formatter_bottom(x, pos):
        val = x / 15.0
        hours = int(val)
        minutes = int((val - hours) * 60)
        seconds = int(((val - hours) * 60 - minutes) * 60)
        return f"{hours}h{abs(minutes):02d}m{abs(seconds):02d}s"

    ax_bottom.xaxis.set_major_formatter(ticker.FuncFormatter(ra_formatter_bottom))
    ax_bottom.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax_bottom.set_xlabel('Right Ascension (J2000)', fontsize=12)
    ax_bottom.set_ylabel('RMS (Jy/beam)', fontsize=12)
    ax_bottom.set_xlim(left_ra, right_ra) 
    
    ax_bottom.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax_bottom.grid(color='grey', ls='dotted', alpha=0.4)

    # 右侧 Dec 切片 (ax_right)
    y_coords_deg = np.zeros(rms_map_raw.shape[0])
    for y in range(rms_map_raw.shape[0]):
        coord = wcs.pixel_to_world(rms_map_raw.shape[1]//2, y)
        y_coords_deg[y] = coord.dec.degree

    for i, x_cut in enumerate(x_cuts):
        profile = rms_map_sliced[:, x_cut]
        valid = ~np.isnan(profile)
        ax_right.plot(profile[valid], y_coords_deg[valid], color=cut_colors[i], ls='--', alpha=0.8, lw=1)

    ax_right.set_xlabel('RMS (Jy/beam)', fontsize=12)
    bottom_dec = coord_bl.dec.degree
    top_dec = coord_tr.dec.degree
    ax_right.set_ylim(bottom_dec, top_dec)
    
    def dec_formatter_right(x, pos):
        deg = int(x)
        min_val = int(abs(x - deg) * 60)
        return f"{deg}°{min_val:02d}'"
    
    ax_right.yaxis.set_major_formatter(ticker.FuncFormatter(dec_formatter_right))
    ax_right.yaxis.tick_right()
    ax_right.yaxis.set_label_position("right")
    ax_right.set_ylabel('Declination (J2000)', fontsize=12)
    
    ax_right.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    
    ax_right.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax_right.grid(color='grey', ls='dotted', alpha=0.4)

    # Colorbar 安置
    plt.subplots_adjust(left=0.08, right=0.82, top=0.92, bottom=0.1)

    cax = fig.add_axes([0.91, 0.1, 0.015, 0.82]) 
    cbar = plt.colorbar(im, cax=cax)
    
    cbar.set_label('RMS Noise (Jy/beam)', fontsize=12)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()

    # 保存 
    save_file = os.path.join(save_path, f'{base_name}_rms_map.png')
    plt.savefig(save_file, dpi=1000, bbox_inches='tight')
    save_file_pdf = os.path.join(save_path, f'{base_name}_rms_map.pdf')
    plt.savefig(save_file_pdf, dpi=1000, bbox_inches='tight')
    print(f"\n结果图已保存为:\n{save_file}\n{save_file_pdf}")

if __name__ == "__main__":
    main()
