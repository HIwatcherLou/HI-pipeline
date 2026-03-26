import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import os
import argparse
import matplotlib.ticker as ticker # 导入刻度管理模块，用于进一步优化刻度

def main():
    # 1. 设置 argparse 参数解析
    parser = argparse.ArgumentParser(description="Calculate and plot RMS spatial distribution for FITS data cube.")
    parser.add_argument('--indir', type=str, required=True, help="Input FITS file path (输入FITS文件的完整路径)")
    parser.add_argument('--outdir', type=str, required=True, help="Output directory path (输出图表和结果的目录)")
    
    args = parser.parse_args()

    filename = args.indir
    save_path = args.outdir

    # 自动提取输入文件的基本名称
    base_name = os.path.splitext(os.path.basename(filename))[0]

    if not os.path.exists(save_path): 
        os.makedirs(save_path)

    if not os.path.exists(filename):
        print(f"找不到文件: {filename}")
        exit()

    print(f"正在分析数据: {filename}")
    with fits.open(filename) as hdul:
        header = hdul[0].header
        data = hdul[0].data  
    
        # 如果是 4D 数据，我们去掉第一个维度，变成 [freq, dec, ra]
        if data.ndim == 4:
            data = data[0]

        # 提取 2D 空间坐标信息 (忽略频率轴)，用于在图上画出真实的 RA 和 DEC
        wcs = WCS(header).celestial

        # 2. 计算 RMS 空间分布
        print("计算 RMS 空间分布中...")
        # 沿着频率轴 (axis=0) 计算标准差。
        rms_map = np.nanstd(data, axis=0)

    # 3. 绘图 
    print("绘制 RMS 空间分布图...")
    
    # (横向, 纵向)
    plt.figure(figsize=(12, 4)) 

    # 使用 WCS 投影，让坐标轴直接显示赤经和赤纬
    ax = plt.subplot(projection=wcs)

    # 3.2 增加对比度：手动设定色彩范围 (使用分位数)
    vmin = np.nanpercentile(rms_map, 2)
    vmax = np.nanpercentile(rms_map, 98)
    print(f"设定色彩范围 (2%-98%分位数): vmin={vmin:.2e}, vmax={vmax:.2e}")

    im = ax.imshow(rms_map, origin='lower', cmap='inferno', 
                   aspect='auto', vmin=vmin, vmax=vmax)

    # 3.4 调整刻度定位器和格式化器：
    # 拉伸后，Dec 的刻度可能会变得混乱。可以手动控制它。
    dec = ax.coords['dec']
    dec.set_major_formatter('dd:mm') # 例如，设置为 度.分 格式
    dec.set_ticks(number=10) # 强制显示大约 10 个刻度，让拉伸区域有参考

    # 3.5 添加 Colorbar (颜色条) 并设置科学计数法
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label('RMS Noise (Jy/beam)', fontsize=12)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()

    # 3.6 设置坐标轴标签和标题
    ax.set_xlabel('Right Ascension (J2000)', fontsize=12)
    ax.set_ylabel('Declination (J2000)', fontsize=12)
    plt.title(f'RMS Spatial Distribution (Stretched) - {base_name}', fontsize=14)

    # 3.7 添加坐标网格 - 调整透明度
    ax.grid(color='white', ls='dotted', alpha=0.3) # 降低网格透明度，防止遮挡细节

    # 3.8 自动调整布局 - bbox_inches='tight' 很重要
    plt.tight_layout()

    # 4. 保存图片 - 增加 dpi 和 bbox_inches
    save_file = os.path.join(save_path, f'{base_name}_rms_map.png')
    plt.savefig(save_file, dpi=1000, bbox_inches='tight')
    save_file_pdf = os.path.join(save_path, f'{base_name}_rms_map.pdf')
    plt.savefig(save_file_pdf, dpi=1000,bbox_inches='tight')
    print(f"\n结果图已保存为: {save_file},{save_file_pdf}")
    
    #保存fits文件，可以选择关闭（注释掉就行）
    fits.writeto(os.path.join(save_path, f'{base_name}_rms_map.fits'), rms_map, wcs.to_header(), overwrite=True)

if __name__ == "__main__":
    main()
