import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import argparse

def main():
    # 1. 设置 argparse 参数解析
    parser = argparse.ArgumentParser(description="Calculate and plot RMS noise per channel for FITS data cube.")
    parser.add_argument('--indir', type=str, required=True, help="Input FITS file path (输入FITS文件的完整路径)")
    parser.add_argument('--outdir', type=str, required=True, help="Output directory path (输出图表和结果的目录)")
    parser.add_argument('--threshold', type=float, default=3.0, help="RMS outlier threshold multiplier (异常值判定的标准差倍数，默认: 3.0)")

    args = parser.parse_args()

    filename = args.indir
    save_path = args.outdir
    threshold_multiplier = args.threshold

    # 自动提取输入文件的基本名称（不带 .fits 后缀），用于后续图片命名
    base_name = os.path.splitext(os.path.basename(filename))[0]

    if not os.path.exists(save_path): 
        os.makedirs(save_path)

    if not os.path.exists(filename):
        print(f"找不到文件: {filename}")
        exit()

    print(f"正在分析数据: {filename}")
    with fits.open(filename) as hdul:
        header = hdul[0].header
        data = hdul[0].data  # 维度 [freq, dec, ra]

        # FITS 公式: Freq = CRVAL3 + (index + 1 - CRPIX3) * CDELT3
        crval3 = header.get('CRVAL3')
        cdelt3 = header.get('CDELT3')
        crpix3 = header.get('CRPIX3', 1)  # 如果没写，默认是 1
        
        # 计算每个频道对应的频率 (单位转为 MHz)
        n_channels = data.shape[0]
        indices = np.arange(n_channels)
        freqs_mhz = (crval3 + (indices + 1 - crpix3) * cdelt3) / 1e6
        
        # 计算每个频道的RMS
        print("计算 RMS 中...")
        rms_array = np.sqrt(np.nanmean(data**2, axis=(1, 2)))

    # 2. 识别“偏离太大”的频率
    median_rms = np.nanmedian(rms_array)
    std_rms = np.nanstd(rms_array)
    threshold_val = threshold_multiplier * std_rms # 动态调整倍数

    # 找出异常点的索引
    bad_indices = np.where(np.abs(rms_array - median_rms) > threshold_val)[0]

    # 3. 打印异常频率列表
    print("\n" + "="*60)
    print(f"{'频道':<10} | {'频率 (MHz)':<15} | {'RMS 值 (Jy/beam)':<20} | {'偏离程度'}")
    print("-" * 60)

    for idx in bad_indices:
        diff = (rms_array[idx] - median_rms) / median_rms * 100
        print(f"{idx:<10} | {freqs_mhz[idx]:<15.4f} | {rms_array[idx]:<20.6e} | {diff:>+7.1f}%")

    print("="*60)
    print(f"统计：在 {n_channels} 个频道中，发现了 {len(bad_indices)} 个异常频率。")

    # 4. 绘图
    plt.figure(figsize=(12, 6))

    # 绘制 RMS 曲线
    plt.plot(freqs_mhz, rms_array, color='blue', linewidth=0.7, label='RMS per Channel')

    # 绘制中位数参考线
    plt.axhline(y=median_rms, color='red', linestyle='--', label=f'Median RMS: {median_rms:.2e}')

    # 设置 X 轴和 Y 轴的标签解释
    plt.xlabel('Frequency (MHz)', fontsize=12)
    plt.ylabel('RMS Noise (Jy/beam)', fontsize=12)
    plt.title(f'RMS Distribution - {base_name}', fontsize=14)

    # 强制 Y 轴使用科学计数法
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # 显示图例
    plt.legend(loc='upper right')

    # 添加网格线，方便观察
    plt.grid(True, linestyle=':', alpha=0.6)

    # 自动调整布局，防止标签被裁剪
    plt.tight_layout()

    # 保存图片 (使用动态提取的文件名)
    save_file = os.path.join(save_path, f'{base_name}_rms.png')
    plt.savefig(save_file, dpi=300)
    print(f"\n结果图已保存为: {save_file}")

if __name__ == "__main__":
    main()
