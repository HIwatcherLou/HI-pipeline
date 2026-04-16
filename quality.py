import argparse
import h5py
import numpy as np
import os
import glob
import warnings
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

def read_hdf5_data(filepath):
    """
    读取 FAST HDF5 数据，并根据 S/is_rfi 手动设置 NaN
    """
    with h5py.File(filepath, 'r') as f:
        if 'S/flux' not in f:
            raise ValueError(f"文件中缺少 'S/flux' 数据集: {filepath}")
        
        flux = f['S/flux'][:].astype(np.float32)
        
        if 'S/is_rfi' in f:
            rfi_mask = f['S/is_rfi'][:]
            flux[:, rfi_mask > 0] = np.nan
            print(f"  [Info] 已根据 S/is_rfi 掩盖了 {np.sum(rfi_mask)} 个数据点")
        else:
            print(f"  [Warning] 文件中未找到 'S/is_rfi'，跳过掩模处理: {filepath}")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            data_2d = np.nanmean(flux, axis=0) 

        if 'S/freq' in f:
            freq = f['S/freq'][:]
        else:
            naxis1 = data_2d.shape[1]
            crval1 = 1.300006628036E+03
            cdelt1 = 7.629394531250E-03
            freq = crval1 + np.arange(naxis1) * cdelt1
            
    return data_2d, freq

def plot_comparison(file1_path, file2_path, output_dir, suffix1, title1, title2, main_title, smooth_width=1):
    basename = os.path.basename(file1_path).replace(suffix1, '')
    print(f"正在绘制并保存: {basename} ...")
    
    data1, freq1 = read_hdf5_data(file1_path)
    data2, freq2 = read_hdf5_data(file2_path)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        spectrum1 = np.nanmean(data1, axis=0)
        spectrum2 = np.nanmean(data2, axis=0)

    # --- 新增功能：频率平滑 (Channel Smoothing) ---
    if smooth_width > 1:
        kernel = np.ones(smooth_width) / smooth_width
        # 使用 mode='same' 保持长度一致，注意 NaN 在卷积后会扩散（这是天文处理中的正常现象）
        spectrum1 = np.convolve(spectrum1, kernel, mode='same')
        spectrum2 = np.convolve(spectrum2, kernel, mode='same')
    # --------------------------------------------

    cmap = plt.get_cmap('hot').copy()
    cmap.set_bad(color='white') 

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"{main_title}: {basename}", fontsize=16)
    
    im1 = axes[0, 0].imshow(data1, aspect='auto', origin='lower', cmap=cmap,
                            extent=[freq1[0], freq1[-1], 0, data1.shape[0]])
    axes[0, 0].set_title(f'{title1} - Waterfall')
    axes[0, 0].set_ylabel('Time (Integration)')
    axes[0, 0].set_xlabel('Frequency (MHz)')
    fig.colorbar(im1, ax=axes[0, 0])

    im2 = axes[0, 1].imshow(data2, aspect='auto', origin='lower', cmap=cmap,
                            extent=[freq2[0], freq2[-1], 0, data2.shape[0]])
    axes[0, 1].set_title(f'{title2} - Waterfall')
    axes[0, 1].set_ylabel('Time (Integration)')
    axes[0, 1].set_xlabel('Frequency (MHz)')
    fig.colorbar(im2, ax=axes[0, 1])

    axes[1, 0].plot(freq1, spectrum1)
    axes[1, 0].set_title(f'{title1} - Averaged Spectrum')
    axes[1, 0].set_ylabel('Intensity (T-flux)')
    axes[1, 0].set_xlabel('Frequency (MHz)')

    axes[1, 1].plot(freq2, spectrum2)
    axes[1, 1].set_title(f'{title2} - Averaged Spectrum')
    axes[1, 1].set_ylabel('Intensity (T-flux)')
    axes[1, 1].set_xlabel('Frequency (MHz)')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_filename = os.path.join(output_dir, f"{basename}_{title1}_vs_{title2}.png".replace(' ', '_'))
    plt.savefig(output_filename, dpi=300)
    plt.close(fig) 

def main():
    parser = argparse.ArgumentParser(description="FAST 数据处理流程比对绘图工具 (修复 RFI Mask + 频率平滑)")
    parser.add_argument('-f', '--file_pattern', type=str, required=True, help='文件通配符路径')
    parser.add_argument('-o', '--output_dir', type=str, default='./', help='输出目录')
    parser.add_argument('--suffix1', type=str, default='_T-flux.hdf5', help='阶段1后缀')
    parser.add_argument('--suffix2', type=str, default='_T-flux-bld_p.hdf5', help='阶段2后缀')
    parser.add_argument('--title1', type=str, default='Before RFI', help='标题1')
    parser.add_argument('--title2', type=str, default='After RFI', help='标题2')
    parser.add_argument('--maintitle', type=str, default='Pipeline Comparison', help='主标题')
    
    # --- 新增参数：平滑宽度 ---
    parser.add_argument('--smooth', type=int, default=1, help='频率平滑宽度 (channel 数量，建议使用奇数)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    files_stage1 = glob.glob(args.file_pattern)
    if not files_stage1:
        print(f"Error: 找不到文件 {args.file_pattern}")
        return
    
    files_stage1.sort()

    print(f"共找到 {len(files_stage1)} 个待处理文件。平滑宽度: {args.smooth}")

    for file1_path in files_stage1:
        file2_path = file1_path.replace(args.suffix1, args.suffix2)
        if not os.path.exists(file2_path):
            print(f"Skip: 找不到对照文件 {file2_path}")
            continue
        # 将 smooth 参数传入绘图函数
        plot_comparison(file1_path, file2_path, args.output_dir, args.suffix1, args.title1, args.title2, args.maintitle, smooth_width=args.smooth)

if __name__ == "__main__":
    main()
