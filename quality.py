import argparse
import h5py
import numpy as np
import os
import glob
import warnings
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

def masked_boxcar_smooth_1d(values, width):
    """
    NaN 感知的 1D 箱型平滑。

    逻辑：
    - 分子卷积：仅累积有限值；
    - 分母卷积：累积有效权重；
    - 输出=分子/分母，分母为 0 时保留 NaN。
    """
    if width <= 1:
        return values

    y = np.asarray(values, dtype=np.float64)
    good = np.isfinite(y)
    kernel = np.ones(int(width), dtype=np.float64) / float(width)

    num = np.convolve(np.where(good, y, 0.0), kernel, mode='same')
    den = np.convolve(good.astype(np.float64), kernel, mode='same')

    out = np.full_like(y, np.nan, dtype=np.float64)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def read_hdf5_data(filepath, pol_stat='median', chunk_size=256):
    """
    读取 FAST HDF5 数据，按共享 S/is_rfi 手动设置 NaN，并分块生成 2D 瀑布图。

    处理策略：
    - 沿极化轴聚合为 (time, channel)；
    - 默认使用 nanmedian 提升对离群值/RFI 边界的稳健性；
    - 按时间分块读取，避免一次性载入超大数组。
    """
    with h5py.File(filepath, 'r') as f:
        if 'S/flux' not in f:
            raise ValueError(f"文件中缺少 'S/flux' 数据集: {filepath}")

        flux_ds = f['S/flux']
        n_pol, n_time, n_chan = flux_ds.shape
        data_2d = np.full((n_time, n_chan), np.nan, dtype=np.float32)

        if 'S/is_rfi' in f:
            rfi_ds = f['S/is_rfi']
            print("  [Info] 检测到共享 S/is_rfi，按分块方式应用掩膜。")
        else:
            rfi_ds = None
            print(f"  [Warning] 文件中未找到 'S/is_rfi'，跳过掩模处理: {filepath}")

        step = max(1, int(chunk_size))
        for start in range(0, n_time, step):
            end = min(n_time, start + step)
            flux_chunk = np.asarray(flux_ds[:, start:end, :], dtype=np.float32)

            if rfi_ds is not None:
                rfi_chunk = np.asarray(rfi_ds[start:end, :], dtype=bool)
                rfi_mask_3d = np.broadcast_to(rfi_chunk[None, :, :], flux_chunk.shape)
                flux_chunk[rfi_mask_3d] = np.nan

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                if pol_stat == 'mean':
                    data_2d[start:end, :] = np.nanmean(flux_chunk, axis=0).astype(np.float32)
                else:
                    data_2d[start:end, :] = np.nanmedian(flux_chunk, axis=0).astype(np.float32)

        if 'S/freq' in f:
            freq = f['S/freq'][:]
        else:
            naxis1 = data_2d.shape[1]
            crval1 = 1.300006628036E+03
            cdelt1 = 7.629394531250E-03
            freq = crval1 + np.arange(naxis1) * cdelt1
            
    return data_2d, freq

def plot_comparison(
    file1_path,
    file2_path,
    output_dir,
    suffix1,
    title1,
    title2,
    main_title,
    smooth_width=1,
    stat='median',
    chunk_size=256,
):
    basename = os.path.basename(file1_path).replace(suffix1, '')
    print(f"正在绘制并保存: {basename} ...")
    
    data1, freq1 = read_hdf5_data(file1_path, pol_stat=stat, chunk_size=chunk_size)
    data2, freq2 = read_hdf5_data(file2_path, pol_stat=stat, chunk_size=chunk_size)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if stat == 'mean':
            spectrum1 = np.nanmean(data1, axis=0)
            spectrum2 = np.nanmean(data2, axis=0)
        else:
            spectrum1 = np.nanmedian(data1, axis=0)
            spectrum2 = np.nanmedian(data2, axis=0)

    # 频率方向平滑：使用 NaN 感知算法，避免 NaN 传播造成伪结构。
    if smooth_width > 1:
        spectrum1 = masked_boxcar_smooth_1d(spectrum1, smooth_width)
        spectrum2 = masked_boxcar_smooth_1d(spectrum2, smooth_width)

    cmap = plt.get_cmap('hot').copy()
    cmap.set_bad(color='white') 

    # 两幅瀑布图使用统一色标，保证前后对比可解释。
    vals1 = data1[np.isfinite(data1)]
    vals2 = data2[np.isfinite(data2)]
    if vals1.size > 0 and vals2.size > 0:
        vals_all = np.concatenate([vals1, vals2])
    elif vals1.size > 0:
        vals_all = vals1
    else:
        vals_all = vals2
    if vals_all.size > 0:
        vmin, vmax = np.nanpercentile(vals_all, [5, 95])
        if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or (vmin >= vmax):
            vmin, vmax = np.nanmin(vals_all), np.nanmax(vals_all)
    else:
        vmin, vmax = -1.0, 1.0

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"{main_title}: {basename}", fontsize=16)
    
    im1 = axes[0, 0].imshow(data1, aspect='auto', origin='lower', cmap=cmap,
                            extent=[freq1[0], freq1[-1], 0, data1.shape[0]],
                            vmin=vmin, vmax=vmax)
    axes[0, 0].set_title(f'{title1} - Waterfall')
    axes[0, 0].set_ylabel('Time (Integration)')
    axes[0, 0].set_xlabel('Frequency (MHz)')
    fig.colorbar(im1, ax=axes[0, 0])

    im2 = axes[0, 1].imshow(data2, aspect='auto', origin='lower', cmap=cmap,
                            extent=[freq2[0], freq2[-1], 0, data2.shape[0]],
                            vmin=vmin, vmax=vmax)
    axes[0, 1].set_title(f'{title2} - Waterfall')
    axes[0, 1].set_ylabel('Time (Integration)')
    axes[0, 1].set_xlabel('Frequency (MHz)')
    fig.colorbar(im2, ax=axes[0, 1])

    axes[1, 0].plot(freq1, spectrum1)
    axes[1, 0].axhline(0.0, color='k', lw=0.8, ls='--', alpha=0.7)
    axes[1, 0].set_title(f'{title1} - Averaged Spectrum')
    axes[1, 0].set_ylabel('Intensity (T-flux)')
    axes[1, 0].set_xlabel('Frequency (MHz)')

    axes[1, 1].plot(freq2, spectrum2)
    axes[1, 1].axhline(0.0, color='k', lw=0.8, ls='--', alpha=0.7)
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
    parser.add_argument('--stat', choices=['median', 'mean'], default='median',
                        help='统计口径：median 更稳健；mean 与旧版本一致')
    parser.add_argument('--chunk-size', type=int, default=256,
                        help='按时间轴分块读取大小，降低内存占用')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    files_stage1 = glob.glob(args.file_pattern)
    if not files_stage1:
        print(f"Error: 找不到文件 {args.file_pattern}")
        return
    
    files_stage1.sort()

    print(f"共找到 {len(files_stage1)} 个待处理文件。平滑宽度: {args.smooth}，统计口径: {args.stat}，chunk-size: {args.chunk_size}")

    for file1_path in files_stage1:
        file2_path = file1_path.replace(args.suffix1, args.suffix2)
        if not os.path.exists(file2_path):
            print(f"Skip: 找不到对照文件 {file2_path}")
            continue
        # 将 smooth 参数传入绘图函数
        plot_comparison(
            file1_path=file1_path,
            file2_path=file2_path,
            output_dir=args.output_dir,
            suffix1=args.suffix1,
            title1=args.title1,
            title2=args.title2,
            main_title=args.maintitle,
            smooth_width=args.smooth,
            stat=args.stat,
            chunk_size=args.chunk_size,
        )

if __name__ == "__main__":
    main()

