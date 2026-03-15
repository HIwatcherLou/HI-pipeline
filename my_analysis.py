import argparse
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import hifast.funcs as hf
import warnings

def plot_data(fpath, output_dir):
    """
    读取 HDF5 文件，自动寻找有效数据最多的连续 20 个积分时间，并画出平均频谱图
    """
    polar = 0
    slice_length = 20  # 我们想要的连续时间长度（20行）
    
    filename = os.path.basename(fpath)
    print(f"正在处理: {filename} ...")
    
    try:
        # 1. 读取全部数据
        S = hf.HFSpec(fpath)
        x = S['freq'][:] 
        full_data = S['DATA'][:, :, polar] # 获取该极化的全部 2D 数据 (时间 x 频率)
        
        # 2. 获取总时间长度
        total_time_steps = full_data.shape[0]
        
        # 3. 核心改进：动态寻找最优切片
        if total_time_steps < slice_length:
            # 如果文件本身的时间还不到 20 个，就直接全取
            i_start, i_end = 0, total_time_steps
        else:
            # 计算每一个时间点（每一行）里有多少个非 NaN 的有效数据
            valid_counts_per_time = np.sum(~np.isnan(full_data), axis=1)
            
            # 使用卷积计算“滑动窗口”内的有效数据总和
            # 这行代码的意思是：把连续 20 行的有效数据个数加起来
            window_sums = np.convolve(valid_counts_per_time, np.ones(slice_length), mode='valid')
            
            # 找到总和最大的那个窗口的起始位置
            best_start = int(np.argmax(window_sums))
            
            i_start = best_start
            i_end = best_start + slice_length

        print(f"  -> 自动选定最优时间切片: {i_start} 到 {i_end}")

        # 4. 提取出这最好的一段，开始压扁（求平均）
        best_slice = full_data[i_start:i_end, :]
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            y = np.nanmean(best_slice, axis=0, dtype=np.float64)

        # 5. 开始画图
        plt.figure(figsize=(10, 4))
        
        if np.isnan(y).all():
            print(f"  -> [警告] {filename} 整个文件几乎全为 NaN。")
            plt.text(0.5, 0.5, 'All Data Masked (NaN)', 
                     horizontalalignment='center', verticalalignment='center', 
                     transform=plt.gca().transAxes, color='red', fontsize=14, fontweight='bold')
            plt.plot(x, np.zeros_like(x), alpha=0) 
        else:
            # 正常画图
            plt.plot(x, y, label=f'Mean Spectrum ({i_start}-{i_end})', alpha=0.7)
            
            # 填补个别残余 NaN 后进行高斯平滑
            y_filled = np.nan_to_num(y, nan=np.nanmean(y))
            plt.plot(x, ndimage.gaussian_filter1d(y_filled, 3), label='Gaussian Smoothed', color='red')
            
            # 设定 Y 轴
            y_min, y_max = np.nanpercentile(y, 0), np.nanpercentile(y, 98)
            if not np.isnan(y_min) and not np.isnan(y_max) and (y_min < y_max):
                plt.ylim(y_min, y_max)

        # 标题中带上自动选取的时间段，方便检查
        plt.title(f"{filename}\n(Auto-selected slice: {i_start} to {i_end})", fontsize=10)
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Intensity')
        plt.minorticks_on()
        plt.legend() if not np.isnan(y).all() else None
        plt.grid(True, linestyle=':', alpha=0.6)
        
        # 替换后缀并保存
        plot_filename = filename.replace('.hdf5', '.png')
        plot_path = os.path.join(output_dir, plot_filename)
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"  -> 成功保存: {plot_path}")

    except Exception as e:
        print(f"  -> [错误] 处理 {filename} 时失败: {e}")

    finally:
        plt.close() 

def main():
    parser = argparse.ArgumentParser(description="批量绘制指定文件夹下所有 HDF5 最终数据的频谱图")
    parser.add_argument('--indir', type=str, required=True, help='输入文件夹路径')
    parser.add_argument('--outdir', type=str, required=True, help='输出文件夹路径')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        print(f"已创建输出目录: {args.outdir}")
    
    search_pattern = os.path.join(args.indir, '*.hdf5')
    hdf5_files = sorted(glob.glob(search_pattern))

    if len(hdf5_files) == 0:
        print(f"警告: 在文件夹 '{args.indir}' 中没有找到任何 .hdf5 文件！请检查路径。")
        return 

    print(f"\n共找到 {len(hdf5_files)} 个 HDF5 文件，准备开始批量绘图...\n")
    print("-" * 50)

    for fpath in hdf5_files:
        plot_data(fpath, args.outdir)

    print("-" * 50)
    print("所有文件处理完毕！")

if __name__ == "__main__":
    main()
