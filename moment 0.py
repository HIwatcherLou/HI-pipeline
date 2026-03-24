import matplotlib
matplotlib.use('Agg')  # 强制使用无头模式，解决服务器报错

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
from astropy.io import fits
from astropy.wcs import WCS
import copy
import warnings
import os
import argparse 
import sys

warnings.filterwarnings('ignore')

class CustomExpNorm(mcolors.Normalize):
    def __init__(self, a, vmin=None, vmax=None, clip=False):
        self.a = a
        super().__init__(vmin, vmax, clip)
    def __call__(self, value, clip=None):
        x = super().__call__(value, clip)
        if np.abs(self.a - 1.0) < 1e-6:
            return x
        return (np.ma.power(self.a, x) - 1) / (self.a - 1)

def get_args():
    parser = argparse.ArgumentParser(description="生成 Moment 0 图像，支持指定频率/通道范围、视觉拉伸及自定义配色")
    
    # 输入与输出路径配置
    parser.add_argument('-i', '--indir', type=str, required=True, help='输入 FITS 文件路径')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='输出目录路径')
    
    # 颜色与背景配置
    parser.add_argument('--cmap', type=str, default='inferno', help='颜色映射 (默认: inferno，可选 viridis, plasma, gray 等)')
    parser.add_argument('--bg', type=str, default='#FFFFFF', help='NaN/无效数据的背景颜色 (默认: #FFFFFF 白色)')

    # 频率范围参数 (单位 MHz)
    parser.add_argument('--freq', type=float, nargs=2, metavar=('START', 'END'), 
                        help='指定频率范围 (单位: MHz)，例如 --freq 1420.1 1420.8')
    
    # 通道范围参数
    parser.add_argument('-c', '--channel', type=int, nargs=2, metavar=('START', 'END'),
                        help='指定通道索引范围 (整数)，例如 --channel 100 200')

    # 视觉拉伸参数组
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--log', action='store_true', help='使用对数 (Log) 视觉拉伸')
    group.add_argument('--sqrt', action='store_true', help='使用开方 (Sqrt) 视觉拉伸')
    group.add_argument('--square', action='store_true', help='使用平方 (Square) 视觉拉伸')
    parser.add_argument('--power', type=float, help='使用指数拉伸 a。')
    
    return parser.parse_args()

def run_clean_plot():
    args = get_args()
    input_fits_file = args.indir
    output_dir = args.outdir
    
    # 1. 读取数据
    try:
        with fits.open(input_fits_file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs_3d = WCS(header)
            
            # 打印基础信息，帮助用户判断通道范围
            print(f"成功读取文件: {input_fits_file}")
            print(f"数据形状 (Shape): {data.shape}  [格式: (光谱, Dec, RA)]")
            print(f"总通道数: {data.shape[0]}")
            
    except FileNotFoundError:
        print(f"错误: 找不到文件 {input_fits_file}")
        return

    # 2. 处理切片 (优先检查 Channel，再检查 Freq)
    channel_slice = slice(None) # 默认全选
    freq_suffix = ""
    
    # --- 分支 A: 使用通道索引 (最优先) ---
    if args.channel:
        c_start, c_end = sorted(args.channel)
        
        # 边界检查
        max_chan = data.shape[0] - 1
        c_start = max(0, min(c_start, max_chan))
        c_end = max(0, min(c_end, max_chan))
        
        channel_slice = slice(c_start, c_end + 1)
        freq_suffix = f"_chan{c_start}-{c_end}"
        
        print(f"通道过滤: 索引 {c_start} 到 {c_end} (共 {c_end - c_start + 1} 个通道)")

    # --- 分支 B: 使用频率范围 ---
    elif args.freq:
        f_start_mhz, f_end_mhz = args.freq
        f_start_hz = f_start_mhz * 1e6
        f_end_hz = f_end_mhz * 1e6
        
        try:
            wcs_spec = wcs_3d.sub(['spectral'])
            if wcs_spec.naxis == 0:
                raise ValueError("无法在 WCS 中找到光谱轴")

            ctype = wcs_spec.wcs.ctype[0]
            is_velocity = 'VEL' in ctype or 'VRAD' in ctype or 'FELO' in ctype
            
            val_start, val_end = f_start_hz, f_end_hz
            
            if is_velocity:
                restfrq = header.get('RESTFRQ', header.get('CRVAL3'))
                if not restfrq: restfrq = wcs_3d.wcs.restfrq
                
                if restfrq:
                    c_ms = 299792458.0
                    val_start = c_ms * (1 - f_start_hz / restfrq)
                    val_end = c_ms * (1 - f_end_hz / restfrq)
                    print(f"  [转换] 检测到速度轴，已将 MHz 转换为速度进行查找。")

            pix_start = float(wcs_spec.all_world2pix([val_start], 0)[0])
            pix_end = float(wcs_spec.all_world2pix([val_end], 0)[0])
            
            idx1, idx2 = sorted([int(round(pix_start)), int(round(pix_end))])
            
            max_chan = data.shape[0] - 1
            idx1 = max(0, min(idx1, max_chan))
            idx2 = max(0, min(idx2, max_chan))
            
            channel_slice = slice(idx1, idx2 + 1)
            freq_suffix = f"_{f_start_mhz:.2f}-{f_end_mhz:.2f}MHz"
            
            print(f"频率过滤: {f_start_mhz} MHz -> {f_end_mhz} MHz")
            print(f"对应通道: {idx1} 到 {idx2}")
            
        except Exception as e:
            print(f"警告: 频率转换失败 ({e})，将使用全频率范围。")
            import traceback
            traceback.print_exc()

    # 3. 自动生成输出路径
    base_name = os.path.basename(input_fits_file)
    name_no_ext = os.path.splitext(base_name)[0]
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output_png_file = os.path.join(output_dir, f"{name_no_ext}{freq_suffix}.png")
    output_fits_file = os.path.join(output_dir, f"{name_no_ext}{freq_suffix}.fits")

    # 4. 物理计算 (Moment 0)
    try:
        rest_freq = header.get('RESTFRQ', header.get('CRVAL3'))
        delta_freq = header.get('CDELT3')
        if delta_freq is None and 'CD3_3' in header: delta_freq = header['CD3_3']
            
        c_kms = 299792.458
        if rest_freq and delta_freq:
            delta_v = abs(delta_freq / rest_freq) * c_kms
        else:
            delta_v = 1.0
            
        original_unit = header.get('BUNIT', 'Jy/beam')
        unit_str = f"{original_unit} km/s"
    except Exception:
        delta_v = 1.0
        original_unit = 'Jy/beam'
        unit_str = "Integrated Intensity"

    # 执行切片和累加
    subset_data = data[channel_slice, :, :]
    
    if subset_data.shape[0] == 0:
        print("错误: 选定的范围内没有数据 (通道数=0)。")
        return

    moment0_data = np.nansum(subset_data, axis=0) * delta_v
    
    all_nans = np.all(np.isnan(subset_data), axis=0)
    moment0_data[all_nans] = np.nan

    # 5. 保存 FITS
    new_header = header.copy()
    new_header['NAXIS'] = 2
    for key in ['NAXIS3', 'CRVAL3', 'CDELT3', 'CRPIX3', 'CTYPE3', 'CUNIT3',
                'PC1_3', 'PC2_3', 'PC3_1', 'PC3_2', 'PC3_3']:
        if key in new_header: del new_header[key]
    
    new_header['BUNIT'] = unit_str
    if args.channel:
        new_header['HISTORY'] = f"Channel range: {args.channel[0]}-{args.channel[1]}"
    if args.freq:
        new_header['HISTORY'] = f"Freq range: {args.freq[0]}-{args.freq[1]} MHz"
    
    print(f"正在保存 FITS 到: {output_fits_file}")
    fits.writeto(output_fits_file, moment0_data, new_header, overwrite=True)

    # 6. 绘图准备
    valid_data = moment0_data[~np.isnan(moment0_data)]
    if len(valid_data) > 0:
        vmin = np.percentile(valid_data, 0.2)
        vmax = np.percentile(valid_data, 99.8)
    else:
        vmin, vmax = 0, 1

    norm = None
    scale_label = "Linear"
    if args.log:
        scale_label = "Log"
        log_vmin = max(vmin, 1e-4 * vmax) if vmax > 0 else 1e-4
        norm = mcolors.LogNorm(vmin=log_vmin, vmax=vmax)
    elif args.sqrt:
        scale_label = "Sqrt"
        norm = mcolors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax)
    elif args.square:
        scale_label = "Square"
        norm = mcolors.PowerNorm(gamma=2.0, vmin=vmin, vmax=vmax)
    elif args.power is not None:
        scale_label = f"Exp (a={args.power})"
        norm = CustomExpNorm(a=args.power, vmin=vmin, vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # 7. 绘图
    wcs_2d = wcs_3d.celestial
    plt.style.use('default') 
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs_2d)

    # 动态应用用户设置的颜色映射和背景色
    cmap = copy.copy(plt.get_cmap(args.cmap))
    cmap.set_bad(color=args.bg) 

    im = ax.imshow(moment0_data, cmap=cmap, norm=norm, origin='lower', interpolation='nearest')

    try:
        lon, lat = ax.coords[0], ax.coords[1]
        lon.set_major_formatter('hh:mm:ss'); lat.set_major_formatter('dd:mm:ss')
        lon.set_axislabel('Right Ascension (J2000)', fontsize=12)
        lat.set_axislabel('Declination (J2000)', fontsize=12)
    except Exception:
        pass

    cax = ax.inset_axes([1.06, 0, 0.03, 1], transform=ax.transAxes)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(f"{unit_str} [{scale_label}]", rotation=270, labelpad=20, fontsize=11)
    
    title = 'Moment 0 Map'
    if args.channel:
        title += f' (Chan: {args.channel[0]}-{args.channel[1]})'
    elif args.freq:
        title += f' ({args.freq[0]}-{args.freq[1]} MHz)'
    ax.set_title(title, pad=15, fontsize=14, fontweight='bold')

    print(f"正在保存图片: {output_png_file}")
    plt.savefig(output_png_file, dpi=200, bbox_inches='tight')
    print("完成。")
    
if __name__ == '__main__':
    run_clean_plot()
