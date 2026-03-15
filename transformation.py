import numpy as np
import h5py
from astropy.io import fits
import os
import glob
import re
import argparse
import sys

def main():
    # ================= 1. 参数解析 =================
    parser = argparse.ArgumentParser(description='Auto Match & Merge Tool for FAST RAW Data')

    # 改动点：这里接收的是“文件夹路径”，而不是文件通配符
    parser.add_argument('input_dir', type=str, 
                        help='Input directory containing the _T.fits files')
    
    parser.add_argument('--outdir', type=str, required=True,
                        help='Directory to store output HDF5 files')
    
    parser.add_argument('--start', type=int, default=1,
                        help='Start chunk number (inclusive)')
    
    parser.add_argument('--stop', type=int, default=9999,
                        help='Stop chunk number (inclusive)')

    parser.add_argument('--frange', type=float, nargs=2, default=None,
                        help='Frequency range to keep (min max), e.g., --frange 1300 1450')

    args = parser.parse_args()

    # ================= 2. 扫描目录 =================
    if not os.path.isdir(args.input_dir):
        print(f"Error: 输入路径不是文件夹 -> {args.input_dir}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # 搜索文件夹下所有以 .fits 结尾的文件 (即使不是 _T 结尾也可以尝试匹配)
    search_pattern = os.path.join(args.input_dir, "*.fits")
    all_files = glob.glob(search_pattern)
    
    if len(all_files) == 0:
        print(f"Error: 在目录 {args.input_dir} 中没有找到 .fits 文件")
        sys.exit(1)

    print(f"扫描到 {len(all_files)} 个文件，正在进行 Beam 分组...")

    # ================= 3. 智能分组 (按 Beam) =================
    beam_groups = {}
    
    # 正则表达式：匹配 ...-M(数字)_W_(数字)_T.fits
    # Group 1: Scan Name (e.g. Dec-0011_...)
    # Group 2: Beam ID (e.g. 01)
    # Group 3: Chunk ID (e.g. 0001)
    # 兼容 _T.fits 和普通 .fits
    filename_pattern = re.compile(r'(.*)-M(\d+)_W_(\d+)(_T)?\.fits$')

    for file_path in all_files:
        filename = os.path.basename(file_path)
        match = filename_pattern.search(filename)
        
        if match:
            scan_name = match.group(1)
            beam_id = match.group(2)      # 01, 19...
            chunk_id = int(match.group(3)) # 1, 2...
            
            # 筛选：只处理在 start 和 stop 范围内的文件
            if args.start <= chunk_id <= args.stop:
                # 使用 Scan + Beam 作为唯一标识键
                key = f"M{beam_id}" 
                
                if key not in beam_groups:
                    beam_groups[key] = {
                        'scan': scan_name,
                        'beam': beam_id,
                        'files': [] 
                    }
                # 存入 (chunk_id, file_path) 以便后续排序
                beam_groups[key]['files'].append((chunk_id, file_path))

    # ================= 4. 循环处理每个 Beam =================
    sorted_beams = sorted(beam_groups.keys()) # 让处理顺序看起来整齐 (M01, M02...)
    
    if not sorted_beams:
        print(f"Warning: 在指定的 Chunk 范围 ({args.start}-{args.stop}) 内没有找到符合命名规则的文件。")
        sys.exit(0)

    print(f"即将处理以下 Beam: {sorted_beams}")
    print("-" * 60)

    for beam_key in sorted_beams:
        info = beam_groups[beam_key]
        scan = info['scan']
        beam = info['beam']
        
        # 核心步骤：按 Chunk ID 排序 (防止 1, 10, 2 这种乱序)
        sorted_files = sorted(info['files'], key=lambda x: x[0])
        
        # 提取排序后的路径和ID
        flist = [x[1] for x in sorted_files]
        chunk_ids = [x[0] for x in sorted_files]
        
        real_start = chunk_ids[0]
        real_end = chunk_ids[-1]

        print(f"正在处理 Beam M{beam} | Chunk: {real_start} -> {real_end} | 共 {len(flist)} 个文件")

        # --- 数据读取逻辑 (复用你提供的逻辑) ---
        freq_obs = []
        mjd_obs = []
        data_obs = []

        try:
            for idx, fpath in enumerate(flist):
                # print(f"  读取: {os.path.basename(fpath)}") # 可以注释掉减少刷屏
                with fits.open(fpath) as hdu:
                    # 假设结构: 0=Freq, 1=MJD, 2=Data
                    freq = hdu[0].data
                    mjd = hdu[1].data
                    data = hdu[2].data
                    
                    if len(flist) == 1:
                        freq_obs = freq
                        mjd_obs = mjd
                        data_obs = data
                    else:
                        freq_obs.append(freq)
                        mjd_obs.append(mjd)
                        data_obs.append(data)

            # 合并
            if len(flist) > 1:
                mjd_obs = np.concatenate(mjd_obs, axis=0)
                data_obs = np.concatenate(data_obs, axis=0)
                if isinstance(freq_obs, list):
                    freq_obs = np.array(freq_obs[0])

            # --- 频率切除 ---
            if args.frange is not None:
                f_min, f_max = args.frange
                freq_mask = (freq_obs >= f_min) & (freq_obs <= f_max)
                if np.sum(freq_mask) == 0:
                    print(f"  [跳过] Beam M{beam}: 频率范围内无数据")
                    continue
                freq_obs = freq_obs[freq_mask]
                data_obs = data_obs[:, freq_mask, :]

            # --- 写入 HDF5 ---
            # 文件名自动带上 chunk 范围
            out_filename = f"{scan}-M{beam}_W_{real_start:04d}-{real_end:04d}_specs_T.hdf5"
            output_path = os.path.join(args.outdir, out_filename)
            
            # 数据转置 (Time, Freq, Pol) -> (Pol, Time, Freq)
            # 根据你的 transformation.py 逻辑
            data_hdf5 = np.transpose(data_obs, (2, 0, 1)) 
            samp_n = len(mjd_obs)
            freq_n = len(freq_obs)

            with h5py.File(output_path, 'w') as f:
                f.create_group('Header')
                g_s = f.create_group('S')
                
                # 写入核心数据
                g_s.create_dataset('Ta', data=data_hdf5)
                g_s.create_dataset('freq', data=freq_obs)
                g_s.create_dataset('mjd', data=mjd_obs)
                
                # 写入占位符 (兼容 hifast 格式)
                g_s.create_dataset('Tcal', data=np.zeros((1, freq_n, 2)))
                g_s.create_dataset('inds_ton', data=np.zeros((6, 2)))
                g_s.create_dataset('is_aband_whole', data=np.zeros((6, 2)))
                g_s.create_dataset('is_delay', data=np.zeros(samp_n))
                g_s.create_dataset('is_on', data=np.full(samp_n, False, dtype=bool))
                g_s.create_dataset('next_to_cal', data=np.zeros(samp_n))
                g_s.create_dataset('pcals_amp_diff_interp_values', data=np.zeros((samp_n, 2)))
                g_s.create_dataset('pcals_merged', data=np.zeros((1, freq_n, 2)))
                g_s.create_dataset('pcals_merged_s', data=np.zeros((1, freq_n, 2)))

                f.create_group('Waterfall').create_dataset('DATA', data=data_hdf5)

            print(f"  -> 生成文件: {out_filename}")

        except Exception as e:
            print(f"  Error processing Beam M{beam}: {e}")
            import traceback
            traceback.print_exc()

    print("-" * 60)
    print("所有任务完成。")

if __name__ == "__main__":
    main()
