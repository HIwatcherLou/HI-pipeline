# FAST Data Processing Toolkit

Welcome to the **FAST Data Processing Toolkit**! This repository provides a comprehensive suite of Python scripts and Jupyter Notebooks designed to streamline and enhance the data reduction and analysis workflow for the Five-hundred-meter Aperture Spherical radio Telescope (FAST).

Built to complement the `hifast` pipeline and the `SoFiA` (Source Finding Application) software, these tools cover the entire lifecycle of radio astronomy data processing—from raw FITS data grouping, observation verification, and chunk merging, to final scientific visualization and automated source cross-matching.

## Workflow & Tool Index

This toolkit is structured around a standard FAST data reduction and analysis workflow:

### 1. Data Preprocessing & Formatting
* **`transformation.py`**: Scans raw chunked FITS files, intelligently groups them by Beam ID, sorts them chronologically, and converts them into `hifast`-compatible HDF5 formats.

### 2. Quality Control & Observation Verification
* **`RA-DEC_total.py`**: Extracts telescope pointing coordinates, filters out unstable slewing periods, and visualizes the actual RA/DEC drift tracking paths to ensure observation accuracy.
* **`RMS_analysis.py`**: Iterates through FITS datacubes to calculate per-channel RMS noise, automatically flagging anomalous channels (RFI) using a robust statistical threshold.
* **`path.ipynb`**: A lightweight utility notebook for quickly inspecting FITS header WCS information (e.g., pixel size in degrees/arcminutes) and plotting RA/DEC drift paths directly from intermediate HDF5 files.

### 3. Data Merging
* **`merge.py`**: Safely concatenates time-chunked HDF5 segments back into continuous, full-beam files. It features dynamic memory resizing and strict preservation of `hifast` soft-links (crucial for Carta compatibility).

### 4. Scientific Visualization & Extraction
* **`my_analysis.py`**: A smart batch-plotter that uses a 1D convolution algorithm to dynamically locate the "cleanest" continuous time slices, extracting and plotting high-quality mean spectra while bypassing RFI.
* **`moment0.py`**: Computes 2D Integrated Intensity (Moment 0) maps from 3D FITS datacubes. Fully optimized for headless servers, supporting flexible spectral slicing and advanced visual contrast stretches (Log, Sqrt, Square, Exp).

### 5. Source Finding Analysis & Cross-Matching (SoFiA Integration)
* **`source.py`**: Reads SoFiA catalog outputs and generates a spatial distribution scatter plot (RA vs. Dec) of detected HI sources. It dynamically scales point sizes by Signal-to-Noise Ratio (SNR) and colors them by frequency, automatically annotating the top candidates.
* **`check_source.ipynb`**: An interactive Jupyter Notebook leveraging `ipywidgets` to quickly scan through 1D spectra (`_spec.txt`) of SoFiA-generated cubelets. Perfect for visually verifying the integrated flux of hundreds of sources via a seamless slider interface.
* **`search.ipynb`**: An automated cross-matching tool. It calculates the HI redshift ($z_{HI}$) for detected sources and queries the NASA/IPAC Extragalactic Database (NED) via `astroquery`. It identifies optical/galaxy matches based on spatial radius and redshift tolerance, outputting a clear Pandas DataFrame summary.

## Installation & Dependencies

Ensure you have the `hifast` pipeline installed in your environment. Clone this repository and install the required packages:

git clone [https://github.com/yourusername/FAST-Data-Toolkit.git](https://github.com/yourusername/FAST-Data-Toolkit.git)
cd FAST-Data-Toolkit
pip install numpy scipy pandas h5py astropy matplotlib ipywidgets notebook astroquery

# FAST 射电天文数据处理工具箱 (Project Introduction)

欢迎使用 FAST 射电天文数据处理工具箱！本仓库提供了一整套基于 Python 的自动化脚本和 Jupyter Notebook，旨在简化并增强 FAST (500米口径球面射电望远镜) 观测数据的后处理与分析工作流。

这套工具链完美契合 `hifast` 处理管线以及 `SoFiA` (自动寻源软件)，覆盖了射电天文数据处理的全生命周期——从原始 FITS 碎片的智能分组与格式转换，到望远镜轨迹核查、HDF5 分块无损合并，再到科学级数据可视化以及最终的寻源结果自动化交叉匹配。

## 工作流与核心工具索引

本工具箱按照标准的 FAST 数据处理与寻源分析逻辑划分为以下五个模块：

### 1. 数据预处理与格式转换
* **`transformation.py`**: 自动扫描原始 FITS 碎片文件，按波束 (Beam) 智能分组、按时间顺序排列，并转化为完全兼容 `hifast` 格式的 HDF5 文件。

### 2. 数据质量控制与观测核查
* **`RA-DEC_total.py`**: 提取望远镜指向数据，智能剔除转场或不稳定的扫描片段，将实际的 RA/DEC 漂移轨迹合并绘制在同一张高对比度图表中。
* **`RMS_analysis.py`**: 遍历 FITS 数据立方体计算逐通道的均方根 (RMS) 噪声，并利用统计阈值自动检测和标记受射频干扰 (RFI) 污染的异常通道。
* **`path.ipynb`**: 轻量级实用工具本。用于快速检查 FITS 头文件中的 WCS 像素尺度（度/角分），并支持直接从 HDF5 过程文件中提取和绘制 RA/DEC 实际指向轨迹。

### 3. 分块数据合并
* **`merge.py`**: 将分块处理后的 HDF5 切片安全、无损地拼接回完整的单波束文件。支持动态扩容写入，并完美保留软链接（这对 Carta 等高级可视化软件至关重要）。

### 4. 科学级可视化与特征提取
* **`my_analysis.py`**: 智能频谱提取与批量绘图工具。采用一维卷积算法动态寻找“最干净”的连续时间切片，完美避开 RFI 密集区，提取高质量平均频谱并应用高斯平滑。
* **`moment0.py`**: 从 3D FITS 数据立方体生成 2D 积分强度图 (Moment 0)。专为无头服务器环境优化，支持灵活的光谱/通道切片及多种非线性视觉拉伸算法（对数、开方、平方、指数）。

### 5. 自动寻源分析与交叉匹配 (SoFiA 结果集成)
* **`source.py`**: 读取 SoFiA 寻源结果目录，绘制空间分布散点图 (RA vs. Dec)。自动根据信噪比 (SNR) 调整点的大小，以颜色映射观测频率，并自动标注出信噪比最高的前 10 个候选源。
* **`check_source.ipynb`**: 基于 `ipywidgets` 的交互式 Jupyter Notebook，用于快速浏览 SoFiA 生成的每一个源的 1D 频谱 (`_spec.txt`)。通过拖动滑动条即可无缝切换查看数百个源的积分通量，极大提升人工核查效率。
* **`search.ipynb`**: 自动化数据库交叉匹配工具。根据观测频率计算中性氢红移 ($z_{HI}$)，利用 `astroquery` 模块自动批量连接 NED (NASA/IPAC Extragalactic Database) 数据库。在指定的搜索半径和红移误差范围内寻找对应天体，并输出直观的匹配结果汇总表。

## 环境配置

> **⚠️ 注意：** > 请确保你的运行环境中已经成功安装了 `hifast` 管线以及其他诸如numpy scipy pandas h5py astropy matplotlib等库。
