# FAST Data Processing Toolkit

Welcome to the **FAST Data Processing Toolkit**! This repository provides a comprehensive suite of Python scripts designed to streamline and enhance the data reduction workflow for the Five-hundred-meter Aperture Spherical radio Telescope (FAST).

Built to complement the `hifast` pipeline, these tools cover the entire lifecycle of radio astronomy data processing—from raw FITS data grouping and preprocessing, to telescope tracking verification, chunk merging, and final scientific visualization (Spectra and Moment 0 maps).

## Workflow & Tool Index

This toolkit is structured around a standard FAST data reduction workflow:

### 1. Data Preprocessing & Formatting
* **`transformation.py`**: Scans raw chunked FITS files, intelligently groups them by Beam ID, sorts them chronologically, and converts them into `hifast`-compatible HDF5 formats.

### 2. Quality Control & Observation Verification
* **`RA-DEC_total.py`**: Extracts telescope pointing coordinates, filters out unstable slewing periods, and visualizes the actual RA/DEC drift tracking paths to ensure observation accuracy.
* **`RMS_analysis.py`**: Iterates through FITS datacubes to calculate per-channel RMS noise, automatically flagging anomalous channels (RFI) using a robust statistical threshold.

### 3. Data Merging
* **`merge.py`**: Safely concatenates time-chunked HDF5 segments back into continuous, full-beam files. It features dynamic memory resizing and strict preservation of `hifast` soft-links (crucial for Carta compatibility).

### 4. Scientific Visualization & Extraction
* **`my_analysis.py`**: A smart batch-plotter that uses a 1D convolution algorithm to dynamically locate the "cleanest" continuous time slices, extracting and plotting high-quality mean spectra while bypassing RFI.
* **`moment0.py`**: Computes 2D Integrated Intensity (Moment 0) maps from 3D FITS datacubes. Fully optimized for headless servers, supporting flexible spectral slicing and advanced visual contrast stretches (Log, Sqrt, Square, Exp).

## Installation & Dependencies

Ensure you have the `hifast` pipeline installed in your environment.
