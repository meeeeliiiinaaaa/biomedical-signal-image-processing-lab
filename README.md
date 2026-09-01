# Biomedical Signal and Image Processing Lab

MATLAB implementations for the **Biomedical Signal and Image Processing Laboratory** course at Sharif University of Technology. This repository contains code, analysis, and results for nine lab assignments covering biomedical signal acquisition, filtering, source separation, pattern detection, and medical image processing.

## Course Overview

The labs progress from foundational signal analysis to advanced source separation and image segmentation techniques, applied to real and simulated biomedical data (EEG, ECG, EMG, EOG, and MRI/CT).

## Labs

| # | Title | Topics |
|---|---|---|
| 1 | **Multi-domain Display of EEG and ECG Signals** | Sampling, DFT, STFT/spectrogram, EEG channel visualization, ECG morphology (P-QRS-T), EMG (healthy vs. neuropathy vs. myopathy), EOG |
| 2 | **EEG Noise and Artifact Removal** | Independent Component Analysis (ICA) for muscle artifact removal, SNR-controlled noise injection, relative RMSE evaluation, artifact rejection on real epileptic EEG |
| 3 | **Maternal-Fetal ECG Separation** | Blind Source Separation (BSS): Singular Value Decomposition (SVD/PCA) vs. ICA, eigenspectrum analysis, correlation-based evaluation of recovered fetal ECG |
| 4 | **EEG Patterns and Detection Methods** | Event-Related Potentials (P300, synchronous averaging), Steady-State Visual Evoked Potentials (SSVEP), Event-Related Synchronization/Desynchronization (ERS/ERD) across frequency bands (delta, theta, alpha, beta) and mental tasks |
| 5 | **ECG Patterns and Cardiac Arrhythmia Detection** | ECG morphology and arrhythmias (AFL, AFib, VFL, VFib), bandpass filtering for noise/baseline wander removal, frequency- and morphology-based feature extraction, automatic ventricular fibrillation detection (accuracy/sensitivity/specificity) |
| 6 | **Brain Source Localization** | Three-shell spherical head model, forward lead-field modeling, dipole simulation, EEG potential topography, inverse problem solving via MNE, parametric dipole localization (genetic algorithm / simulated annealing) |
| 7 | **Introduction to Image Processing in MATLAB** | Image representation and coordinate systems, 2D Fourier Transform, convolution theorem, frequency-domain filtering, image rotation/shift properties, edge detection (Sobel, Canny), gradient-based operators |
| 8 | **Medical Image Denoising** | Noise models (Gaussian, Rayleigh, Gamma, exponential, uniform, salt & pepper, Rician, speckle), median/mean/Gaussian smoothing filters, image deblurring and deconvolution, Gradient Descent-based restoration, anisotropic diffusion filtering |
| 9 | **Medical Image Segmentation** | Thresholding, Region Growing, K-means and Fuzzy C-Means (FCM) clustering for multi-modal MRI (T1, T2, PD) segmentation, comparison of hard vs. soft clustering assignments |


## Requirements

- MATLAB (R2020a or later recommended)
- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox (for `fcm`, clustering functions)

Some labs use provided helper functions and datasets supplied with the original course materials (e.g., `disp_eeg.m`, `ica.m`, `plotEEG.m`, `shell3_ForwardModel.m`, `dv3plot.m`, `ch3plot.m`). These are included alongside each lab's code where applicable.


## Datasets

Signals used across the labs include:
- EEG recordings from epileptic patients (ictal/interictal), simulated EEG with added noise, and multi-class mental-task EEG
- ECG recordings from the MIT-BIH Arrhythmia and Malignant Ventricular Arrhythmia databases
- Simulated maternal/fetal ECG and noise components
- EMG recordings (healthy, neuropathy, myopathy) from PhysioNet's EMG database
- T1/T2/PD-weighted MRI and CT images for segmentation and denoising experiments

Dataset files are not redistributed in this repository; see each lab's instructions for the original source.

## Contributors

- Reza Sanaeinejad
