ADC Analysis Project

Description:
This project simulates and analyzes ADC (Analog-to-Digital Converter) data for a CMS-like detector.
The goal is to study pedestal, intrinsic noise, and common-mode noise effects, and to perform per-cell statistical analysis.

Datasets:

ADC1: Pedestal with mean = 120

ADC2: Pedestal + intrinsic noise (σ = 4)

ADC3: Pedestal + intrinsic noise + common-mode noise (σ = 2)

Pedestal subtraction and common-mode subtraction are applied to study corrected ADC distributions.

Analysis:

Compute per-cell mean and standard deviation

Generate histogram and line plots for individual cells

Apply common-mode subtraction to remove event-level correlated noise

Create observational summary tables with min, max, mean, and standard deviation for each cell


code is available on colab for viewers ...
https://colab.research.google.com/drive/1uXZkUw42riwqpAY9CBEYA2vfz8m1bguS?usp=sharing
