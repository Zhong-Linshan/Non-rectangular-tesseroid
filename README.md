# Revisiting terrain representation in gravity forward modeling using non-rectangular tesseroids

This repository contains MATLAB code for the forward modelling of the
vertical gravitational attraction of a tesseroid (spherical prism).
Three families of integration rules are implemented and compared:

- **Gauss–Legendre quadrature (GLQ)**
- **Closed Newton–Cotes quadrature (CNCR)**
- **Open Newton–Cotes quadrature (ONCR)**

Both synthetic and real (SRTM) topographies are used to evaluate the
accuracy and computational performance of each method.  An analytical
solution obtained by integrating the kernel over longitude is used as
a reference.

The code is associated with the paper submitted to *Computers &
Geosciences*.  Every function and script is thoroughly commented and
ready for reuse.

## File inventory

### Core integration functions (accept all angles in **degrees**)

| File | Description |
|------|-------------|
| `GLQ_L.m` | Gauss–Legendre quadrature over a tesseroid.  Nodes and weights are provided by `nodes.m`. |
| `NCR_Close.m` | Closed Newton–Cotes quadrature (equally spaced nodes, endpoints included).  Uses an internal weight table for n = 2,…,7. |
| `NCR_Open.m` | Open Newton–Cotes quadrature (nodes at midpoints, endpoints excluded).  Internal weight table for n = 2,…,7. |
| `nodes.m` | Computes Gauss–Legendre nodes and weights on [-1, 1] via the Golub–Welsch algorithm.  This is the only external node generator; the Newton–Cotes routines use their own local `nodes` subfunctions. |

### Analytical reference function

| File | Description |
|------|-------------|
| `Vnorth_pole.m` | Vertical gravity of a tesseroid – analytical integration over longitude, leaving a 1D integral over co‑latitude.  **Expects angles in radians.**  Contains the private helper `Vvr`. |

### Scripts for synthetic topography experiments

| Script | Purpose |
|--------|---------|
| `script_synthetic_topo.m` | Generates a multi‑peak synthetic elevation model, applies progressive Gaussian filtering, and plots the surfaces (Figure 2). |
| `script_monte_carlo_comparison.m` | Monte‑Carlo accuracy test: 30 random realisations of the synthetic topography, compares six quadrature strategies (GLQ/CNCR/ONCR with variable or constant top) against the `Vnorth_pole` reference (Figure 3). |

### Scripts for real SRTM experiments

| Script | Purpose |
|--------|---------|
| `script_srtm_four_regions.m` | Loads 3″ SRTM tiles for Grand Canyon, Qinghai–Xizang Plateau, Mount Qomolangma, and Sichuan Basin.  Computes elevation statistics and (optionally) runs gravity forward modelling on 800×800 grids (Figures 4 & 5). |
| `script_srtm_convergence_speed.m` | Uses a single 1°×1° SRTM block (Qinghai–Xizang) to study the convergence and timing of GLQ with orders from 2 to 600, and also measures the fixed‑order Newton–Cotes rules. |

### Data files (not included – see data section)

- `srtm_15_05.tif`
- `srtm_57_06.tif`
- `srtm_54_07.tif`
- `srtm_58_06.tif`
- `srtm_57_07.tif`
- `srtm_58_07.tif`

## System requirements

- MATLAB R2016b or later (implicit array expansion is used).
- Image Processing Toolbox (for `imread` of GeoTIFF files; scripts that use SRTM data require it).
- Optional: `tight_subplot` from File Exchange (used in the Monte‑Carlo comparison script for figure layout).

The code was developed and tested under MATLAB R2020b.

## Installation

1. Clone or download this repository.
2. Ensure all `.m` files are in the same folder or on your MATLAB path.
3. If you plan to run the real‑data scripts, download the required SRTM
   tiles (3 arc‑second, GeoTIFF format) from a source such as
   [EarthExplorer](https://earthexplorer.usgs.gov/) or the SRTM
   downloader of your choice, and place them in the same folder as the
   scripts.

## Running the code

### Synthetic topography experiments

1. Open `script_synthetic_topo.m` and run it.  
   It creates four figures showing the effect of Gaussian filtering.
2. Open `script_monte_carlo_comparison.m` and run it.  
   This will take a few minutes and produce a plot of log‑relative errors
   for six quadrature configurations.

### Real SRTM experiments

1. Make sure the required `.tif` files are present.
2. Run `script_srtm_four_regions.m`.  
   By default, only the topography visualisation and the analytical
   reference computation are active; to run the full quadrature comparison,
   uncomment the relevant lines (see in‑script comments).
3. Run `script_srtm_convergence_speed.m`.  
   This performs a detailed timing and accuracy study for the GLQ method
   at many quadrature orders and records the results for the Newton–Cotes
   methods.

## Important notes on angle units

- `GLQ_L`, `NCR_Close`, `NCR_Open` expect **all angular inputs in degrees**
  and convert them to radians internally.
- `Vnorth_pole` expects **all angular inputs in radians**.  
  The scripts handle the conversions explicitly (`deg2rad`/`rad2deg`) where
  needed; please keep this distinction in mind when writing your own
  top‑level scripts.

## Core function signatures

```matlab
v = GLQ_L(x, x1, x2, y, y1, y2, z, z1, H, t, w, str)
v = NCR_Close(x, x1, x2, y, y1, y2, z, z1, H, n, str)
v = NCR_Open(x, x1, x2, y, y1, y2, z, z1, H, n, str)
v = Vnorth_pole(r, r1, r2, phi1, phi2, lam1, lam2)
[t, w] = nodes(n)
