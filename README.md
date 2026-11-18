# spatial-stdp

MATLAB source code for Figures 1 and 2 of **Formation and maintenance of phase
organization in spatially embedded spiking networks**, by Xiuying Zhou and Dong Yu.

## Scope / 内容

This repository contains MATLAB source files and this README only. The code
implements spatial spike-timing-dependent plasticity, delayed spiking-network
simulation, phase classification, and analysis and plotting for Figures 1 and 2.

本仓库仅含 MATLAB 源码与 README，不附带数据、绘图模板或生成结果。
运行现有入口前，需要另行准备下面列出的输入文件；仅下载本仓库不能直接重现图件。

## Requirements

- MATLAB R2024b is the tested baseline.
- Statistics and Machine Learning Toolbox is required for phase classification.
- Parallel Computing Toolbox is optional for a parallel full parameter scan;
  the scan can also run serially.

## Code layout

```text
run_figures12.m                 Main entry point
function/plot_figure1.m         Figure 1 plotting
grid/reproduce_grid.m          Archived voting, pilot and full parameter scan
grid/grid_simulate_vote.m      One simulation and state classification
grid/lib/                      LIF simulator, STDP and phase classification
reps/reproduce_reps.m           Figure 1(b-g) representative trajectories
fig2/reproduce_figure2.m        Paired S/N formation analysis and Figure 2
fig2/private/modularity_dir.m   Third-party directed modularity function
```

## Required inputs, supplied separately

The current entry points expect the following paths relative to the repository
root. These files are **not included** in this source-only repository.

| Path | Purpose |
| --- | --- |
| `function/figure1_palette.mat` | Figure 1 colors |
| `grid/data/scan_inputs.mat` | Fixed coordinates, model parameters and grid settings |
| `grid/data/trajectory_votes.csv` | Archived per-realization state votes |
| `grid/data/map_reference.csv` | Manuscript map counts used for validation |
| `reps/data/panel_b.mat` through `panel_g.mat` | Example parameters, seeds, spikes and validation hashes |
| `fig2/data/paired_records.mat` | Paired S/N inputs and archived spike/weight records |
| `fig2/data/layout.fig` | Figure 2 graphics layout |
| `fig2/data/reference_panels.mat` | Figure 2 panel validation arrays |

The fixed spatial input used in the study is a 100-neuron subset of FlyWire v783
coordinates, converted from nanometers to micrometers, from the archive associated
with Zhang et al., *Geometric Scaling Law in Real Neuronal Networks*, Physical
Review Letters 133, 138401 (2024), DOI 10.1103/PhysRevLett.133.138401.
Upstream coordinate repository: https://github.com/xinyacheung/geoscaling.
That upstream archive does not by itself provide all of the study-specific input
files listed above.

## Running, after the required inputs are available

Open the repository folder in MATLAB and use:

```matlab
run_figures12                 % Recompute from archived votes, spikes and weights
run_figures12('resimulate')   % Reintegrate examples and paired S/N trajectories
run_figures12('pilot')        % Reintegrate one boundary-audited realization
run_figures12('full')         % Reintegrate the full parameter grid and examples
```

All modes still require their corresponding input/configuration files above.
The default output location is `output/<mode>/`. A second argument selects a
different writable output folder. For a selected grid subset or worker count:

```matlab
addpath('grid');
reproduce_grid('full','output/grid_subset',1,[4 2414]);
```

The full grid is a long calculation. A historical estimate on the original
machine was 35-45 hours with 24 processes; runtime depends on hardware and the
number of tie-breaking realizations. The source-only packaging step did not run
a new simulation. The unchanged code was previously checked with the complete
input package; this source-only distribution is not a self-contained data release.

## Third-party code

`fig2/private/modularity_dir.m` is from the Brain Connectivity Toolbox, with its
original author header preserved. The pinned distribution is
https://github.com/brainlife/BCT/tree/b7f59ebcf0b4f4979ab3914e30c61eaf345b11d5.
Its distribution's license notice is included below in full, so a separate
license file is not needed in this source-only upload. The notice applies to
that third-party component; no new license for the authors' other code is
selected here.

```text
MIT License

Copyright (c) 2020 brainlife.io

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
