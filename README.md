# DEGREE-CORRECTED RICCI CURVATURE FOR NETWORKS

This is a repository for a paper titled *Degree-corrected Ricci Curvature for Networks*.

## Repository Structure
```
dcrc/
├── data/                # Raw real datasets
├── code/                
│   ├── algorithms/      # Folder of community detection algorithm files
│   ├── matlab/          # Folder for conducting community detection and evaluation
│   ├── prev_curv/       # Folder for BFC calculation codes
│   ├── simulation/      # Folder for simulated experiments in the paper
│   └── real_data/       # Folder for real datat applications
└── README.md            
```

## Package requirement
### Python
- Python >= 3.9
- [PyTorch](https://pytorch.org/get-started/locally/) >= 2.1.2
- [Networkx](https://networkx.org/documentation/stable/install.html)
- [matplotlib](https://matplotlib.org) >= 3.8.2
- [numpy](https://numpy.org) >= 1.26.2
- [pandas](https://pandas.pydata.org) >= 2.1.4
- [seaborn](https://seaborn.pydata.org) >= 0.13.0
- torchVision
- torchaudio
### R
-R 4.4.0
- ggplot2
- igraph
- R.matlab
- readr
- dplyr
- tidyr
- stringr
- broom
- rstatix
- ggpubr
- gridExtra
- showtext


  
## References
- Yudong Chen, Xiaodong Li, and Jiaming Xu. Convexified modularity maximization for degree-corrected stochastic block models. Annals of Statistics, 46(4):1573–1602, 2018. doi: 10.1214/
17-AOS1595.
- Zhuang Ma, Zongming Ma, and Hongsong Yuan. Universal latent space model fitting for large networks with edge covariates. Journal of Machine Learning Research, 21(4):1–67, 2020. URL
http://jmlr.org/papers/v21/17-470.html.
- Yuan Zhang, Elizaveta Levina, and Ji Zhu. Detecting overlapping communities in networks using spectral methods. SIAM Journal on Mathematics of Data Science, 2(2):265–283, 2020.
- Tai Qin and Karl Rohe. Regularized spectral clustering under the degree-corrected stochastic block model. Advances in neural information processing systems, 26, 2013.
- Jiashun Jin. Fast community detection by score. Annals of Statistics, 43(1):57–89, February 2015. doi: 10.1214/14-AOS1265.
- Jiashun Jin, Zheng Tracy Ke, and Shengming Luo. Improvements on score, especially for weak signals. Sankhya A, 84(1):127–162, March 2021. ISSN 0976-8378. doi: 10.1007/s13171-020-00240-1. URL http://dx.doi.org/10.1007/s13171-020-00240-1.

