# RARE
**RARE** is a multivariable Mendelian randomization method, and is the short name of 'MVMR incorporating **R**are variants **A**ccounting for multiple **R**isk factors and shared horizontal pl**E**iotopy'.

RARE is the only method (2024-07-31) that accounts for the impact of rare variants in causal inference while simultaneously considers UHP and CHP.

image/Github_RARE.jpg

# Installation
Install this tool by use of the 'devtools' package. Note that RARE partly depends on the C++ languange, thus you should appropriately set 
Rtools and X code for Windows, Mac OS/X, and Linux, respectively.
```
install.packages( 'devtools' )  
library( devtools )  
install_github( 'Hide-in-lab/RARE@main', force  = T )
```
# References
Yu Cheng<sup>1</sup>, Xinjia Ruan<sup>1</sup>, Tiantian Liu<sup>#</sup>, Fangrong Yan<sup> #</sup>. **Accounting for the impact of rare variants on causal inference with RARE: A novel multivariable Mendelian randomization method**

# Development
This package is developed and maitained by Yu Cheng (yucheng.cpu@foxmail.com) and Xinjia RUan (ruan.cpu@foxmail.com). Cite the code: ___________.
