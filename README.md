# cFIT

R package for cFIT.

cFIT (common Factor Space Integration & Transfer) is a tool for data integration and transfer for scRNAseq data. It is applied to data from multiple labs and experimental conditions, technologies and even species. The proposed method models the shared information between various data sets by a common factor space, while allowing for unique distortions and shift per batch. The model parameters are learned under and iterative non-negative matrix factorization (NMF) framework and then used for synchronized intefration from across-domain assays. In addtion, the model enables transferring via low-rank matrix from more informative data to allow for precise identification of data of lower quality.


## Citation

A preprint of the manuscript is available [here](https://www.biorxiv.org/content/10.1101/2020.08.31.276345v1)
```
@article{peng2021integration,
  title={Integration and transfer learning of single-cell transcriptomes via cFIT},
  author={Peng, Minshi and Li, Yue and Wamsley, Brie and Wei, Yuting and Roeder, Kathryn},
  journal={Proceedings of the National Academy of Sciences},
  volume={118},
  number={10},
  year={2021},
  publisher={National Acad Sciences}
}
```

## Installation
This package can be installed through `devtools` in R:
```{r}
library("devtools")
devtools::install_github("pengminshi/cFIT")
```

## Examples
Please follow the [vignette](https://htmlpreview.github.io/?https://github.com/pengminshi/cFIT/blob/master/vignettes/vignette.html) for examples of using this R package on simulated and real data sets.
