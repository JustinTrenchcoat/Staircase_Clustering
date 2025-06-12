# Staircase Clustering Detection Algorithm

Written by Mikhail Schee for:
Mikhail Schee, Erica Rosenblum, Jonathan M. Lilly, and Nicolas Grisouard (2024) "Unsupervised Clustering Identifies Thermohaline Staircases in the Canada Basin of the Arctic Ocean," Environmental Data Science, 3:e13, 1-19, doi:[10.1017/eds.2024.13](https://doi.org/10.1017/eds.2024.13)

## License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions in source code must retain the accompanying copyright notice, this list of conditions, and the following disclaimer.
    2. Redistributions in binary form must reproduce the accompanying copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3. Names of the copyright holders must not be used to endorse or promote products derived from this software without prior written permission from the copyright holders.
    4. If any files are modified, you must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Contact

* Corresponding Author: Mikhail Schee (he/him)
* Email: [mikhail.schee@alumni.utoronto.ca](mailto:mikhail.schee@alunmi.utoronto.ca)
* GitHub: https://github.com/scheemik/Staircase_Clustering_Detection_Algorithm

## Summary

This repository contains the code used by the above study to apply the Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN) clustering algorithm to data from Ice Tethered Profilers.

A detailed explanation of how the code was used to make each plot in the study can be found in the accompanying Jupyter notebook `Create_Figures.ipynb`.

## Acknowledgements

We acknowledge fruitful discussions with Maike Sonnewald and Carine van der Boog. M.G.S. acknow- ledges fruitful discussions that occurred during the Kavli Institute of Theoretical Physics program on Layering in Atmospheres, Oceans and Plasmas (supported by the National Science Foundation under grant no. NSF PHY-1748958). E.R. is grateful to the researchers, staff, and students of the Centre for Earth Observation Science for the support received during the preparation of this manuscript. J. M. L. thanks the Department of Physics at the University of Toronto, where he was a visiting scientist during the time this work was carried out.

This repository includes the [orthoregress code](https://gist.github.com/robintw/d94eb527c44966fbc8b9) from [Robin Wilson](https://blog.rtwilson.com/orthogonal-distance-regression-in-python/).

## Data availability statement

The Ice-Tethered Profiler data were collected and made available by the Ice-Tethered Profiler Program (Toole et al., 2011; Krishfield et al., 2008) based at the Woods Hole Oceanographic Institution https://www2.whoi.edu/site/itp/

## Funding statement

M.G.S. and N.G. were supported by the Natural Sciences and Engineering Research Council of Canada (NSERC) (funding reference numbers RGPIN-2015-03684 and RGPIN-2022-04560). J.M.L. was supported by grant number 2049521 from the Physical Oceanography program of the United States National Science Foundation. ER was supported by the NSERC Postdoctoral Fellowship award, the NSERC Canada-150 Chair (Award G00321321), and the NSERC Discovery Grant (funding reference number RGPIN-2024-05545).


## Development Log:

Below is a comprehensive development log with dates and brief comments.
2025/06/12:
- Tried to follow `Create_Figures.ipynb`'s instruction but unable to produce figure in step 2. Investigation turned out that the dataframe was empty.
- isolated Figure 2 section from `Create_Figures.ipynb` into `experiment.ipynb`. Reproduced figure 2 successfully.
- ITP3 clustering takes a long time and takes up most of computing power, might need to reconsider its performance when put to actual use
- Will fully assess its performance on ITP2