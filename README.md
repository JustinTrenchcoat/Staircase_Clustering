# Staircase Clustering Detection Algorithm

Written by Mikhail Schee for:
Mikhail Schee, Erica Rosenblum, Jonathan M. Lilly, and Nicolas Grisouard (2023) "Unsupervised Clustering Identifies Thermohaline Staircases in the Canada Basin of the Arctic Ocean"

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
* Email: mikhail.schee@mail.utoronto.ca
* GitHub: https://github.com/scheemik/Staircase_Clustering_Detection_Algorithm

## Summary

This repository contains the code used by the above study to apply the Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN) clustering algorithm to data from Ice Tethered Profilers.

A detailed explanation of how the code was used to make each plot in the study can be found in the accompanying Jupyter notebook `Create_Figures.ipynb`.

## Acknowledgements

The Ice-Tethered Profiler data were collected and made available by the Ice-Tethered Profiler Program (Toole et al., 2011; Krishfield et al., 2008) based at the Woods Hole Oceanographic Institution https://www2.whoi.edu/site/itp/

This repository includes the [orthoregress code](https://gist.github.com/robintw/d94eb527c44966fbc8b9) from [Robin Wilson](https://blog.rtwilson.com/orthogonal-distance-regression-in-python/).

This research was supported in part by the National Science Foundation under Grant No. NSF PHY-1748958. 

We acknowledge fruitful discussions with Maike Sonnewald and Carine van der Boog.

M.S. and N.G. were supported by the Natural Sciences and Engineering Research Council of Canada (NSERC) [funding reference numbers RGPIN-2015-03684 and RGPIN-2022-04560]. J.M.L. was supported by grant number 2049521 from the Physical Oceanography program of the United States National Science Foundation.

