# Third-party notices

This analysis code is GPL-3.0-or-later; the full license is in LICENSE.
The archived R source excerpts and report screenshots derive from susieR.
No association statistics or individual genotypes are redistributed.

### susieR

Gaussian SuSiE and SuSiE-RSS code is adapted from
[susieR 0.16.6](https://github.com/stephenslab/susieR/tree/8e56a8e038e989856d106d9ca5175cc664fea9d2),
commit `8e56a8e038e989856d106d9ca5175cc664fea9d2`, under BSD-3-Clause:

Copyright (c) 2017-2022, Gao Wang, Peter Carbonetto, Yuxin Zou, Kaiqian Zhang,
Matthew Stephens.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of Gao Wang, Peter Carbonetto, Yuxin Zou, Kaiqian
   Zhang, Matthew Stephens, nor the names of its contributors may be used
   to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


## R optimizer attribution

The archived native timing patch is part of prusie's GPL-3.0-or-later source.
The prusie dependency adapts R's `optimize.c` under GPL-2.0-or-later:
Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka;
Copyright (C) 2003–2004 The R Foundation;
Copyright (C) 1998–2023 The R Core Team.
The optimizer implementation and compiled package wheel are not redistributed
in this analysis repository. Their notices remain in the package distribution.

## Data providers and plotting libraries

FinnGen, GTEx and IGSR data access and attribution are described in README.md.
NumPy, SciPy, Matplotlib, Pillow, PyYAML, pandas, PyArrow, pysam and pgenlib
are installed dependencies, not vendored source or bundled binaries.
The GPL license text and the susieR BSD notice are retained in full.
