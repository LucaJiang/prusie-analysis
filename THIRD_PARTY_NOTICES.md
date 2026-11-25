# Third-party notices

The analysis code in this repository is licensed under GPL-3.0-or-later; see
[LICENSE](LICENSE) for the full license text.

### susieR

The archived R source files under `02_benchmark/archive/` are copies of
[susieR 0.16.6](https://github.com/stephenslab/susieR/tree/8e56a8e038e989856d106d9ca5175cc664fea9d2)
at commit `8e56a8e038e989856d106d9ca5175cc664fea9d2`, distributed under the
BSD-3-Clause license:

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

The archived native timing patch (`02_benchmark/archive/native_timing.patch`)
applies to `prusie` (GPL-3.0-or-later), which adapts Brent's univariate
optimization routine from R's `optimize.c` under GPL-2.0-or-later:
Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka;
Copyright (C) 2003–2004 The R Foundation;
Copyright (C) 1998–2023 The R Core Team.

## Data providers and external packages

Data access terms and attributions for FinnGen, GTEx, and IGSR 1000 Genomes are
documented in [README.md](README.md). Python and system dependencies (NumPy,
SciPy, Matplotlib, Pillow, PyYAML, pandas, PyArrow, pysam, and pgenlib) are
installed via standard package managers as listed in `environment/`.
