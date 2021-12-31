# ForRF
Time iterative deconvolution in fortran with Intel complier. The Intel MKL library is used for FFT

1. change directory of Intel complier in the `Makefile`
2. `make`
3. `./make_rf`

## Comparison with Seispy
<img src=https://user-images.githubusercontent.com/7437523/147818243-8d5bbfd8-0ae4-4d8c-a3fc-a0e6677d6701.png width=60% />

## Defect
The amplitude at peaks may not be exactly consistent with Seispy.
