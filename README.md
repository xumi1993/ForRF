# ForRF
Time iterative deconvolution in Fortran with CUDA acceleration.

## Requirements
- Fortran Compiler (gfortran/ifort)
- CUDA Toolkit (for GPU support)
- CMake >= 3.10

## Build
```bash
mkdir build && cd build
cmake .. && make
```

## Usage
Run the test program with an optional argument to select the device (default is GPU).

```bash
cd src
# Run on GPU (default)
./forrf_test gpu

# Run on CPU
./forrf_test cpu
```

## Comparison with Seispy
<img src=https://user-images.githubusercontent.com/7437523/147818243-8d5bbfd8-0ae4-4d8c-a3fc-a0e6677d6701.png width=60% />

## Comparison between CPU and GPU
<img src="https://github.com/user-attachments/assets/d6a90119-cd19-4f6d-a20a-b2ac8b4ef6f8" />

Performance comparison on a test dataset:
- CPU: ~3.07s
- GPU: ~0.23s (~13x speedup)

## Defect
The amplitude at peaks may not be exactly consistent with Seispy.
