# Introduction to High Performance Computing Wintger 20/21

Exercises for LEcture

## Usage 

### ex02

```bash
cd ex02
make doc             # builds latex document
make <part_name>     # builds src/<part_name>.cpp into <part_name>.out
make sync            # pushes src/ and scripts to hpc
make run/<part_name> # builds and enqueues into slurm
make pull            # pulls out/ from hpc
```
