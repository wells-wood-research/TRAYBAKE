# External dependency

`setup_environment.sh` checks out [QupKake](https://github.com/Shualdon/QupKake)
at a pinned upstream commit in the ignored `dependencies/QupKake/` directory,
applies `patches/qupkake-traybake-compat.patch`, and installs it into the
project-local Conda environment. Keep TRAYBAKE-specific modifications in that
patch instead of committing a modified upstream checkout. QupKake's generated
files are kept with each run under `test-results/` or the output directory
supplied by the user.
