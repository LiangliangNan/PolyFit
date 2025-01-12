# PolyFit Python Bindings

Python bindings for the PolyFit reconstruction method.

More information about PolyFit:
https://github.com/LiangliangNan/PolyFit

### Build PolyFit bindings

Make sure you have [Python](https://www.python.org/downloads/) installed.
Run CMake and then build. After building PolyFit, you
can find the Python bindings module `PyPolyFit` (`PyPolyFit.pyd` on Windows, `PyPolyFit.so` on macOS/Linux) in
`YOUR_BUILD_DIRECTORY/lib`.

### Use PolyFit bindings in Python code

After building the Python bindings, you should already be able to use the Python bindings module in your Python code.
See the [examples](./Examples).

**Note**: the Python version used for running the code should match the version used for building the bindings.
