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


### (Optional) Installation

If you want to avoid specifying the path of the generated bindings in your Python code, you can create a wheel (`.whl`)
file for the bindings and install it globally. This will make PoliFit work like other Python packages such as [Numpy](https://numpy.org/).

Below are the detailed steps for creating the wheel file and installation.

- Step 1: Install the build Tool (only if not installed)

    - `pip install build`

- Step 2: Build the Wheel

    - `cd YOUR_BUILD_DIRECTORY/lib/python`

    - `python -m build`

- Step 3: Install the Wheel

    - `pip install dist/polyfit-2.6.0-py3-none-any.whl` # Update the wheel file name accordingly

  You can verify the installation by running the following command:

    - `pip show polyfit`
