# PolyFit Python Bindings

Python bindings for the PolyFit reconstruction method.

More information about PolyFit:
https://github.com/LiangliangNan/PolyFit

### Build PolyFit bindings

Before building the Python bindings, ensure that [Python](https://www.python.org/downloads/) is installed on your system.
Then, run CMake to configure the project and build PolyFit (see detailed instructions [here](../../ReadMe.md)). 
After building PolyFit, the Python bindings module, `PyPolyFit`, will be located in your build directory at:
- Windows: `YOUR_BUILD_DIRECTORY/lib/PyPolyFit.pyd`
- macOS/Linux: `YOUR_BUILD_DIRECTORY/lib/PyPolyFit.so`

### Use PolyFit bindings in Python code

After successfully building the Python bindings, you can directly import the PyPolyFit module into your Python code.
To get started, refer to the [examples](./Examples) directory for code snippets demonstrating typical usage.

**Note**: the Python version used for running your code should match the version used for building the bindings.


### (Optional) Installing PolyFit Python bindings

If you prefer not to build PolyFit from source on each machine, you can create a wheel installer for easier 
distribution and installation. This will allow you to install PolyFit like any other Python package (e.g., [Numpy](https://numpy.org/)).

Follow these steps to build and install the wheel file:

- Step 1: Install the build tool (if not already installed)

    - `pip install build`

- Step 2: Navigate to the bindings directory and build the wheel

    - `cd YOUR_BUILD_DIRECTORY/lib/python`

    - `python -m build`

- Step 3: Install the wheel

    - `pip install dist/polyfit-1.6.0-py3-none-any.whl` # Update the wheel file name accordingly

  You can verify the installation by running the following command:

    - `pip show polyfit`
