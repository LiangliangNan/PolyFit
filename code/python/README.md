## PolyFit Python Bindings

Python bindings for the PolyFit reconstruction method.

For more information about PolyFit, visit:
https://github.com/LiangliangNan/PolyFit

### Obtain/Build PolyFit bindings

Prebuilt Python bindings are available for download on the [Releases](https://github.com/LiangliangNan/PolyFit/releases) page.
Simply Download the appropriate wheel file for your platform and Python version, and install it using `pip`:
```bash
pip install polyfit-1.6.0-py3-none-any.whl # <-- Use your actual wheel file name
```

#### Build from source

If prebuilt bindings are not available for your platform or specific Python version, you can build the bindings from the
source code. Here’s how:
- Ensure you have [Python](https://www.python.org/downloads/) installed on your system.
- Run CMake to configure the project and build PolyFit (see detailed instructions [here](../../ReadMe.md)). 
- After building, the Python bindings module will be located in `YOUR_BUILD_DIRECTORY/lib/python/polyfit`.

#### Use PolyFit bindings in Python code

After successfully building the Python bindings, you can directly import the PyPolyFit module into your Python code.
To get started, refer to the [examples](./Examples) directory for code snippets demonstrating typical usage.

**Note**: the Python version used for running your code should match the version used for building the bindings.


### (Optional) Create distributable installer for PolyFit Python bindings

If you prefer not to build PolyFit from source on each machine, you can create a wheel (`.whl`) installer for easier 
distribution and installation. This will allow you to install PolyFit like any other Python package (e.g., [Numpy](https://numpy.org/)).

Follow these steps to create and install the wheel file (after you have successfully built the bindings following the steps above):

- Step 1: Install the `build` tool (if not already installed)

    - `pip install build`

- Step 2: Navigate to the bindings directory and build the wheel

    - `cd YOUR_BUILD_DIRECTORY/lib/python`

    - `python -m build`

You can then find the generated wheel file in the `dist` directory. Users can use this wheel file to install PolyFit Python bindings without needing to build from source:

`pip install dist/polyfit-1.6.0-py3-none-any.whl` # <-- Use your actual wheel file name

You can verify the installation by running the following command:

`pip show polyfit`
