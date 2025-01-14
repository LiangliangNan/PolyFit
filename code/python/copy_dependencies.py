import os
import shutil
import sys
import subprocess
from typing import Set, List

# List of system paths to exclude from the PATH variable
EXCLUDED_PATHS_AND_FILES = [
    "c:\\windows",  # Example exclusion for Windows
    "api-ms-", "VCRUNTIME", "MSVCP",
    "/usr/lib",      # Common exclusion for Unix-like systems
    "/System",       # macOS system paths
    # Add more paths here as needed
]

def find_dependencies(target_binary: str, platform: str) -> Set[str]:
    """
    Finds runtime dependencies of a given binary using platform-specific tools.
    """
    dependencies = set()
    try:
        if platform == 'win32':  # Windows
            result = subprocess.run(['dumpbin', '/DEPENDENTS', target_binary], capture_output=True, text=True, check=True)
            dependencies.update(line.strip() for line in result.stdout.splitlines() if line.strip().endswith(".dll"))

        elif platform == 'darwin':  # macOS
            result = subprocess.run(['otool', '-L', target_binary], capture_output=True, text=True, check=True)
            dependencies.update(line.split()[0] for line in result.stdout.splitlines()[1:])

        elif platform.startswith('linux'):  # Linux
            result = subprocess.run(['ldd', target_binary], capture_output=True, text=True, check=True)
            dependencies.update(line.split("=>")[1].split("(")[0].strip() for line in result.stdout.splitlines() if "=>" in line)

    except Exception as e:
        print(f"Error finding dependencies for {target_binary}: {e}", file=sys.stderr)

    return dependencies

def find_in_search_paths(filename: str, search_paths: List[str]) -> str:
    """
    Search for a file in the given search paths.
    """
    for directory in search_paths:
        filepath = os.path.join(directory, filename)
        if os.path.exists(filepath):
            return filepath
    return None

def is_excluded_path_or_file(filepath: str, excluded_list: List[str]) -> bool:
    """
    Checks if the given filepath contains any of the strings in the excluded list. The check is case-insensitive.
    """
    filepath_lower = filepath.lower()
    return any(excl.lower() in filepath_lower for excl in excluded_list)

def copy_dependencies(target_binary: str, destination: str, resolved_dependencies: Set[str] = None):
    """
    Recursively copies the target binary and its transitive runtime dependencies to the destination directory,
    avoiding system libraries and excluded paths.
    """
    if resolved_dependencies is None:
        resolved_dependencies = set()

    platform = sys.platform
    os.makedirs(destination, exist_ok=True)

    # Add directory of the target binary to search paths
    search_paths = []
    if platform == 'win32':
        search_paths = os.environ['PATH'].split(os.pathsep)
    elif platform == 'darwin':
        search_paths = ['/System/Library/Frameworks']
    elif platform.startswith('linux'):
        search_paths = os.environ.get('LD_LIBRARY_PATH', '').split(os.pathsep)

    # Ensure the target binary is processed only once
    if target_binary in resolved_dependencies:
        return
    resolved_dependencies.add(target_binary)

    # Copy the main target binary
    try:
        shutil.copy2(target_binary, destination)
        print(f"Copied: {target_binary} -> {destination}")
    except Exception as e:
        print(f"Error copying {target_binary}: {e}", file=sys.stderr)

    # Find and process dependencies
    dependencies = find_dependencies(target_binary, platform)
    for dep in dependencies:
        if dep in resolved_dependencies:
            continue  # Skip already processed dependencies

        dep_path = find_in_search_paths(os.path.basename(dep), search_paths)
        if dep_path and not is_excluded_path_or_file(dep_path, EXCLUDED_PATHS_AND_FILES):
            copy_dependencies(dep_path, destination, resolved_dependencies)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python copy_runtime_dependencies.py <target_binary> <destination>", file=sys.stderr)
        sys.exit(1)

    target_binary = sys.argv[1]
    destination = sys.argv[2]

    if not os.path.exists(target_binary):
        print(f"Error: Target binary {target_binary} does not exist.", file=sys.stderr)
        sys.exit(1)

    copy_dependencies(target_binary, destination)
