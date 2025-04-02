import os
import sys
import shutil
import subprocess

def get_dependencies(lib_path):
    """Returns a list of dependent shared libraries for a given library."""
    result = subprocess.run(["otool", "-L", lib_path], capture_output=True, text=True)
    lines = result.stdout.split("\n")[1:]  # Skip the first line (library itself)
    dependencies = []

    for line in lines:
        parts = line.strip().split(" ")
        if parts:
            dep_path = parts[0]
            if not dep_path.startswith("/usr/lib") and not dep_path.startswith("/System/"):
                dependencies.append(dep_path)

    return dependencies

def resolve_rpath(lib_path, dep):
    """Finds the real path for an @rpath dependency using otool -l."""
    result = subprocess.run(["otool", "-l", lib_path], capture_output=True, text=True)
    lines = result.stdout.split("\n")
    rpaths = []

    # Parse LC_RPATH entries
    for i in range(len(lines)):
        if "LC_RPATH" in lines[i]:
            path_line = lines[i + 2].strip().split("path ")[-1].split(" (offset")[0]
            rpaths.append(path_line)

    # Try to resolve @rpath/libXXX.dylib using known rpaths
    for rpath in rpaths:
        possible_path = os.path.join(rpath, os.path.basename(dep))
        if os.path.exists(possible_path):
            return possible_path

    return None  # Not found

def copy_and_fix_dependencies(source_lib, dest_dir):
    """Recursively copies the library and its dependencies to dest_dir and fixes paths."""
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    lib_dir = os.path.dirname(source_lib)  # Look for missing @rpath libraries here
    copied_libs = {}
    queue = [source_lib]
    processed_libs = set()  # Track already processed libraries
    copied_lib_paths = set()  # Track already copied libraries
    resolved_deps = set()  # Track already resolved dependencies

    while queue:
        lib_path = queue.pop()
        if lib_path in processed_libs or lib_path in resolved_deps or not os.path.isfile(lib_path):  # Skip if already processed, resolved or is a directory
            continue

        lib_name = os.path.basename(lib_path)
        dest_path = os.path.join(dest_dir, lib_name)

        if not os.path.exists(lib_path):
            # 🔹 Check if missing @rpath libraries exist in the same folder as source_lib
            alternative_path = os.path.join(lib_dir, lib_name)
            if os.path.exists(alternative_path) and os.path.isfile(alternative_path):
                print(f"🔍 Found missing @rpath library: {lib_path} -> {alternative_path}")
                lib_path = alternative_path
            else:
                print(f"⚠️ Warning: {lib_path} not found. Trying to resolve @rpath.")
                resolved_path = resolve_rpath(source_lib, lib_path)
                if resolved_path:
                    print(f"✅ Resolved {lib_path} -> {resolved_path}")
                    lib_path = resolved_path
                else:
                    print(f"❌ Could not resolve {lib_path}. Skipping.")
                    continue

        try:
            shutil.copy2(lib_path, dest_path)
            os.chmod(dest_path, 0o755)  # Ensure writable permissions
            copied_lib_paths.add(lib_path)
        except PermissionError:
            print(f"⚠️ Permission denied when copying {lib_path}. Try running with sudo.")
            continue

        copied_libs[lib_path] = dest_path
        processed_libs.add(lib_path)  # Mark this library as processed
        resolved_deps.add(lib_path)  # Mark this library as resolved

        # 🔹 Get dependencies and add to queue
        dependencies = get_dependencies(lib_path)
        for dep in dependencies:
            if dep.startswith("@rpath"):
                resolved_dep = resolve_rpath(lib_path, dep)
                if resolved_dep:
                    print(f"🔄 Resolved @rpath dependency: {dep} -> {resolved_dep}")
                    if resolved_dep not in resolved_deps and resolved_dep not in processed_libs and resolved_dep not in copied_lib_paths:
                        queue.append(resolved_dep)  # Add unresolved dependency to queue
                else:
                    print(f"⚠️ Could not resolve @rpath dependency: {dep}")
            else:
                if dep not in resolved_deps and dep not in processed_libs and dep not in copied_lib_paths:
                    queue.append(dep)  # Add resolved dependency to queue

    # 🔹 FIX: Ensure all dependencies are updated to @loader_path
    for original_path, new_path in copied_libs.items():
        new_name = os.path.basename(new_path)
        subprocess.run(["install_name_tool", "-id", f"@loader_path/{new_name}", new_path])

        for dep in get_dependencies(new_path):
            dep_name = os.path.basename(dep)
            if dep in copied_libs:
                new_dep_path = f"@loader_path/{dep_name}"
                print(f"🔄 Updating dependency: {new_path} -> {new_dep_path}")
                subprocess.run(["install_name_tool", "-change", dep, new_dep_path, new_path])

            # 🔹 Ensure @rpath dependencies are replaced
            elif dep.startswith("@rpath"):
                resolved_dep = resolve_rpath(original_path, dep)
                if resolved_dep and os.path.exists(resolved_dep) and os.path.isfile(resolved_dep):
                    print(f"🔄 Replacing @rpath: {new_path} -> {resolved_dep}")
                    subprocess.run(["install_name_tool", "-change", dep, f"@loader_path/{os.path.basename(resolved_dep)}", new_path])

    print("✅ All dependencies copied and fixed successfully.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python copy_runtime_dependencies.py <target_binary> <destination_directory>", file=sys.stderr)
        sys.exit(1)

    target_binary = sys.argv[1]
    destination_directory = sys.argv[2]
    copy_and_fix_dependencies(target_binary, destination_directory)