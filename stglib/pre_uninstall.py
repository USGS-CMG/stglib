import os
from pathlib import Path

def remove_from_path():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    path_to_remove = str(Path(package_dir).parents[3].joinpath('scripts'))

    # Remove the package directory from the PATH environment variable
    paths = os.environ['PATH'].split(':')
    paths = [path for path in paths if path != path_to_remove]
    os.environ['PATH'] = ':'.join(paths)
