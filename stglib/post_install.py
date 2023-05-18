import os
from pathlib import Path

def update_path():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    new_path = str(Path(package_dir).parents[3].joinpath('scripts'))

    # Modify the PATH environment variable
    os.environ['PATH'] = f"{new_path}:{os.environ['PATH']}"
