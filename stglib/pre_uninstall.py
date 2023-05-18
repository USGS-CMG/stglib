import os

def remove_from_path():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    path_to_remove = os.path.join(package_dir, '..', '..', '..', scripts')

    # Remove the package directory from the PATH environment variable
    paths = os.environ['PATH'].split(':')
    paths = [path for path in paths if path != path_to_remove]
    os.environ['PATH'] = ':'.join(paths)
