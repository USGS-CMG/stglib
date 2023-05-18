import os

def update_path():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    new_path = os.path.join(package_dir, '..', '..', '..', 'scripts')

    # Modify the PATH environment variable
    os.environ['PATH'] = f"{new_path}:{os.environ['PATH']}"
