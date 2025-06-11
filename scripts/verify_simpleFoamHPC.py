import os
import shutil
import subprocess
import hashlib
import sys

def checksum_file(file_path):
    """Returns the SHA256 checksum of a file"""
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as f:
        while chunk := f.read(8192):
            sha256.update(chunk)
    return sha256.hexdigest()

def run_allrun_script(directory):
    """Executes the Allrun shell script in the provided directory"""
    allrun_path = os.path.join(directory, 'Allrun')
    if os.path.exists(allrun_path) and os.access(allrun_path, os.X_OK):
        subprocess.run([allrun_path], cwd=directory, check=True)
    else:
        print(f"Error: 'Allrun' script not found or not executable in {directory}")
        sys.exit(1)

def main(project_dir):
    tmp_dir = os.path.join(os.getcwd(), 'tmp')

    # Check if project directory exists
    if not os.path.isdir(project_dir):
        print(f"Error: Project directory {project_dir} does not exist")
        sys.exit(1)

    # Create tmp directory if it doesn't exist
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Define source and destination paths
    simple_car_src = os.path.join(project_dir, 'tutorials', 'incompressible', 'simpleFoam', 'simpleCar')
    simple_car_dst = os.path.join(tmp_dir, 'simpleCar')
    simple_car_hpc_dst = os.path.join(tmp_dir, 'simpleCarHPC')

    # Ensure source directory exists
    if not os.path.isdir(simple_car_src):
        print(f"Error: Source directory {simple_car_src} does not exist")
        sys.exit(1)

    # Copy simpleCar to tmp/simpleCar and tmp/simpleCarHPC
    shutil.copytree(simple_car_src, simple_car_dst)
    shutil.copytree(simple_car_src, simple_car_hpc_dst)

    # Run Allrun script in both directories
    try:
        run_allrun_script(simple_car_dst)
        run_allrun_script(simple_car_hpc_dst)
    except subprocess.CalledProcessError:
        print("Error: The Allrun script failed to execute.")
        shutil.rmtree(tmp_dir)  # Clean up the tmp directory
        sys.exit(1)

    # Checksum comparison
    file1 = os.path.join(simple_car_dst, '208', 'U')
    file2 = os.path.join(simple_car_hpc_dst, '208', 'U')

    if os.path.exists(file1) and os.path.exists(file2):
        checksum1 = checksum_file(file1)
        checksum2 = checksum_file(file2)

        if checksum1 == checksum2:
            print("The checksum of '208/U' files are identical.")
        else:
            print("Warning: The checksum of '208/U' files are different!")
    else:
        print(f"Error: One or both of the files '208/U' do not exist in the expected locations.")

    # Clean up the tmp directory
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <project_directory>")
        sys.exit(1)

    project_directory = sys.argv[1]
    main(project_directory)
