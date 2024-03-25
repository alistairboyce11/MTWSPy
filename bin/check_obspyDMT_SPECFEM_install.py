import importlib

def check_package(package_name):
    try:
        importlib.import_module(package_name)
        print(f"{package_name} is installed and available.")
        return True
    except ImportError:
        print(f"{package_name} is not installed or not available.")
        return False

def main():
    # These are the pacakges required to produces processed SPECFEM and obspyDMT outputs...
    packages_to_check = ['obspy','numpy', 'scipy', 'pandas', 'matplotlib', 'sys', 'glob', 'shutil', 'os', 'time', 'warnings', 'datetime', 'inspect', 'yaml']

    all_packages_installed = all(check_package(package) for package in packages_to_check)

    if all_packages_installed:
        print(' ')
        print("All required packages are installed.")
        ERROR_CODE=0
    else:
        print(' ')
        print("Some packages are missing. Please install them before running SPECFEM/obspyDMT.")
        print('Make sure Python can import the following: ')
        print(str(packages_to_check))
        ERROR_CODE=1

    return ERROR_CODE

if __name__ == "__main__":
    out=main()
    print(out)
