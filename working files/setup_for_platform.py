import platform
import os
def build_for_linux():
    print("Starting build for linux")
    os.system("python3 linux/setup.py build_ext --inplace")
    os.system("mv *.so midasml2/")
def build_for_windows():
    print("Starting build for Windows")
    os.system("python win64/setup.py build_ext --inplace")
    os.system("move *.pyd midasml2/")
    os.system("copy win64\\lib_win64\\*.dll midasml2")
    
    
def build_for_mac():
    print("Starting build for Mac")
    os.system("rm -R build")
    os.system("rm midasml2/*.so")
    os.system("python3 mac/setup.py build_ext --inplace")
    os.system("mv *.so midasml2/")
    print("______________________Running example___________")
    os.system("python3 examples.py")
if platform.system() == 'Linux':
    build_for_linux()
elif platform.system() == 'Windows':
    build_for_windows()
else:
    build_for_mac()
