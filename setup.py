import os
# Always prefers setuptools over distutils
from setuptools import setup, find_packages

_VERSION = "1.0.0"


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as _in:
        return _in.read()


if __name__ == "__main__":
    setup(name="midasmlpy",
        version=_VERSION,
        description="Python wrapper for midasml",
        long_description=read('README.md'),
        author="Jonas Striaukas, Kris Stern, and Marcus Egelund-Müller",
        author_email="jonas.striaukas@gmail.com",
        url="https://github.com/jstriaukas/midasmlpy",
        install_requires=[
            "numpy>=2.0.0",
            "scipy>=1.14.0",
            "scikit-learn>=1.5.1",
            "pandas>=2.2.2",
            "openpyxl>=3.1.5"
        ],
        python_requires=">=3.10.0",
        setup_requires=["setuptools"],
        packages=find_packages(include=["midasmlpy"]),
        package_data={"midasmlpy.src.sparseglf90": ["sparsegllog_module.cpython-311-darwin.so"]},
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Operating System :: MacOS',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
            'Topic :: Scientific/Engineering'
        ],
    )
