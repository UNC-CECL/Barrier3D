from setuptools import setup, find_packages


setup(
    name="barrier3d",
    version="2.0",
    description="A spatially explicit exploratory model of barrier island evolution",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Ian Reeves",
    author_email="ian.reeves@whoi.edu",
    url="https://github.com/UNC-CECL/Barrier3D",
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    install_requires=open("requirements.txt", "r").read().splitlines(),
    setup_requires=[],
    include_package_data=True,
    packages=find_packages(),
    entry_points={"console_scripts": ["b3d=barrier3d.cli:barrier3d"]},
)
