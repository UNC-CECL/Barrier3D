from setuptools import setup, find_packages

# import versioneer


setup(
    name="barrier3d",
    description="A spatially explicit exploratory model of barrier island evolution in three dimensions",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Ian Reeves",
    author_email="irbreeves@gmail.com",
    url="https://github.com/anardek/Barrier3D",
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    install_requires=open("requirements.txt", "r").read().splitlines(),
    setup_requires=[],
    packages=find_packages(),
    # entry_points={
    #     "console_scripts": ["compact=compaction.cli:compact"],
    #     "landlab.components": ["Compact=compaction.landlab:Compact"],
    # },
    # version=versioneer.get_version(),
    # cmdclass=versioneer.get_cmdclass(),
)
