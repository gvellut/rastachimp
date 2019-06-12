from setuptools import find_packages, setup

with open("rastachimp/__init__.py") as f:
    for line in f:
        if line.find("__version__") >= 0:
            version = line.split("=")[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            break

with open("README.md") as f:
    readme = f.read()

with open("requirements.txt") as f:
    requirements = f.readlines()

with open("requirements-dev.txt") as f:
    requirements_dev = f.readlines()

setup_args = dict(
    name="rastachimp",
    version=version,
    description="Tools for topological processing of polygons vectorized from rasters",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/gvellut/rastachimp",
    author="Guilhem Vellut",
    author_email="g@vellut.com",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: GIS",
    ],
    keywords="gis raster vector data",
    packages=find_packages(exclude=["docs", "tests"]),
    install_requires=requirements,
    extras_require={"dev": requirements_dev},
    project_urls={
        "Bug Reports": "https://github.com/gvellut/rastachimp/issues",
        "Source": "https://github.com/gvellut/rastachimp",
    },
)

setup(**setup_args)
