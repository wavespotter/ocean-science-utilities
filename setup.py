import setuptools

with open("README.md", "r") as file:
    readme_contents = file.read()

setuptools.setup(
    name="ocean-science-utilities",
    version="0.0.0",
    license="Apache 2 License",
    install_requires=[
        "pysofar>=0.1.13",
        "numpy",
        "pandas",
        "scipy",
        "numba",
        "xarray",
        "requests",
        "numba",
        "numba_progress",
    ],
    description="Python package to interact with Sofar wave data",
    long_description=readme_contents,
    long_description_content_type="text/markdown",
    author="Pieter Bart Smit",
    author_email="sofaroceangithubbot@gmail.com",
    url="https://github.com/sofarocean/ocean-science-utilities.git",
    package_dir={"": "src"},
    packages=setuptools.find_packages("src"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    project_urls={"Sofar Ocean Site": "https://www.sofarocean.com"},
    include_package_data=True,
    package_data={"": ["*.json"]},
)