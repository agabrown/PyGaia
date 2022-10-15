from setuptools import find_packages, setup

setup(
    name="PyGaia",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "": ["LICENSE", "AUTHORS.md", "HISTORY.md", "INSTALL.md"],
        "pygaia.errors": ["data/*"],
    },
    include_package_data=True,
)
