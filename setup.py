from setuptools import find_packages, setup

setup(
    name="crystalprint",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "trimesh",
        "pymatgen",
        "pyglet<2",
    ],
    include_package_data=True,
    description="A Python package to generate 3D meshes from CIF files with bonds and atoms.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://rub.de",
    author="Pierre-Paul De Breuck",
    author_email="pierre-paul.debreuck@rub.de",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
