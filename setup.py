from setuptools import setup, find_packages

setup(
    name="WonderPeaks",
    version="0.1.6",
    description="Peak and UTR calling functions for ChIP-seq and RNA-seq in Fungi with WonderPeaks",
    author="Megan E. Garber",
    author_email="mgarber21@gmail.com",
    url="https://github.com/mgarber21/WonderPeaks",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,  # Enables MANIFEST.in
    install_requires=[],  # no dependencies
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
)