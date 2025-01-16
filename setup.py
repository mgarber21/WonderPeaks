from setuptools import setup, find_packages

# Dynamically read requirements.txt
with open("requirements.txt", "r") as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="WonderPeaks",
    version="0.1.0",
    packages=find_packages(),
    install_requires=requirements,  # Dependencies are loaded dynamically
)