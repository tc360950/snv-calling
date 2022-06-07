import setuptools

setuptools.setup(
    name="generator-py",
    version="1.0.0",
    author="Tomasz Cakala",
    author_email="tc360950@gmail",
    description="Utility wrapper for CONET",
    url="https://github.com/tc360950/CONET",
    classifiers=[
        "Programming Language :: Python :: 3"
    ],
    packages=["generator"],
    install_requires=[
      'numpy>=1.22.3',
      'pandas>=1.4.2',
      'networkx>=2.8'
    ],
    python_requires=">=3.7",
)
