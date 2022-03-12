import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='Lypunov-Exponent-and-Correlation-Dimension',
    version='1.0.0',
    author='Sayan Nag',
    description='Installation of Lypunov-Exponent-and-Correlation-Dimension package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/sayannag/Lypunov-Exponent-and-Correlation-Dimension',
    project_urls = {
        "Bug Tracker": "https://github.com/sayannag/Lypunov-Exponent-and-Correlation-Dimension"/issues"
    },
    license='MIT',
    packages=['Lypunov-Exponent-and-Correlation-Dimension'],
    install_requires=['requests'],
)
