import ez_setup
ez_setup.use_setuptools()


from setuptools import setup, find_packages

setup(
    name="bio_utils",
    version="0.0.1",
    packages=find_packages(),

    install_requires=['gseapy>=0.9.2'],

)
