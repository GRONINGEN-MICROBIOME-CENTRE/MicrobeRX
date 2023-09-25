from setuptools import setup, find_packages
import codecs
import os
import versioneer

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

DESCRIPTION = "MicrobeRX is A tool for enzymatic reaction-based metabolite prediction in the gut microbiome."
# LONG_DESCRIPTION = ''

# Setting up
setup(
    name="MicrobeRX",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="GRONINGEN-MICROBIOME-CENTRE (Angel J. Ruiz-Moreno)",
    author_email="<angel.j.ruiz.moreno@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=[],
    keywords=[
        "python",
        "metabolite",
        "prediction",
        "microbiome",
        "cheminformatics",
        "metabolism",
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires=">=3.8",
    url="https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX",
    download_url="https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX",
    include_package_data=True,
    package_data={"": ["DataBase/*.tsv.gz"]},
)
