'''
Standard setup.py to upload the code on pypi.
First, remember to do
    bash generate_documentation.sh -all
to be sure docstrings and readme are up-to-date. Then
    python setup.py sdist
    python setup.py sdist upload
'''

from setuptools import setup

# Long description from readme
def readme():
    os.system('pandoc --from=markdown --to=rst --output=README.rst README.md') # Convert md to rst
    with open('README.rst') as f:
        return f.read()

# Extract version
def get_version():
    with open('gwdet/gwdet.py') as f:
        for line in f.readlines():
            if "__version__" in line:
                return line.split('"')[1]

setup(
    name='gwdet',
    version=get_version(),
    description='Detectability of gravitational-wave signals from compact binary coalescences',
    long_description=readme(),
    classifiers=[
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    keywords='gravitational-wave, black-hole binary',
    url='https://github.com/dgerosa/gwdet',
    author='Davide Gerosa',
    author_email='dgerosa@caltech.edu',
    license='MIT',
    packages=['gwdet'],
    install_requires=['numpy','scipy','scipy','matplotlib','astropy','pathos'],
    include_package_data=True,
    zip_safe=False,
)
