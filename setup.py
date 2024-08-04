import sys
import setuptools


long_description = '''
NA
'''

if sys.version_info < (3, 6):
    sys.exit('Python>=3.6 is required by Texar.')

setuptools.setup(
    name="mhodlr",
    version="0.0.1",
    url="https://github.com/chenxinye/mhodlr/",

    description="Matlab code for simulating mixed-precision HODLR matrix",
    long_description=long_description,
    license='BSD 3-Clause License',

    packages=['mhodlr'],
    platforms='any',

    install_requires=[
        'regex>=2018.01.10',
        'numpy<1.17.0',
        'pathlib>=1.0',
        'pyyaml',
        'requests',
        'funcsigs>=1.0.2',
        'sentencepiece>=0.1.8',
        'packaging'
    ],
    package_data={
        "mhodlr": [
            "mhodlr",
        ]
    },
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',    ],
)
