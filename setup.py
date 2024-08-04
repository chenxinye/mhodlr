import sys
import setuptools


long_description = '''
NA
'''
setuptools.setup(
    name="mhodlr",
    version="0.0.1",
    url="https://github.com/chenxinye/mhodlr/",

    description="Matlab code for simulating mixed-precision HODLR matrix",
    long_description=long_description,
    license='BSD 3-Clause License',

    packages=['mhodlr'],
    platforms='any',

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
