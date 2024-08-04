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

    description="Toolkit for Machine Learning and Text Generation",
    long_description=long_description,
    license='BSD 3-Clause License',

    packages=setuptools.find_packages(),
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
    extras_require={
        'tensorflow-cpu': [
            'tensorflow>=1.10.0,<2.0',
            'tensorflow-probability>=0.3.0,<0.8.0'
        ],
        'tensorflow-gpu': [
            'tensorflow-gpu>=1.10.0,<2.0',
            'tensorflow-probability>=0.3.0,<0.8.0'
        ]
    },
    package_data={
        "texar": [
            "../bin/utils/multi-bleu.perl",
        ]
    },
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
