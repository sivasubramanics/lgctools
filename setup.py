from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.readlines()

long_description = 'Python tool to process and perform various \
    LGC related tasks, mainly in focus to breeding applications.'

setup(
    name='lgctools',
    version='1.0.0',
    author='S.Sivasubramani',
    author_email='c.s.sivasubramani@gmail.com',
    url='https://github.com/siavsubramanics/lgctools',
    description='Processes the LGC files for tasks involved in purity check.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    packages=find_packages(),
    entry_points={
            'console_scripts': [
                'lgctools = lgctools.main:main'
            ]
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    keywords='geeksforgeeks gfg article python package vibhu4agarwal',
    install_requires=requirements,
    zip_safe=False
)
