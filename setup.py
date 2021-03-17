import setuptools

setuptools.setup(
    name="mssuite", 
    version="0.1",
    author="Kevin Klann",
    author_email="klann@em.uni-frankfurt.de",
    description="Collection of data analysis tools",
    
    packages=setuptools.find_packages(),
    package_data={'':['../data/*']},
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)