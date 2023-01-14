from setuptools import setup, find_packages


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="AnnoSpat",
    version="1.0.0",
    author="Aanchal Mongia",
    author_email="aanchal.mongia@pennmedicine.upenn.edu",
    description="A tool for annotating and inferring patterns in single cells from spatial proteomics",
   long_description_content_type="text/markdown",
    url="https://github.com/aanchalMongia/AnnoSpat",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
#     python_requires=">=3.6",
#     packages=setuptools.find_packages(),
        entry_points={
            "console_scripts": ["AnnoSpat=AnnoSpat.run:main"]
    },
    #package_dir={"": "AnnoSpat"},
    packages=["AnnoSpat",'AnnoSpat.dependencies'],
    
    include_package_data= True,
    package_data={  
        '': ['data/*.csv', 'output_AnnoSpat/*.csv' ,  'models/model_ELM.pkl', 'models/model_kmeans.pkl',   'sample_Rscript.R', 'pp.R']
    },
        #'AnnoSpat.data': ['signatures_T1D.csv'],
        #'AnnoSpat.output_AnnoSpat': ['trte_labels_ELM_IMC_T1D_AnnoSpat.csv'],
        #'AnnoSpat.output_AnnoSpat': ['trte_labels_numericLabels_ELM_IMC_T1D_AnnoSpat.csv'],
   # },
    
    #packages=['AnnoSpat.src.data','AnnoSpat.src.dependencies','AnnoSpat.src.models'],
    python_requires=">=3.6",

    install_requires=["typer","numpy","scikit-learn","pandas"],
    
#       packages=['AnnoSpat'],
#       package_dir={'AnnoSpat': 'src/mypkg'},
#       package_data={'AnnoSpat': ['data/*.dat']},

    
)

