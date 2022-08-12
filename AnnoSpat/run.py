#!/usr/bin/env python
# coding: utf-8

# In[ ]:




import typer
from AnnoSpat import classify_IMCcells

import os
from pathlib import Path

app = typer.Typer()

@app.command('generateLabels')
def anno(
    data_path: str = typer.Option(..., "--inputfile",'-i', help=""),
    signature_path: str = typer.Option(..., "--markerfile","-m", help=""),
    op_dir: str = typer.Option(..., "--outputdir", "-o",help=""),
    p1: str = typer.Option(..., "--firstprotein","-f",help=""),
    p2: str = typer.Option(..., "--lastprotein","-l",help=""),
    roi_name_col: str = typer.Option(..., "--roicol","-r",help=""),
    
    toDrop: str = typer.Option('', "--colsToNeglect","-n",help=""),
    disease_status_col: str = typer.Option('', "--diseasecol","-d",help=""),
    suffix: str = typer.Option('_IMC_T1D_AnnoSpat', "--suffix","-s",help=""),
    classifier: str = typer.Option('ELM', "--classifier","-c",help=""),
    th_str: str = typer.Option('[99.999, 70]', "--thresholdmaxandunknown", "-t",help=""),
    adaptive_th_vec_str: str = typer.Option('[99.5,99.5,99.5,99.5,99.9,  99,99.5,  99,99,  99.5,99.9,  99.9,99.9,99.9,  99.5,99.5]', "--thresholdadaptive","-a",help=""),
    fileseparator: str = typer.Option(',', "--fileseparator","-b",help="")
    
    
    
):
    """Generate cell type annotations
    """
    print("=====Function chosen: ANNOTATION=====")

    

    data_path = (Path(os.getcwd()) / data_path).resolve() 
    signature_path = (Path(os.getcwd()) / signature_path).resolve()
    op_dir = str( (Path(os.getcwd()) / op_dir).resolve() )
    
#    print("PATHS = ",data_path, signature_path, op_dir)
    
    #compulsory
    # data_path = '' #'data/mgDF.csv'
    # signature_path= '' #'data/signatures_T1D.csv'
    # op_dir = 'output_AnnoSpat'
    # p1 = 'HLA.ABCp' # first protein/column in protein expression file
    # p2 = 'Ghrelinp' # last protein/column in protein expression file
    # roi_name_col = 'TIFFfilename'

    #optional
    
#     suffix = '_IMC_T1D_AnnoSpat'
#     disease_status_col = '' # disease_status_col = 'Status'
#     classifier = 'ELM'#
#     th=[99.9, 99.999, 70] 
#     adaptive_th_vec=[99.5,99.5,99.5,99.5,99.9,  99,99.5,  99,99,  99.5,99.9,  99.9,99.9,99.9,  99.5,99.5] # 2 adaptive_th_vec

    #fixed
    sampling_ratio=0.5 # train test ratio (50% train ROIs labelled using semi supervised clustering, remaining test ROIs labelled by training classifier )

        
    # input modifications
    th_str_temp = th_str[1: len(th_str) -1].split(',')
    th = [float(item) for item in th_str_temp]

    adaptive_th_vec_str_temp = adaptive_th_vec_str[1: len(adaptive_th_vec_str) -1].split(',')
    adaptive_th_vec = [float(item) for item in adaptive_th_vec_str_temp]

           
    classify_IMCcells.anno(data_path, signature_path, op_dir, p1, p2, toDrop, roi_name_col,disease_status_col,classifier,th,adaptive_th_vec,   sampling_ratio,suffix,fileseparator)


    
@app.command('findPattern')
def spatial(
    ip_roi: str = typer.Option(..., "--input",'-i',help="Input ROI file "),
    out: str = typer.Option(..., "--out",'-o',help="Path of directory to save output files"),
    labels: str = typer.Option(..., "--labels",'-l',help="Predicted cell labels (optional)"),
    marks: str = typer.Option(..., "--marks",'-m',help="Comma separated list of Cell types")
):
    """ Finds patterns in across cell types using generated annotations
    """
    print("=====Function chosen: SPATIAL ANALYSIS=====")
    classify_IMCcells.spatial(marks, labels, out, ip_roi)
    


def main():
    #print("Please specify command: generateLabels / findPattern")
    app()

if __name__ == "__main__":
    app()

    
    #sudo Rscript pp.R --marks "Alpha_Cells,B_Cells" --labels /mnt/data1/gw/research/sc_integration/labels/labels_aanchal_imc.csv --out out_dir/ --input /mnt/data1/gw/research/sc_integration/data/imc/split_wide_data/20181130_HPAP003_1_C_ROI4.tif_mat.csv
