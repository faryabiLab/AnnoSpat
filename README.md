# AnnoSpat automatically annotates cell types and quantifies inter-cellular arrangements from spatial proteomics
Aanchal Mongia, Diane Sauders, Alvin C. Powers, Ali Naji, Gregory W. Schwartz and Robert B. Faryabi

## Introduction
`Annospat` is a tool for annotating and inferring patterns in single cells from spatial
proteomics. It uses neural network and point process algorithms to automatically identify cell types and quantify cell-cell spatial relationships in the absence of manual annotation. Using the canonical markers for each protein in the antibody panel , it predicts the cell type labels by rendering the similarity between the cells in the proteomic space.




## Implementation
`Annospat` can be used as a python library or as a command line tool using the argument `generateLabels` which takes as input the path to the raw proteomics matrix (with cells on rows and proteins as columns) and a signature file holding the canonical protein markers with other data related inputs. 
The tool is available for installation at [https://pypi.org/project/AnnoSpat/1.0.0/](https://pypi.org/project/AnnoSpat/1.0.0/)


## Installation
Dependencies: Python 3.6+, pip 
Libraries->typer, numpy, scikit-learn, pandas 

> pip install AnnoSpat==1.0.0


From your working directory, execute:
```
mkdir outputdir
```
<!--
 AnnoSpat generateLabels -i /mnt/data2/aanchal/data/IMC_T1D/raw_data/mgDF.csv -m /mnt/data2/aanchal/data/IMC_T1D/signatures_T1D_withImmuneCelltypes_withNegMarkers_d.csv -o delete_outputdir -f 'HLA.ABC' -l 'Ghrelin' -r 'TIFFfilename' -t '[99.9,99.999,70]' -a '[99.5,99.5,99.5,99.5,99.9,99,99.5,99,99,99.5,99.9,99.9,99.9,99.9,99.5,99.9]' -d 'Status'

python3 run.py generateLabels -i /mnt/data2/aanchal/data/IMC_T1D/raw_data/mgDF.csv -m /mnt/data2/aanchal/data/IMC_T1D/signatures_T1D_withImmuneCelltypes_withNegMarkers_d.csv -o delete_outputdir -f 'HLA.ABC' -l 'Ghrelin' -r 'Status'
-->
```
AnnoSpat generateLabels -i <path_to_proteomics_matrix> -m <path_to_marker_file> -o outputdir -f <first_protein_name> -l <last_protein_name> -r <name_of_col_holding_ROInames>
```
Please replace the arguments to --inputfile/-i argument, --markerfile/-m argument and other arguments as per your own paths to proteomics and marker files and data.

```
Usage: AnnoSpat generateLabels [OPTIONS]

  Generate cell type annotations

Options:
  -i, --inputfile TEXT            [required] #path ot input proteomics file
  -m, --markerfile TEXT           [required] #path to marker file
  -o, --outputdir TEXT            [required] #path to output dir
  -f, --firstprotein TEXT         [required] #first protein to pick in the proteomics file
  -l, --lastprotein TEXT          [required] #last protein to pick in the proteomics file
  -r, --roicol TEXT               [required] #ROI column
  -n, --colsToNeglect TEXT        [default: ] #proteins to negelct in the proteomics file
  -d, --diseasecol TEXT           [default: ] # T1D/control
  -s, --suffix TEXT               [default: _Data_AnnoSpat] #suffix to be given
  -c, --classifier TEXT           [default: ELM] #classifier to use 
  -t, --thresholdmaxandunknown TEXT
                                  [default: [99.999, 70]] #thresholds for each protein default: b/w 99 and 99.9
  -a, --thresholdadaptive TEXT    [default: [99.5,99.5,99.5,99.5,99.9,
                                  99,99.5,  99,99,  99.5,99.9,
                                  99.9,99.9,99.9,  99.5,99.5]] #adaptive thresholds for each protein

  -b, --fileseparator TEXT        [default: ,] #file spearator    
```
NOTE: The number of values in  "adaptive_th_vec_str" is assumed to be the same as no of cell types as each value by definition correspond to the threshold to be taken for each cell type.

<!--
## Edit source 



> git clone https://github.com/faryabiLab/AnnoSpat.git

> cd AnnoSpat


> mkdir outputdir

> python3 run.py generateLabels -i <path_to_proteomics_matrix> -m <path_to_marker_file> -o outputdir -f <first_protein_name> -l <last_protein_name> -r <name_of_col_holding_ROInames>
-->




Once the cell type annotations are obtained, The neighborhood analysis with point process can be implemented as shown in as described at: [Too-Many-Cells](https://gregoryschwartz.github.io/too-many-cells/#spatial)

