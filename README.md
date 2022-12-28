# AnnoSpat automatically annotates cell types and quantifies inter-cellular arrangements from spatial proteomics
Aanchal Mongia, Diane Sauders, Alvin C. Powers, Ali Naji, Gregory W. Schwartz and Robert B. Faryabi

## Introduction
`Annospat` is a tool for annotating and inferring patterns in single cells from spatial
proteomics. It uses neural network and point process algorithms to automatically identify cell types and quantify cell-cell spatial relationships in the absence of manual annotation. Using the canonical markers for each protein in the antibody panel , it predicts the cell type labels by rendering the similarity between the cells in the proteomic space.

## Dependencies
Python 3.6 and libraries:
> numpy, scikit-learn, pandas 


## Implementation
To implement `Annospat`, the script `classify_IMCcells.py` takes as input the path to the raw proteomics matrix (with cells on rows and proteins as columns) and a signature file holding the canonical protein markers. 

Sample run:
```bash
pip install AnnoSpat
```
```
AnnoSpat generateLabels --inputfile data/<file.csv> --markerfile data/signatures_T1D.csv --output dir output_Annospat
```
```

Usage: AnnoSpat generateLabels [OPTIONS]

  Generate cell type annotations

Options:
  -i, --inputfile TEXT            [required] #path ot input proteomics file
  -m, --markerfile TEXT           [required] #path to markers file
  -o, --outputdir TEXT            [required] #path to output dir
  -f, --firstprotein TEXT         [required] #first protein to pick in the proteomics file
  -l, --lastprotein TEXT          [required] #last protein to pick in the proteomics file
  -r, --roicol TEXT               [required] #ROI column
  -n, --colsToNeglect TEXT        [default: ] #proteins to negelct in the proteomics file
  -d, --diseasecol TEXT           [default: ] # T1D/control
  -s, --suffix TEXT               [default: _IMC_T1D_AnnoSpat] #suffix to be given
  -c, --classifier TEXT           [default: ELM] #classifier to use 
  -t, --thresholdmaxandunknown TEXT
                                  [default: [99.999, 70]] #thresholds for each protein default: b/w 99 and 99.9
  -a, --thresholdadaptive TEXT    [default: [99.5,99.5,99.5,99.5,99.9,
                                  99,99.5,  99,99,  99.5,99.9,
                                  99.9,99.9,99.9,  99.5,99.5]] #adaptive thresholds for each protein

  -b, --fileseparator TEXT        [default: ,] #file spearator    
    
```
