<img src="https://raw.githubusercontent.com/JEFworks-Lab/STalign/main/STalign_logos_fin.png">

`STalign` aligns spatial transcriptomics (ST) tissue sections to each other and to 3D atlases like the Allen Brain Atlas using algorithms that build upon diffeomorphic metric mapping. 

More information regarding the overall approach, methods and validations can be found in our publication:
<a href="https://www.nature.com/articles/s41467-023-43915-7">
<b>STalign: Alignment of spatial transcriptomics data using diffeomorphic metric mapping</b>
Kalen Clifton*, Manjari Anant*, Gohta Aihara, Lyla Atta, Osagie K Aimiuwu, Justus M Kebschull, Michael I Miller, Daniel Tward^, Jean Fan^
</a>

## Overview

<img src="https://jef.works/assets/papers/STalign_anim_small.gif">

STalign enables:
- alignment of single-cell resolution ST datasets within technologies
- alignment single-cell resolution ST datasets to histology
- alignment of ST datasets across technologies
- alignment of ST datasets to a 3D common coordinate framework 

## Installation & Import

### Installation using pip

This installation method is intended for users who sets up a Python environment without `pipenv`.

```
pip install --upgrade "git+https://github.com/JEFworks-Lab/STalign.git"
```

*All dependencies will be installed into your selected environment with the above command. Dependencies can be found in the requirements.txt file.*

### Installation using Pipfile from source

This installation method is intended for users who sets up a Python environment with `pipenv`. `pipenv` allows users to create and activate a virtual environment with all dependencies within the Python project. For more information and installation instructions for `pipenv`, see https://pipenv.pypa.io/en/latest/.

Fork and `git clone` the `STalign` github repository.

From the base directory of your local `STalign` git repo, create a `Pipfile.lock` file from `Pipfile` using:

```
pipenv install requests
```

> **_NOTE:_** Since `Pipfile.lock` is platform-dependent and different across operating systems, do not commit `Pipfile.lock` to the git repo if contributing to `STalign` or collaborating with other people.

Activate the virtual environment using:

```
pipenv shell
```

Deactivate the virtual environment using:

```
exit
```

### Import
To import STalign into your Python script, use: 

```python
from STalign import STalign
```



## Introduction to `curve_annotator.py`**

`STalign` relies on an iterative gradient descent method to align two images. However, this process can be slow, and alignment accuracy may depend on the initialization. To assist with this, the `curve_annotator.py` script is provided, allowing users to manually select landmark points to help initialize the alignment. The tutorial is in this link: https://jef.works/STalign/notebooks/xenium-xenium-alignment.html

### **How to Use `curve_annotator.py`**

To use this script, the images must first be saved as `.npz` files. This can be done using `numpy` as shown below:

```python
import numpy as np

# Optional: Save images as .npz files for landmark point selection
np.savez('../xenium_data/Xenium_Breast_Cancer_Rep1', x=XI, y=YI, I=I)
np.savez('../xenium_data/Xenium_Breast_Cancer_Rep2', x=XJ, y=YJ, I=J)
```

This will generate `.npz` files:

- `Xenium_Breast_Cancer_Rep1.npz`
- `Xenium_Breast_Cancer_Rep2.npz`

Once the `.npz` files are created, run the `curve_annotator.py` script to manually select landmark points:

```sh
python curve_annotator.py Xenium_Breast_Cancer_Rep1.npz
python curve_annotator.py Xenium_Breast_Cancer_Rep2.npz
```

This script will launch a graphical user interface (GUI) that enables users to select landmark points interactively. The selected points will be saved in the following files:

- `Xenium_Breast_Cancer_Rep1_points.npy`
- `Xenium_Breast_Cancer_Rep2_points.npy`

To read the saved landmark points, use:

```python
import numpy as np

# Load saved landmark points
pointsIList = np.load('../xenium_data/Xenium_Breast_Cancer_Rep1_points.npy', allow_pickle=True).tolist()
pointsJList = np.load('../xenium_data/Xenium_Breast_Cancer_Rep2_points.npy', allow_pickle=True).tolist()

print(pointsIList)
print(pointsJList)
```

These files store selected points in dictionary format, mapping each index to its corresponding coordinates.

The corresponding datasets' code has been uploaded in the five folders. Since the landmarks for each code are different, `pointsIList` and `pointsJList` can be used directly.

## **Saving 3D Results**

Each 3D result is saved separately in individual CSV files for better organization and accessibility. These files store 3D coordinate values and relevant data.

## Appendix 
1. Executed on GPU
2. The unit of runtime is minute(min) + seconds (s).
