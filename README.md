# Automatic Astrometry and Photometry Fitting for Intensity VLBI Maps (AAP-VLBI-Fit)

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)  
![License](https://img.shields.io/badge/License-MIT-green)  

A Python-based tool for fitting Gaussian functions in Very Long Baseline Interferometry (VLBI) intensity maps. This routine retrieves both **astrometry and photometry** for each fitted component, offering an approach distinct from traditional methods that focus solely on astrometry.  

This tool was developed as part of the research described in:  
📄 **Amador, A. (2025), "AAP-VLBI-Fit: Gaussian Fitting in Intensity VLBI Maps," *Ark-Xiv*.** [DOI:]  

---

## **Features**
✅ Fits components in VLBI intensity maps  
✅ Extracts **astrometric and photometric** properties  
✅ Implements **adaptive thresholding** for source detection  
✅ Provides **residual analysis** and model comparison  
✅ Includes **Jupyter notebooks** for demonstration and usage  

---

## **Installation**  
### **1. Clone the Repository**  
```bash
git clone https://github.com/your-username/vlbi-auto-photometry-tool.git
cd vlbi-auto-photometry-tool
```
### **2. Install Dependencies**  
```bash
pip install -r requirements.txt
```
Or manually install:
```bash 
pip install numpy astropy photutils matplotlib
```

### **Usage**  
The core module with the fitting functions are declared in 'main_function.py'. However any of those function can be called as:
```bash 
from main.main_functions import any_function  # Example function
fit_function("../your_path/your_map.fits")
```
Jupyter Notebook Examples
To explore interactive examples, check the notebooks files:
main_usage.ipynb: Guide on using main_functions.py
demo.ipynb: Lightweight demo with a small VLBI map
#### Important Notes
Please provide the files as follows: Each map (.fits/.IMAP) in its own folder, so that you have a list with the path+file for each map as well as another list with only the path where each map is located.
We strongly suggest only fitting maps from the same source automatically. This is because the same detection limit may not necessarily be appropriate for each source.

### **Repository Structure** 
```bash 
aap-vlbi-fit-tool/
├── main/ 				    # Core Python module
│   └── fit_components.py   # VLBI component fitting functions
├── notebooks/ 			    # Jupyter notebooks for demos and usage
│   ├── main_usage.ipynb    # Main usage
│   └── demo.ipynb 		    # Demo
├── data/ 				    # Sample VLBI maps for the demo (optional)
│   └── sources/maps.fits 	# Example files
├── README.md 		        # Documentation
├── LICENSE 			    # License (MIT)
├── requirements.txt 		# Required Python packages
└── .gitignore 			    # Ignored files
```

### **Citing This Work**
If you use this software in your research, please cite:
📄 Amador-Portes et al., (2025b), "AAP-VLBI-Fit: Gaussian Fitting in Intensity VLBI Maps," Ark-Xiv. [DOI:]
Alternatively, you can use this BibTeX entry:
@article{Amador2025,
}

### **License**
This project is licensed under the MIT License. See the LICENSE file for details

### **Contributing**
Contributions are greatly welcome! Please fork the repository and submit a pull request.
For any questions or issues, feel free to open a GitHub issue or contact me at alfreportess9730@gmail.com