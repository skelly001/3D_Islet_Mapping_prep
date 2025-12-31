# ==============================================================================
# ROI RGB Histogram Extraction for Cell Type Classification
# ==============================================================================
# Script: 1a_1-ROI_hist_and_save.py
# Description: Fiji/ImageJ Python script that iterates through all ROIs in the
#              RoiManager and extracts RGB color histograms for each segmented
#              cell. Saves histogram data as CSV files for downstream cell type
#              classification based on immunofluorescence channel intensities.
#
# Input: - data/1-ImageJ_cell_ROI/cropped.tif (IF image)
#        - ROI.zip (loaded in RoiManager)
# Output: - data/1-ImageJ_cell_ROI/ROI_RGB_values/*.csv
# ==============================================================================

# Run script in Fiji (ImageJ)
from ij.plugin.frame import RoiManager
from ij import IJ
from ij import WindowManager
import os
imp = IJ.getImage();
rm = RoiManager.getRoiManager();
# Assume RoiManager is opened

for roi in RoiManager.getInstance():
	m = rm.getRoiIndex(roi);
	rm.select(m);
	IJ.run(imp, "Color Histogram", "");
	root = "data/1-ImageJ_cell_ROI/ROI_RGB_values"
	m = str(m)
	suffix = "csv"
	path = os.path.join(root, m + '.' + suffix)
	IJ.saveAs("Results", path);
	WindowManager.closeAllWindows()
	
	imagePath = 'data/1-ImageJ_cell_ROI/cropped.tif'
	imp = IJ.openImage(imagePath)
	imp.show(imagePath)
