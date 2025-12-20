from ij.plugin.frame import RoiManager
from ij import IJ
from ij import WindowManager
import os
imp = IJ.getImage();
rm = RoiManager.getRoiManager();
# Assume a RoiManager is opened

for roi in RoiManager.getInstance():
	#m = getIndexesAsString()
	m = rm.getRoiIndex(roi);
	rm.select(m);
	IJ.run(imp, "Color Histogram", "");
	root = "C:/Users/kell343/OneDrive - PNNL/Documents/11 HuBMAP/Protein_Data/Image/Immunostaining/Fluorescence/Custom_composites/Fiji_mask/Slide61_mask_cell_boundaries/Imagej_py_Slide61_cropped_colors"
	m = str(m)
	suffix = "csv"
	path = os.path.join(root, m + '.' + suffix)
	IJ.saveAs("Results", path);
	WindowManager.closeAllWindows()
	
	imagePath = 'C:/Users/kell343/OneDrive - PNNL/Documents/11 HuBMAP/Protein_Data/Image/Immunostaining/Fluorescence/Custom_composites/Fiji_mask/Slide61_mask_cell_boundaries/Slide61-islet23_with_cell_boundaries_and_edges/Slide61-islet23_with_cell_bound_edge_crop/cropped - Copy.tif'
	imp = IJ.openImage(imagePath)
	imp.show(imagePath)
	
	
	
	
	#dataset = ij.io().open('C:/Users/kell343/OneDrive - PNNL/Documents/11 HuBMAP/Protein_Data/Image/Immunostaining/Fluorescence/Custom_composites/Fiji_mask/Slide61_mask_cell_boundaries/Slide61-islet23_with_cell_boundaries_and_edges/Slide61-islet23_with_cell_bound_edge_crop/cropped - Copy.tif')
	#ij.py.show(dataset)
	#print(WindowManager.getTempCurrentImage())
	#print(fail)
	#print(WindowManager.getNonImageWindows())
	#print(WindowManager.getNonImageTitles())
	#print(WindowManager.getAllNonImageWindows())
	#IJ.selectWindow("Results")
	#WindowManager.getWindow("Histogram of cropped")
	#imp.close()
	#imp.close()
	
