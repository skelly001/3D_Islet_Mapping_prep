// ==============================================================================
// Segmented Cell ROI to Pixel Color Histogram Extraction
// ==============================================================================
// Script: 2a_1-cell_ROI_to_pixel_IJ_macro_lang.ijm
// Description: Will measure the histogram (256 bins) of an 8-bit color image
//              of MS pixels for all *segmented cell* ROIs and save the results table.
//              Assign segmented cell ROIs to pixels based on dominant color in the
//              ROI. All pixels have been assigned a different color. 
//              Note: Color pixel image (data/2-ROI_to_pixel/Slide61_mask_8bit_color.png)
//              and segmented cell ROIs (data/1-ImageJ_cell_ROI/ROI.zip) must 
//              already be open in ImageJ.
//
// Input: - data/2-ROI_to_pixel/Slide61_mask_8bit_color.png
//        - data/1-ImageJ_cell_ROI/ROI.zip
// Output: - data/2-ROI_to_pixel/1-ImageJ_cell_ROI_to_pixel_histograms/*.txt
// ==============================================================================

path = "data/2-ROI_to_pixel/1-ImageJ_cell_ROI_to_pixel_histograms/"

m = roiManager('count');
for (i = 0; i < m; i++) {
    roiManager('select', i);
    n = i;
    nBins = 256;
    run("Clear Results");
    row = 0;
    getHistogram(values, counts, nBins);
    for (j=0; j<nBins; j++) {
        setResult("Value", row, values[j]);
        setResult("Count", row, counts[j]);
        row++;
     }
    updateResults();
    saveAs("results", path+n+".txt");
    close("Results");
}
