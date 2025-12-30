//This macro will measure the histogram (256 bins) of an 8-bit color image of pixels for all *pixel* ROIs and save the results table.
//Goal: Assign a color bin to each pixel to use as a key when matching segmented cell ROIs to pixels. All pixels have been assigned a different color.
//Note: Color pixel image ("data/2-ROI_to_pixel/Slide61_mask_8bit_color.png") & pixel ROIs ("data/2-ROI_to_pixel/Slide61_Pixel_Flu_RoiSet.zip") 
//must already be open in ImageJ

path = "data/2-ROI_to_pixel/2-ImageJ_pixel_ROI_histograms/"

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


