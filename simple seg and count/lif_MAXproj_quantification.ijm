// Macro to process .lif files with Z-stacks
// Creates max projections, quantifies cells in channel 1
// Exports ROI masks and results CSV

// Clean up any open images
run("Close All");
roiManager("reset");
run("Clear Results");

// Prompt user to select .lif file
lifPath = File.openDialog("Select .lif file");
lifDir = File.getParent(lifPath);
lifName = File.getName(lifPath);
lifBaseName = replace(lifName, ".lif", "");

// Create output directory for ROI masks
outputDir = lifDir + File.separator + lifBaseName + "_output" + File.separator;
File.makeDirectory(outputDir);

// Get number of series (images) in the .lif file
run("Bio-Formats Macro Extensions");
Ext.setId(lifPath);
Ext.getSeriesCount(seriesCount);

// Initialize results arrays
imageNames = newArray(seriesCount);
cellCounts = newArray(seriesCount);
gfapAreas = newArray(seriesCount);
gfapPercentages = newArray(seriesCount);

// Process each series
for (i = 0; i < seriesCount; i++) {
    
    // Open the current series
    run("Bio-Formats Importer", "open=[" + lifPath + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (i+1));
    
    // Get the current image title
    fullTitle = getTitle();
    print("Processing: " + fullTitle);
    
    // Extract image name (after " - ")
    // Format is typically "filename.lif - ImageName"
    if (indexOf(fullTitle, " - ") > 0) {
        imageName = substring(fullTitle, indexOf(fullTitle, " - ") + 3);
    } else {
        imageName = fullTitle;
    }
    
    // Check if image name matches criteria
    // Must contain "s100b" and must NOT contain "_Lng", "_DCV", or "_HD"
    imageNameLower = toLowerCase(imageName);
    hasS100b = indexOf(imageNameLower, "s100b") >= 0;
    hasLng = indexOf(imageName, "_Lng") >= 0;
    hasDCV = indexOf(imageName, "_DCV") >= 0;
    hasHD = indexOf(imageName, "_HD") >= 0;
    
    if (!hasS100b || hasLng || hasDCV || hasHD) {
        print("Skipping: " + imageName + " (does not match filter criteria)");
        close("*");
        imageNames[i] = "";
        cellCounts[i] = 0;
        gfapAreas[i] = 0;
        gfapPercentages[i] = 0;
        continue;
    }
    
    // Store clean image name
    imageNames[i] = imageName;
    
    // Get image properties
    getDimensions(width, height, channels, slices, frames);
    
    // Create max projection if Z-stack
    if (slices > 1) {
        run("Z Project...", "projection=[Max Intensity]");
        selectWindow(fullTitle);
        close();
        selectWindow("MAX_" + fullTitle);
        rename(imageName + "_MAX");
    } else {
        rename(imageName + "_MAX");
    }
    
    currentImage = getTitle();
    
    // Store max projection image name for later use
    maxImage = currentImage;
    
    // === PROCESS CHANNEL 1: Cell counting ===
    
    // Extract channel 1 if multiple channels exist
    if (channels > 1) {
        selectWindow(maxImage);
        run("Duplicate...", "duplicate channels=1");
        ch1Image = getTitle();
        rename(imageName + "_ch1");
    } else {
        selectWindow(maxImage);
        rename(imageName + "_ch1");
    }
    
    workingImage = getTitle();
    
    // Ensure 8-bit
    run("8-bit");
    
    // Apply Otsu thresholding
    setAutoThreshold("Otsu dark");
    run("Convert to Mask");
    
    // Clean up the mask
    run("Despeckle");
    run("Fill Holes");
    
    // Watershed to separate touching cells
    run("Watershed");
    
    // Save the mask/ROI image
    saveAs("Tiff", outputDir + imageName + "_mask.tif");
    maskImage = getTitle();
    
    // Make sure we're working with a binary mask
    run("8-bit");
    setThreshold(1, 255);
    run("Convert to Mask");
    
    // Clear previous results
    run("Clear Results");
    
    // Analyze particles to count cells with area filter (50-150 µm²)
    // ImageJ will interpret the size in calibrated units if the image has calibration
    run("Analyze Particles...", "size=50-150 circularity=0.0-1.0 show=Overlay display exclude clear add");
    
    // Get the count from Results table (number of rows)
    cellCount = nResults;
    cellCounts[i] = cellCount;
    
    print("  Channel 1 - Cells: " + cellCount);
    
    // Close channel 1 images
    roiManager("reset");
    close("*");
    
    // === PROCESS CHANNEL 2: GFAP area measurement ===
    
    // Re-open the max projection and extract channel 2
    if (channels > 1) {
        // Re-open the series to get the max projection again
        run("Bio-Formats Importer", "open=[" + lifPath + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (i+1));
        fullTitle = getTitle();
        
        // Create max projection if Z-stack
        if (slices > 1) {
            run("Z Project...", "projection=[Max Intensity]");
            selectWindow(fullTitle);
            close();
            maxTitle = "MAX_" + fullTitle;
            selectWindow(maxTitle);
        }
        
        // Extract channel 2
        run("Duplicate...", "duplicate channels=2");
        ch2Image = getTitle();
        rename(imageName + "_ch2");
        close(maxTitle);
        
        selectWindow(imageName + "_ch2");
        
        // Ensure 8-bit
        run("8-bit");
        
        // Apply Triangle thresholding
        setAutoThreshold("Triangle dark");
        run("Convert to Mask");
        
        // Clean up the mask (no watershed for GFAP)
        run("Despeckle");
        run("Fill Holes");
        
        // Save the GFAP mask
        saveAs("Tiff", outputDir + imageName + "_ch2_GFAP_mask.tif");
        
        // Measure area
        // Get image dimensions
        getDimensions(w, h, ch, sl, fr);
        totalPixels = w * h;
        
        // Measure histogram to count white pixels
        getHistogram(values, counts, 256);
        whitePixels = counts[255];
        
        // Calculate area and percentage
        gfapArea = whitePixels;
        gfapPercent = (whitePixels / totalPixels) * 100;
        
        gfapAreas[i] = gfapArea;
        gfapPercentages[i] = gfapPercent;
        
        print("  Channel 2 - GFAP Area: " + gfapArea + " pixels (" + d2s(gfapPercent, 2) + "%)");
        
        // Close channel 2 images
        close("*");
    } else {
        // No channel 2, set to 0
        gfapAreas[i] = 0;
        gfapPercentages[i] = 0;
        print("  Channel 2 - Skipped (single channel image)");
    }
    
    print("");
    
}

// Create final results table (only for processed images)
run("Clear Results");
resultRow = 0;
for (i = 0; i < seriesCount; i++) {
    if (imageNames[i] != "") {
        setResult("Image", resultRow, imageNames[i]);
        setResult("Cell_Count", resultRow, cellCounts[i]);
        setResult("GFAP_Area_pixels", resultRow, gfapAreas[i]);
        setResult("GFAP_Percent", resultRow, gfapPercentages[i]);
        resultRow++;
    }
}
updateResults();

// Save results to CSV
saveAs("Results", outputDir + lifBaseName + "_cell_counts.csv");

// Close Summary window if it exists
if (isOpen("Summary")) {
    selectWindow("Summary");
    run("Close");
}

print("\n=== Processing Complete ===");
processedCount = 0;
for (i = 0; i < seriesCount; i++) {
    if (imageNames[i] != "") {
        processedCount++;
    }
}
print("Processed " + processedCount + " of " + seriesCount + " images (filtered for s100b, excluding _Lng, _DCV, and _HD)");
print("Results saved to: " + outputDir);
print("CSV file: " + lifBaseName + "_cell_counts.csv");

// Show output folder
File.openDialog("Processing complete! Results saved in output folder.");
