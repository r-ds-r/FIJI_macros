// Macro to classify cells based on nuclear content
// Batch processes .tif files in a folder with user-defined parameters

// ===== GET USER INPUTS =====
Dialog.create("Cell-Nucleus Classification Settings");
Dialog.addNumber("Cell Channel:", 4);
Dialog.addNumber("Nucleus Channel:", 3);
Dialog.addChoice("Cell Threshold Method:", newArray("Default", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"), "Otsu");
Dialog.addChoice("Nucleus Threshold Method:", newArray("Default", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"), "Triangle");
Dialog.show();

cellChannel = Dialog.getNumber();
nucleusChannel = Dialog.getNumber();
cellThreshold = Dialog.getChoice();
nucleusThreshold = Dialog.getChoice();

// ===== SELECT INPUT FOLDER =====
inputDir = getDirectory("Choose a folder containing .tif files");
outputDir = inputDir + "Results" + File.separator;
File.makeDirectory(outputDir);

// Get list of .tif files
fileList = getFileList(inputDir);
tifFiles = newArray();
for (i = 0; i < fileList.length; i++) {
    if (endsWith(fileList[i], ".tif") || endsWith(fileList[i], ".tiff")) {
        tifFiles = Array.concat(tifFiles, fileList[i]);
    }
}

if (tifFiles.length == 0) {
    exit("No .tif files found in the selected folder!");
}

print("\\Clear");
print("===== BATCH PROCESSING START =====");
print("Found " + tifFiles.length + " .tif file(s)");
print("Cell Channel: " + cellChannel);
print("Nucleus Channel: " + nucleusChannel);
print("Cell Threshold: " + cellThreshold);
print("Nucleus Threshold: " + nucleusThreshold);
print("Output Directory: " + outputDir);
print("");

// ===== PROCESS EACH FILE =====
setBatchMode(true);

for (fileIdx = 0; fileIdx < tifFiles.length; fileIdx++) {
    // Open image
    open(inputDir + tifFiles[fileIdx]);
    originalImage = getTitle();
    baseName = File.nameWithoutExtension;
    
    print("Processing [" + (fileIdx+1) + "/" + tifFiles.length + "]: " + originalImage);
    
    // ===== STEP 1: Detect Cells =====
    selectWindow(originalImage);
    run("Duplicate...", "title=cells duplicate channels=" + cellChannel);
    run("8-bit");
    
    // Apply user-selected threshold
    setAutoThreshold(cellThreshold + " dark");
    setOption("BlackBackground", true);
    
    // Detect cells
    run("Analyze Particles...", "size=50.00-5000.00 show=Masks display");
    rename("Cell_Masks");
    run("Invert LUT");
    
    // Add cell masks to ROI Manager
    selectWindow("Cell_Masks");
    roiManager("reset");
    run("Analyze Particles...", "size=50.00-5000.00 add");
    run("Select None");
    
    // Store number of cells detected
    nCells = roiManager("count");
    print("  - Total cells detected: " + nCells);
    
    // ===== STEP 2: Detect Nuclei =====
    selectWindow(originalImage);
    run("Duplicate...", "title=nuclei duplicate channels=" + nucleusChannel);
    run("8-bit");
    
    // Threshold with user-selected method and convert to mask
    setAutoThreshold(nucleusThreshold + " dark");
    run("Convert to Mask");
    
    // Process nuclei: fill holes, watershed, analyze particles
    run("Fill Holes");
    run("Watershed");
    run("Analyze Particles...", "size=10-Infinity show=Nothing add");

// Store nuclei ROIs separately (they start after cell ROIs in ROI Manager)
nNuclei = roiManager("count") - nCells;
    print("  - Total nuclei detected: " + nNuclei);
    
    // ===== STEP 3: Classify Cells as Positive or Negative =====
    positiveCells = newArray();
    negativeCells = newArray();
    nPositive = 0;
    nNegative = 0;
    
    selectWindow(originalImage);
    for (i = 0; i < nCells; i++) {
        hasNucleus = false;
        
        // Check if this cell fully contains any nucleus
        for (j = nCells; j < roiManager("count"); j++) {
            // Get nucleus area
            roiManager("select", j);
            getStatistics(nucleusArea);
            
            // Get intersection area between cell and nucleus
            roiManager("select", newArray(i, j));
            roiManager("AND");
            
            // Check if nucleus is FULLY inside cell (intersection = full nucleus area)
            if (selectionType() != -1) {
                getStatistics(overlapArea);
                // Use tolerance for rounding errors (99.9% overlap = full containment)
                if (overlapArea >= nucleusArea * 0.999) {
                    hasNucleus = true;
                    break; // Found a fully enclosed nucleus, stop checking
                }
            }
        }
        
        // Classify based on nucleus presence
        if (hasNucleus) {
            positiveCells = Array.concat(positiveCells, i);
            nPositive++;
        } else {
            negativeCells = Array.concat(negativeCells, i);
            nNegative++;
        }
    }
    
    run("Select None");
    
    // ===== STEP 4: Save ROI Sets =====
    // First, save all ROIs (cells + nuclei) for reference
    roiManager("deselect");
    roiManager("save", outputDir + baseName + "_All_ROIs_with_Nuclei.zip");
    
    // Remove all nuclei ROIs (indices nCells and above) before saving cell-only ROIs
    nucleiToDelete = newArray();
    for (i = nCells; i < roiManager("count"); i++) {
        nucleiToDelete = Array.concat(nucleiToDelete, i);
    }
    if (nucleiToDelete.length > 0) {
        roiManager("select", nucleiToDelete);
        roiManager("delete");
    }
    
    // Now ROI Manager only contains cell ROIs (0 to nCells-1)
    // Save all cell ROIs
    roiManager("deselect");
    roiManager("save", outputDir + baseName + "_All_Cell_ROIs.zip");
    
    // Save positive cells - delete from highest index to lowest to preserve indices
    if (nPositive > 0 && nNegative > 0) {
        // Sort negative cells in descending order and delete one by one
        Array.sort(negativeCells);
        for (i = negativeCells.length - 1; i >= 0; i--) {
            roiManager("select", negativeCells[i]);
            roiManager("delete");
        }
        // Save remaining (positive) cells
        roiManager("deselect");
        roiManager("save", outputDir + baseName + "_Positive_Cells_ROIs.zip");
        
        // Restore all cells for negative save
        roiManager("reset");
        roiManager("open", outputDir + baseName + "_All_Cell_ROIs.zip");
        
        // Sort positive cells in descending order and delete one by one
        Array.sort(positiveCells);
        for (i = positiveCells.length - 1; i >= 0; i--) {
            roiManager("select", positiveCells[i]);
            roiManager("delete");
        }
        // Save remaining (negative) cells
        roiManager("deselect");
        roiManager("save", outputDir + baseName + "_Negative_Cells_ROIs.zip");
    } else if (nPositive > 0) {
        // All cells are positive
        roiManager("deselect");
        roiManager("save", outputDir + baseName + "_Positive_Cells_ROIs.zip");
    } else if (nNegative > 0) {
        // All cells are negative
        roiManager("deselect");
        roiManager("save", outputDir + baseName + "_Negative_Cells_ROIs.zip");
    }
    
    // ===== STEP 5: Print Summary for This File =====
    print("  - Positive cells (with nucleus): " + nPositive);
    print("  - Negative cells (without nucleus): " + nNegative);
    print("  - ROIs saved to: " + outputDir);
    print("");
    
    // Clean up
    close("*");
    roiManager("reset");
}

setBatchMode(false);

// ===== FINAL SUMMARY =====
print("===== BATCH PROCESSING COMPLETE =====");
print("Processed " + tifFiles.length + " file(s)");
print("Results saved to: " + outputDir);
