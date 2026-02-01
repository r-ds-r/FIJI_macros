// Macro: Single Image Analysis with Cellpose-SAM Segmentation and Sholl Analysis
// Process a single image from .lif file with:
//   - ROI Marker channel (e.g., sox9, S100b, etc.) - segmented with cellpose-sam
//   - TdTomato channel (optional)
//   - GFAP channel for Sholl analysis (optional)
// Outputs: ROI+ cell count, TdTomato overlap, and SNT Sholl analysis
// Channel configuration is user-selectable at runtime

// ===== SETUP =====

// Clean up
run("Close All");
roiManager("reset");
run("Clear Results");

// Get user inputs
lifPath = File.openDialog("Select .lif file");
lifDir = File.getParent(lifPath);
lifName = File.getName(lifPath);
lifBaseName = replace(lifName, ".lif", "");

// Ask for cellpose-sam model path
modelPath = File.openDialog("Select trained cellpose-sam model file");

// Get analysis timestamp
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second);
month++; // Month is 0-indexed
analysisDate = "" + year + "-" + IJ.pad(month, 2) + "-" + IJ.pad(dayOfMonth, 2);
analysisTime = "" + IJ.pad(hour, 2) + ":" + IJ.pad(minute, 2) + ":" + IJ.pad(second, 2);
analysisTimestamp = analysisDate + " " + analysisTime;

// Get model filename
modelName = File.getName(modelPath);

// Ask user to specify which channel contains each marker
Dialog.create("Channel Configuration");
Dialog.addMessage("Please specify which channel contains each marker:");
Dialog.addChoice("ROI Marker Channel (for segmentation):", newArray("1", "2", "3", "4"), "1");
Dialog.addChoice("TdTomato Channel:", newArray("None", "1", "2", "3", "4"), "2");
Dialog.addChoice("GFAP Channel (for Sholl analysis):", newArray("None", "1", "2", "3", "4"), "3");
Dialog.addMessage("Select 'None' if that marker is not present in your images.");
Dialog.show();

roiChannel = parseInt(Dialog.getChoice());
tdTomatoChannelStr = Dialog.getChoice();
gfapChannelStr = Dialog.getChoice();

// Parse channel selections
if (tdTomatoChannelStr == "None") {
    tdTomatoChannel = -1;
    hasTdTomato = false;
} else {
    tdTomatoChannel = parseInt(tdTomatoChannelStr);
    hasTdTomato = true;
}

if (gfapChannelStr == "None") {
    gfapChannel = -1;
    hasGFAP = false;
} else {
    gfapChannel = parseInt(gfapChannelStr);
    hasGFAP = true;
}

// Create output directories
outputDir = lifDir + File.separator + lifBaseName + "_single_output" + File.separator;
tempDir = outputDir + "temp" + File.separator;
roiDir = outputDir + "rois" + File.separator;
shollDir = outputDir + "sholl_analysis" + File.separator;
File.makeDirectory(outputDir);
File.makeDirectory(tempDir);
File.makeDirectory(roiDir);
File.makeDirectory(shollDir);

// Copy model file to output directory for reproducibility
modelCopyPath = outputDir + modelName;
File.copy(modelPath, modelCopyPath);

print("\n=== Single Image Analysis with Cellpose-SAM ===");
print("Analysis Date/Time: " + analysisTimestamp);
print("LIF file: " + lifName);
print("Output directory: " + outputDir);
print("Model: " + modelName);
print("Model archived to: " + modelCopyPath);
print("Channel Configuration:");
print("  ROI Marker: Channel " + roiChannel);
if (hasTdTomato) {
    print("  TdTomato: Channel " + tdTomatoChannel);
} else {
    print("  TdTomato: Not present");
}
if (hasGFAP) {
    print("  GFAP: Channel " + gfapChannel);
} else {
    print("  GFAP: Not present");
}
print("");

// ===== STEP 1: OPEN IMAGE WITH BIO-FORMATS DIALOG =====

print("Opening Bio-Formats importer dialog...");
print("Please select the image you want to analyze\n");

// Open Bio-Formats with dialog to let user choose the series
run("Bio-Formats Importer", "open=[" + lifPath + "] color_mode=Default view=Hyperstack stack_order=XYCZT");

// Get the opened image
fullTitle = getTitle();
print("--- Processing Image: " + fullTitle + " ---");

// Extract image name
if (indexOf(fullTitle, " - ") > 0) {
    imageName = substring(fullTitle, indexOf(fullTitle, " - ") + 3);
} else {
    imageName = fullTitle;
}
imageName = replace(imageName, "/", "_"); // Remove problematic characters

// Get image properties and calibration
getDimensions(width, height, channels, slices, frames);
getVoxelSize(voxelWidth, voxelHeight, voxelDepth, unit);
print("Dimensions: " + width + "x" + height + ", Channels: " + channels + ", Slices: " + slices);
print("Voxel size: XY=" + d2s(voxelWidth, 4) + " " + unit + ", Z=" + d2s(voxelDepth, 4) + " " + unit);

// Calculate anisotropy for Cellpose 3D
anisotropy = voxelDepth / voxelWidth;
print("Calculated anisotropy: " + d2s(anisotropy, 2));

// Rename for easier handling
rename(imageName + "_original");
originalImage = getTitle();

// ===== STEP 2: EXTRACT AND SAVE ROI MARKER CHANNEL FOR CELLPOSE =====

print("Extracting ROI marker channel (channel " + roiChannel + ") - full z-stack...");
selectWindow(originalImage);
if (channels > 1) {
    run("Duplicate...", "duplicate channels=" + roiChannel);
    roiMarkerImage = getTitle();
    rename(imageName + "_roiMarker");
} else {
    run("Duplicate...", "duplicate");
    rename(imageName + "_roiMarker");
}

roiMarkerImage = getTitle();

// ===== STEP 2B: RUN CELLPOSE SEGMENTATION ON ROI MARKER =====

print("Running Cellpose segmentation on ROI marker channel...");

// Run Cellpose using BIOP wrapper
selectWindow(roiMarkerImage);

run("Cellpose Advanced", 
    "diameter=0 " +
    "cellproba_threshold=0.0 " +
    "flow_threshold=0.4 " +
    "anisotropy=" + d2s(anisotropy, 2) + " " +
    "diam_threshold=12.0 " +
    "model=" + modelPath + " " +
    "nuclei_channel=0 " +
    "cyto_channel=0 " +
    "dimensionmode=3D " +
    "stitch_threshold=-1.0 " +
    "omni=false " +
    "cluster=false " +
    "additional_flags=");

// Cellpose creates a label image
cellposeMask = getTitle();
rename(imageName + "_roiMarker_mask");
maskImage = getTitle();

// Save the mask
saveAs("Tiff", tempDir + imageName + "_roiMarker_cp_masks.tif");
maskImage = getTitle();
print("Cellpose segmentation complete");

// Close roi marker image
selectWindow(roiMarkerImage);
close();

// ===== STEP 3: EXTRACT CHANNEL (TdTomato) =====

if (hasTdTomato && channels >= tdTomatoChannel) {
    print("Extracting TdTomato channel (channel " + tdTomatoChannel + ") - full z-stack...");
    selectWindow(originalImage);
    run("Duplicate...", "duplicate channels=" + tdTomatoChannel);
    tdTomatoImage = getTitle();
    rename(imageName + "_tdTomato");
    print("TdTomato channel extracted");
} else {
    if (!hasTdTomato) {
        print("TdTomato channel not configured, skipping TdTomato analysis");
    } else {
        print("Warning: TdTomato channel " + tdTomatoChannel + " not found in image, skipping TdTomato analysis");
    }
    tdTomatoImage = "";
}

// ===== STEP 3B: EXTRACT CHANNEL (GFAP) FOR SHOLL ANALYSIS =====

if (hasGFAP && channels >= gfapChannel) {
    print("Extracting GFAP channel (channel " + gfapChannel + ") - full z-stack...");
    selectWindow(originalImage);
    run("Duplicate...", "duplicate channels=" + gfapChannel);
    gfapImage = getTitle();
    rename(imageName + "_GFAP");
    
    // Save GFAP channel for later Sholl analysis
    gfapPath = tempDir + imageName + "_GFAP.tif";
    saveAs("Tiff", gfapPath);
    print("Saved GFAP image: " + gfapPath);
    close();
} else {
    if (!hasGFAP) {
        print("GFAP channel not configured, skipping Sholl analysis");
    } else {
        print("Warning: GFAP channel " + gfapChannel + " not found in image, skipping Sholl analysis");
    }
}

// Close original multi-channel image
selectWindow(originalImage);
close();

// ===== STEP 5: PROCESS CELLPOSE MASK AND QUANTIFY =====

print("\n=== Processing Cellpose Mask ===\n");

// Load cellpose mask
maskPath = tempDir + imageName + "_roiMarker_cp_masks.tif";
open(maskPath);
maskImage = getTitle();
print("Loaded cellpose 3D mask");

// Get ROIs from cellpose mask
roiManager("reset");
run("Clear Results");

// Convert 3D mask to ROIs
setOption("BlackBackground", false);
run("Analyze Particles...", "size=10-Infinity stack add");

roiMarkerCount = roiManager("count");
print("Found " + roiMarkerCount + " ROI marker+ cells");

// Save ROI marker ROIs
if (roiMarkerCount > 0) {
    roiManager("save", roiDir + imageName + "_roiMarker_ROIs.zip");
}

// Store mask for boundary checking
selectWindow(maskImage);
run("Duplicate...", "duplicate");
boundaryMask = getTitle();
rename(imageName + "_boundary_check");

// Get centroids and boundary info for each ROI
run("Set Measurements...", "centroid bounding stack redirect=None decimal=3");
roiManager("Measure");

close(maskImage);
close(boundaryMask);

// ===== STEP 6: MEASURE OVERLAP WITH TdTomato =====

if (hasTdTomato) {
    // Re-open the original z-stack
    run("Bio-Formats Importer", "open=[" + lifPath + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
    fullTitle = getTitle();
    getDimensions(width, height, channels, slices, frames);
    
    if (channels >= tdTomatoChannel) {
        run("Duplicate...", "duplicate channels=" + tdTomatoChannel);
        tdTomatoImage = getTitle();
        rename(imageName + "_tdTomato_3D");
        
        // Close original
        selectWindow(fullTitle);
        close();
    }
    
    selectWindow(imageName + "_tdTomato_3D");
    
    // Threshold TdTomato channel (3D stack)
    run("8-bit");
    setAutoThreshold("Triangle dark stack");
    run("Convert to Mask", "method=Triangle background=Dark black");
    
    // Clean up the mask
    for (s = 1; s <= nSlices; s++) {
        setSlice(s);
        run("Despeckle", "slice");
        run("Fill Holes", "slice");
    }
    
    // Save TdTomato mask
    saveAs("Tiff", outputDir + imageName + "_tdTomato_mask_3D.tif");
    tdTomatoMask = getTitle();
    
    print("Created TdTomato 3D mask");
    
    // Measure overlap in 3D and collect cell info
    if (roiMarkerCount > 0) {
        overlapCount = 0;
        
        // Arrays to store cell status
        cellIsTdTomato = newArray(roiMarkerCount);
        cellInBounds = newArray(roiMarkerCount);
        cellCentroidX = newArray(roiMarkerCount);
        cellCentroidY = newArray(roiMarkerCount);
        cellCentroidZ = newArray(roiMarkerCount);
        
        getDimensions(imgWidth, imgHeight, imgChannels, imgSlices, imgFrames);
        
        for (j = 0; j < roiMarkerCount; j++) {
            roiManager("select", j);
            
            cellCentroidX[j] = getResult("X", j);
            cellCentroidY[j] = getResult("Y", j);
            cellCentroidZ[j] = getResult("Slice", j);
            
            Roi.getBounds(roiX, roiY, roiWidth, roiHeight);
            
            // Check XY bounds
            margin = 5;
            if (roiX < margin || roiY < margin || 
                (roiX + roiWidth) > (imgWidth - margin) || 
                (roiY + roiHeight) > (imgHeight - margin)) {
                cellInBounds[j] = 0;
            } else {
                cellInBounds[j] = 1;
            }
            
            // Check Z bounds
            Stack.getPosition(channel, slice, frame);
            minSlice = imgSlices;
            maxSlice = 1;
            
            for (s = 1; s <= imgSlices; s++) {
                Stack.setSlice(s);
                if (Roi.contains(roiX + roiWidth/2, roiY + roiHeight/2)) {
                    if (s < minSlice) minSlice = s;
                    if (s > maxSlice) maxSlice = s;
                }
            }
            
            if (minSlice <= 2 || maxSlice >= (imgSlices - 1)) {
                cellInBounds[j] = 0;
            }
            
            cellCentroidZ[j] = floor((minSlice + maxSlice) / 2);
            
            // Measure TdTomato overlap
            run("Clear Results");
            run("Measure Stack");
            
            hasOverlap = false;
            for (k = 0; k < nResults; k++) {
                meanVal = getResult("Mean", k);
                if (meanVal > 0) {
                    hasOverlap = true;
                    break;
                }
            }
            
            if (hasOverlap) {
                overlapCount++;
                cellIsTdTomato[j] = 1;
            } else {
                cellIsTdTomato[j] = 0;
            }
        }
        
        overlapPercent = (overlapCount / roiMarkerCount) * 100;
        
        print("Overlap: " + overlapCount + "/" + roiMarkerCount + " ROI marker+ cells (" + d2s(overlapPercent, 1) + "%)");
        
        // Count cells fully in bounds
        inBoundsCount = 0;
        for (j = 0; j < roiMarkerCount; j++) {
            if (cellInBounds[j] == 1) inBoundsCount++;
        }
        print("Cells fully in bounds: " + inBoundsCount + "/" + roiMarkerCount);
    }
    
    close(tdTomatoMask);
    
    // ===== STEP 7: SHOLL ANALYSIS USING SNT =====
    
    if (roiMarkerCount > 0 && inBoundsCount > 0 && hasGFAP) {
        print("\n--- Starting Sholl Analysis ---");
        
        gfapPath = tempDir + imageName + "_GFAP.tif";
        
        if (File.exists(gfapPath)) {
            open(gfapPath);
            gfapStack = getTitle();
            print("Loaded GFAP channel for Sholl analysis");
            
            // Prepare GFAP for sensitive tracing
            run("Duplicate...", "duplicate");
            gfapProcessed = getTitle();
            
            print("Enhancing GFAP signal for fine process detection...");
            
            run("Tubeness", "sigma=1.5 use");
            run("Subtract Background...", "rolling=15 stack");
            
            processedGfapPath = tempDir + imageName + "_GFAP_processed.tif";
            saveAs("Tiff", processedGfapPath);
            gfapForSholl = getTitle();
            
            tdTomatoPosCount = 0;
            tdTomatoNegCount = 0;
            
            // Perform Sholl analysis on each cell that's fully in bounds
            for (j = 0; j < roiMarkerCount; j++) {
                if (cellInBounds[j] == 1) {
                    
                    cellGroup = "";
                    if (cellIsTdTomato[j] == 1) {
                        cellGroup = "TdTomato_positive";
                        tdTomatoPosCount++;
                    } else {
                        cellGroup = "TdTomato_negative";
                        tdTomatoNegCount++;
                    }
                    
                    print("  Analyzing cell " + (j+1) + " - " + cellGroup);
                    print("    Centroid: X=" + d2s(cellCentroidX[j], 1) + ", Y=" + d2s(cellCentroidY[j], 1) + ", Z=" + cellCentroidZ[j]);
                    
                    selectWindow(gfapForSholl);
                    makePoint(cellCentroidX[j], cellCentroidY[j]);
                    Stack.setSlice(cellCentroidZ[j]);
                    
                    print("    Running automated tracing...");
                    run("SNT: Start Tracing...", "");
                    wait(500);
                    
                    run("SNT: Set Tracing Parameters", 
                        "xy_resolution=" + d2s(voxelWidth, 4) + " " +
                        "z_resolution=" + d2s(voxelDepth, 4) + " " +
                        "use_preprocessed_image");
                    
                    run("SNT: Automated Tracing", 
                        "use_seeds max_path_length=150");
                    
                    wait(1000);
                    
                    tracesFile = shollDir + imageName + "_cell" + (j+1) + "_" + cellGroup + ".traces";
                    run("SNT: Save Traces", "file=[" + tracesFile + "]");
                    print("    Saved reconstruction: " + tracesFile);
                    
                    swcFile = shollDir + imageName + "_cell" + (j+1) + "_" + cellGroup + ".swc";
                    run("SNT: Export SWC", "file=[" + swcFile + "]");
                    
                    print("    Running Sholl analysis on reconstruction...");
                    run("SNT: Sholl Analysis", 
                        "starting=5 ending=150 radius_step=0 " +
                        "infer integration=[Mean value] " +
                        "polynomial=[Best fitting degree] " +
                        "most normalizer=Volume " +
                        "show save do");
                    
                    if (isOpen("Sholl Results")) {
                        selectWindow("Sholl Results");
                        shollCellFile = shollDir + imageName + "_cell" + (j+1) + "_" + cellGroup + "_sholl.csv";
                        saveAs("Results", shollCellFile);
                        run("Close");
                    }
                    
                    if (isOpen("Sholl Plot")) {
                        selectWindow("Sholl Plot");
                        close();
                    }
                    
                    run("SNT: Clear All Paths");
                    run("SNT: Quit");
                }
            }
            
            print("Sholl analysis complete:");
            print("  TdTomato+ cells analyzed: " + tdTomatoPosCount);
            print("  TdTomato- cells analyzed: " + tdTomatoNegCount);
            
            close(gfapStack);
            close(gfapForSholl);
        }
    }
    
} else {
    // No TdTomato - still do bounds checking and Sholl
    print("TdTomato channel not configured, skipping overlap analysis");
    
    if (roiMarkerCount > 0) {
        cellIsTdTomato = newArray(roiMarkerCount);
        cellInBounds = newArray(roiMarkerCount);
        cellCentroidX = newArray(roiMarkerCount);
        cellCentroidY = newArray(roiMarkerCount);
        cellCentroidZ = newArray(roiMarkerCount);
        
        run("Bio-Formats Importer", "open=[" + lifPath + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
        fullTitle = getTitle();
        getDimensions(imgWidth, imgHeight, imgChannels, imgSlices, imgFrames);
        
        for (j = 0; j < roiMarkerCount; j++) {
            cellIsTdTomato[j] = 0;
            
            roiManager("select", j);
            cellCentroidX[j] = getResult("X", j);
            cellCentroidY[j] = getResult("Y", j);
            cellCentroidZ[j] = getResult("Slice", j);
            
            Roi.getBounds(roiX, roiY, roiWidth, roiHeight);
            margin = 5;
            
            if (roiX < margin || roiY < margin || 
                (roiX + roiWidth) > (imgWidth - margin) || 
                (roiY + roiHeight) > (imgHeight - margin)) {
                cellInBounds[j] = 0;
            } else {
                cellInBounds[j] = 1;
            }
            
            minSlice = imgSlices;
            maxSlice = 1;
            for (s = 1; s <= imgSlices; s++) {
                Stack.setSlice(s);
                if (Roi.contains(roiX + roiWidth/2, roiY + roiHeight/2)) {
                    if (s < minSlice) minSlice = s;
                    if (s > maxSlice) maxSlice = s;
                }
            }
            
            if (minSlice <= 2 || maxSlice >= (imgSlices - 1)) {
                cellInBounds[j] = 0;
            }
            
            cellCentroidZ[j] = floor((minSlice + maxSlice) / 2);
        }
        
        close(fullTitle);
        
        inBoundsCount = 0;
        for (j = 0; j < roiMarkerCount; j++) {
            if (cellInBounds[j] == 1) inBoundsCount++;
        }
        print("Cells fully in bounds: " + inBoundsCount + "/" + roiMarkerCount);
        
        overlapCount = 0;
    }
}

// ===== STEP 8: SAVE RESULTS =====

print("\n\n=== Saving Results ===");

run("Clear Results");
setResult("Image", 0, imageName);
setResult("ROI_Marker_Cells", 0, roiMarkerCount);
if (hasTdTomato && roiMarkerCount > 0) {
    setResult("TdTomato_Positive", 0, overlapCount);
    setResult("Overlap_Percent", 0, d2s(overlapPercent, 2));
} else {
    setResult("TdTomato_Positive", 0, 0);
    setResult("Overlap_Percent", 0, "N/A");
}
setResult("Analysis_Date", 0, analysisDate);
setResult("Analysis_Time", 0, analysisTime);
setResult("Cellpose_Model", 0, modelName);
updateResults();

saveAs("Results", outputDir + imageName + "_results.csv");

print("Results saved to: " + outputDir + imageName + "_results.csv");
print("\n=== ANALYSIS COMPLETE ===");
print("ROIs saved to: " + roiDir);
if (hasGFAP) {
    print("Sholl analysis results saved to: " + shollDir);
}

// Summary
print("\n--- SUMMARY ---");
if (hasTdTomato && roiMarkerCount > 0) {
    print(imageName + ": " + roiMarkerCount + " ROI marker+ cells, " + overlapCount + " TdTomato+ (" + d2s(overlapPercent, 1) + "%)");
} else {
    print(imageName + ": " + roiMarkerCount + " ROI marker+ cells");
}

close("*");
roiManager("reset");

showMessage("Analysis Complete!", 
    "Processed 1 image: " + imageName + "\n \n" +
    "Results saved to:\n" + outputDir + "\n \n" +
    "Check the Log window for detailed results.");
