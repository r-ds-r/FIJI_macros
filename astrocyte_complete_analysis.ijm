// Macro: Multi-channel Z-stack Analysis with Cellpose-SAM Segmentation and Sholl Analysis
// Process .lif files with:
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
outputDir = lifDir + File.separator + lifBaseName + "_output" + File.separator;
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

print("\n=== Multi-channel Z-stack Analysis with Cellpose-SAM ===");
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

// ===== STEP 1: OPEN AND PROCESS IMAGES =====

run("Bio-Formats Macro Extensions");
Ext.setId(lifPath);
Ext.getSeriesCount(seriesCount);

print("Found " + seriesCount + " images in .lif file\n");

// Initialize results arrays
imageNames = newArray(seriesCount);
sox9Counts = newArray(seriesCount);
tdTomatoCounts = newArray(seriesCount);
overlapCounts = newArray(seriesCount);
overlapPercentages = newArray(seriesCount);

// Process each series
for (i = 0; i < seriesCount; i++) {
    
    print("\n--- Processing Image " + (i+1) + "/" + seriesCount + " ---");
    
    // Open the current series
    run("Bio-Formats Importer", "open=[" + lifPath + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (i+1));
    
    fullTitle = getTitle();
    print("Opened: " + fullTitle);
    
    // Extract image name
    if (indexOf(fullTitle, " - ") > 0) {
        imageName = substring(fullTitle, indexOf(fullTitle, " - ") + 3);
    } else {
        imageName = fullTitle;
    }
    imageName = replace(imageName, "/", "_"); // Remove problematic characters
    imageNames[i] = imageName;
    
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
    // Make sure BIOP Cellpose is installed: Help > Update > Manage Update Sites > PTBIOP
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
    cellposeMask = getTitle(); // Should be something like "Label Image"
    rename(imageName + "_roiMarker_mask");
    maskImage = getTitle();
    
    // Save the mask
    saveAs("Tiff", tempDir + imageName + "_roiMarker_cp_masks.tif");
    maskImage = getTitle();
    print("Cellpose segmentation complete");
    
    // Close roi marker image
    selectWindow(roiMarkerImage);
    close();
    
    // Keep this image open for later
    
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
    
    // Close all images for now, we'll reopen after cellpose
    close("*");
    
}

print("\n=== All images processed and segmented ===");

// ===== STEP 5: PROCESS CELLPOSE MASKS AND QUANTIFY OVERLAP =====

print("\n=== Processing Cellpose Masks ===\n");

for (i = 0; i < seriesCount; i++) {
    
    imageName = imageNames[i];
    print("\n--- Analyzing " + imageName + " ---");
    
    // Load cellpose mask (3D) that was already generated
    maskPath = tempDir + imageName + "_roiMarker_cp_masks.tif";
    
    if (!File.exists(maskPath)) {
        print("ERROR: Cellpose mask not found: " + maskPath);
        print("Skipping this image");
        sox9Counts[i] = 0;
        tdTomatoCounts[i] = 0;
        overlapCounts[i] = 0;
        overlapPercentages[i] = 0;
        continue;
    }
    
    open(maskPath);
    maskImage = getTitle();
    print("Loaded cellpose 3D mask");
    
    // Get ROIs from cellpose mask
    roiManager("reset");
    run("Clear Results");
    
    // Convert 3D mask to ROIs
    // Cellpose masks have unique integer value for each cell
    // For 3D, we need to use 3D Objects Counter or stack-based analysis
    setOption("BlackBackground", false);
    run("Analyze Particles...", "size=10-Infinity stack add");
    
    roiMarkerCount = roiManager("count");
    sox9Counts[i] = roiMarkerCount;
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
        // Re-open the original z-stack to get TdTomato channel
        run("Bio-Formats Importer", "open=[" + lifPath + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (i+1));
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
        
        // Clean up the mask (applied to stack)
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
                
                // Get centroid and bounds from Results table (measured earlier)
                cellCentroidX[j] = getResult("X", j);
                cellCentroidY[j] = getResult("Y", j);
                cellCentroidZ[j] = getResult("Slice", j); // Will use middle slice of ROI
                
                // Get bounding box to check if cell is fully in bounds
                Roi.getBounds(roiX, roiY, roiWidth, roiHeight);
                
                // Check XY bounds (allow small margin from edge)
                margin = 5; // pixels
                if (roiX < margin || roiY < margin || 
                    (roiX + roiWidth) > (imgWidth - margin) || 
                    (roiY + roiHeight) > (imgHeight - margin)) {
                    cellInBounds[j] = 0; // Touches edge
                } else {
                    cellInBounds[j] = 1; // Fully in bounds
                }
                
                // Check Z bounds - need to check across all slices where ROI exists
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
                    cellInBounds[j] = 0; // Touches Z boundary
                }
                
                cellCentroidZ[j] = floor((minSlice + maxSlice) / 2); // Center Z
                
                // Measure TdTomato overlap across all slices
                run("Clear Results");
                run("Measure Stack");
                
                // Check if any slice has signal
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
            
            overlapCounts[i] = overlapCount;
            overlapPercent = (overlapCount / roiMarkerCount) * 100;
            overlapPercentages[i] = overlapPercent;
            
            print("Overlap: " + overlapCount + "/" + roiMarkerCount + " ROI marker+ cells (" + d2s(overlapPercent, 1) + "%)");
            
            // Count cells fully in bounds
            inBoundsCount = 0;
            for (j = 0; j < roiMarkerCount; j++) {
                if (cellInBounds[j] == 1) inBoundsCount++;
            }
            print("Cells fully in bounds: " + inBoundsCount + "/" + roiMarkerCount);
            
        } else {
            overlapCounts[i] = 0;
            overlapPercentages[i] = 0;
        }
        
        close(tdTomatoMask);
        
        // ===== STEP 7: SHOLL ANALYSIS USING SNT =====
        
        if (roiMarkerCount > 0 && inBoundsCount > 0 && hasGFAP) {
            print("\n--- Starting Sholl Analysis ---");
            
            // Load GFAP channel for Sholl analysis
            gfapPath = tempDir + imageName + "_GFAP.tif";
            
            if (File.exists(gfapPath)) {
                open(gfapPath);
                gfapStack = getTitle();
                print("Loaded GFAP channel for Sholl analysis");
                
                // Prepare GFAP for sensitive tracing (preserve fine processes)
                // Use Tubeness filter to enhance processes
                run("Duplicate...", "duplicate");
                gfapProcessed = getTitle();
                
                print("Enhancing GFAP signal for fine process detection...");
                
                // Apply Tubeness filter to enhance filamentous structures
                run("Tubeness", "sigma=1.5 use");
                
                // Gentle background subtraction
                run("Subtract Background...", "rolling=15 stack");
                
                // Save processed GFAP
                processedGfapPath = tempDir + imageName + "_GFAP_processed.tif";
                saveAs("Tiff", processedGfapPath);
                gfapForSholl = getTitle();
                
                // Initialize Sholl results storage
                shollResults_TdTomato_pos = "";
                shollResults_TdTomato_neg = "";
                
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
                        
                        // Run SNT Auto-Tracing first to create reconstructions
                        selectWindow(gfapForSholl);
                        
                        // Create a point ROI at the centroid for tracing origin
                        makePoint(cellCentroidX[j], cellCentroidY[j]);
                        Stack.setSlice(cellCentroidZ[j]);
                        
                        // Run SNT automated tracing
                        // This creates actual path reconstructions that can be saved as .traces
                        print("    Running automated tracing...");
                        run("SNT: Start Tracing...", "");
                        wait(500); // Give SNT time to initialize
                        
                        // Set tracing parameters for fine processes
                        run("SNT: Set Tracing Parameters", 
                            "xy_resolution=" + d2s(voxelWidth, 4) + " " +
                            "z_resolution=" + d2s(voxelDepth, 4) + " " +
                            "use_preprocessed_image");
                        
                        // Auto-trace from the centroid point
                        run("SNT: Automated Tracing", 
                            "use_seeds max_path_length=150");
                        
                        wait(1000); // Give time for tracing to complete
                        
                        // Save the traced paths as .traces file
                        tracesFile = shollDir + imageName + "_cell" + (j+1) + "_" + cellGroup + ".traces";
                        run("SNT: Save Traces", "file=[" + tracesFile + "]");
                        print("    Saved reconstruction: " + tracesFile);
                        
                        // Also save as SWC format for compatibility
                        swcFile = shollDir + imageName + "_cell" + (j+1) + "_" + cellGroup + ".swc";
                        run("SNT: Export SWC", "file=[" + swcFile + "]");
                        
                        // Now run Sholl Analysis on the traced paths
                        print("    Running Sholl analysis on reconstruction...");
                        run("SNT: Sholl Analysis", 
                            "starting=5 ending=150 radius_step=0 " +
                            "infer integration=[Mean value] " +
                            "polynomial=[Best fitting degree] " +
                            "most normalizer=Volume " +
                            "show save do");
                        
                        // Get the Sholl results table
                        if (isOpen("Sholl Results")) {
                            selectWindow("Sholl Results");
                            
                            // Save individual cell Sholl results
                            shollCellFile = shollDir + imageName + "_cell" + (j+1) + "_" + cellGroup + "_sholl.csv";
                            saveAs("Results", shollCellFile);
                            
                            // Also append to group files
                            if (cellIsTdTomato[j] == 1) {
                                // Add to TdTomato positive group
                                if (shollResults_TdTomato_pos == "") {
                                    shollResults_TdTomato_pos = shollDir + imageName + "_TdTomato_positive_all.csv";
                                    saveAs("Results", shollResults_TdTomato_pos);
                                } else {
                                    // Append to existing file
                                    File.append("", shollResults_TdTomato_pos); // Add separator
                                    File.append("Cell " + (j+1), shollResults_TdTomato_pos);
                                }
                            } else {
                                // Add to TdTomato negative group
                                if (shollResults_TdTomato_neg == "") {
                                    shollResults_TdTomato_neg = shollDir + imageName + "_TdTomato_negative_all.csv";
                                    saveAs("Results", shollResults_TdTomato_neg);
                                } else {
                                    File.append("", shollResults_TdTomato_neg);
                                    File.append("Cell " + (j+1), shollResults_TdTomato_neg);
                                }
                            }
                            
                            run("Close");
                        }
                        
                        // Close any Sholl plots
                        if (isOpen("Sholl Plot")) {
                            selectWindow("Sholl Plot");
                            close();
                        }
                        
                        // Clear SNT for next cell
                        run("SNT: Clear All Paths");
                        run("SNT: Quit");
                    }
                }
                
                print("Sholl analysis complete:");
                print("  TdTomato+ cells analyzed: " + tdTomatoPosCount);
                print("  TdTomato- cells analyzed: " + tdTomatoNegCount);
                
                close(gfapStack);
                close(gfapForSholl);
                
            } else {
                print("Warning: GFAP channel not found, skipping Sholl analysis");
            }
        } else {
            if (!hasGFAP) {
                print("GFAP channel not configured, skipping Sholl analysis");
            }
        }
        
    } else {
        print("TdTomato channel not configured, skipping overlap analysis");
        
        // Still need to initialize arrays for cells without TdTomato
        if (roiMarkerCount > 0) {
            cellIsTdTomato = newArray(roiMarkerCount);
            cellInBounds = newArray(roiMarkerCount);
            cellCentroidX = newArray(roiMarkerCount);
            cellCentroidY = newArray(roiMarkerCount);
            cellCentroidZ = newArray(roiMarkerCount);
            
            // All cells are TdTomato negative and need bounds checking
            // Load image to get dimensions
            run("Bio-Formats Importer", "open=[" + lifPath + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (i+1));
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
                
                // Check Z bounds
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
            
            // Run Sholl analysis if GFAP is available
            if (hasGFAP && roiMarkerCount > 0 && inBoundsCount > 0) {
                print("\n--- Starting Sholl Analysis (No TdTomato grouping) ---");
                
                gfapPath = tempDir + imageName + "_GFAP.tif";
                
                if (File.exists(gfapPath)) {
                    open(gfapPath);
                    gfapStack = getTitle();
                    print("Loaded GFAP channel for Sholl analysis");
                    
                    // Prepare GFAP for sensitive tracing (preserve fine processes)
                    run("Duplicate...", "duplicate");
                    gfapProcessed = getTitle();
                    
                    print("Enhancing GFAP signal for fine process detection...");
                    
                    // Apply Tubeness filter to enhance filamentous structures
                    run("Tubeness", "sigma=1.5 use");
                    
                    // Gentle background subtraction
                    run("Subtract Background...", "rolling=15 stack");
                    
                    // Save processed GFAP
                    processedGfapPath = tempDir + imageName + "_GFAP_processed.tif";
                    saveAs("Tiff", processedGfapPath);
                    gfapForSholl = getTitle();
                    
                    cellsAnalyzed = 0;
                    
                    // Perform Sholl analysis on each cell that's fully in bounds
                    for (j = 0; j < roiMarkerCount; j++) {
                        if (cellInBounds[j] == 1) {
                            
                            cellsAnalyzed++;
                            
                            print("  Analyzing cell " + (j+1));
                            print("    Centroid: X=" + d2s(cellCentroidX[j], 1) + ", Y=" + d2s(cellCentroidY[j], 1) + ", Z=" + cellCentroidZ[j]);
                            
                            // Run SNT Auto-Tracing first to create reconstructions
                            selectWindow(gfapForSholl);
                            
                            // Create a point ROI at the centroid for tracing origin
                            makePoint(cellCentroidX[j], cellCentroidY[j]);
                            Stack.setSlice(cellCentroidZ[j]);
                            
                            // Run SNT automated tracing
                            print("    Running automated tracing...");
                            run("SNT: Start Tracing...", "");
                            wait(500);
                            
                            // Set tracing parameters for fine processes
                            run("SNT: Set Tracing Parameters", 
                                "xy_resolution=" + d2s(voxelWidth, 4) + " " +
                                "z_resolution=" + d2s(voxelDepth, 4) + " " +
                                "use_preprocessed_image");
                            
                            // Auto-trace from the centroid point
                            run("SNT: Automated Tracing", 
                                "use_seeds max_path_length=150");
                            
                            wait(1000);
                            
                            // Save the traced paths as .traces file
                            tracesFile = shollDir + imageName + "_cell" + (j+1) + ".traces";
                            run("SNT: Save Traces", "file=[" + tracesFile + "]");
                            print("    Saved reconstruction: " + tracesFile);
                            
                            // Also save as SWC format for compatibility
                            swcFile = shollDir + imageName + "_cell" + (j+1) + ".swc";
                            run("SNT: Export SWC", "file=[" + swcFile + "]");
                            
                            // Now run Sholl Analysis on the traced paths
                            print("    Running Sholl analysis on reconstruction...");
                            run("SNT: Sholl Analysis", 
                                "starting=5 ending=150 radius_step=0 " +
                                "infer integration=[Mean value] " +
                                "polynomial=[Best fitting degree] " +
                                "most normalizer=Volume " +
                                "show save do");
                            
                            // Get the Sholl results table
                            if (isOpen("Sholl Results")) {
                                selectWindow("Sholl Results");
                                
                                // Save individual cell Sholl results
                                shollCellFile = shollDir + imageName + "_cell" + (j+1) + "_sholl.csv";
                                saveAs("Results", shollCellFile);
                                
                                run("Close");
                            }
                            
                            // Close any Sholl plots
                            if (isOpen("Sholl Plot")) {
                                selectWindow("Sholl Plot");
                                close();
                            }
                            
                            // Clear SNT for next cell
                            run("SNT: Clear All Paths");
                            run("SNT: Quit");
                        }
                    }
                    
                    print("Sholl analysis complete: " + cellsAnalyzed + " cells analyzed");
                    
                    close(gfapStack);
                    close(gfapForSholl);
                    
                } else {
                    print("Warning: GFAP channel not found, skipping Sholl analysis");
                }
            }
        }
        
        overlapCounts[i] = 0;
        overlapPercentages[i] = 0;
        close("*");
    }
    
    close("*");
    roiManager("reset");
}

// ===== STEP 8: SAVE RESULTS =====

print("\n\n=== Saving Results ===");

run("Clear Results");
for (i = 0; i < seriesCount; i++) {
    setResult("Image", i, imageNames[i]);
    setResult("ROI_Marker_Cells", i, sox9Counts[i]);
    setResult("TdTomato_Positive", i, overlapCounts[i]);
    setResult("Overlap_Percent", i, d2s(overlapPercentages[i], 2));
    setResult("Analysis_Date", i, analysisDate);
    setResult("Analysis_Time", i, analysisTime);
    setResult("Cellpose_Model", i, modelName);
}
updateResults();

// Save results
saveAs("Results", outputDir + lifBaseName + "_roiMarker_tdTomato_results.csv");

print("Results saved to: " + outputDir + lifBaseName + "_roiMarker_tdTomato_results.csv");
print("\n=== ANALYSIS COMPLETE ===");
print("Total images processed: " + seriesCount);
print("ROIs saved to: " + roiDir);
if (hasGFAP) {
    print("Sholl analysis results saved to: " + shollDir);
}

// Summary
print("\n--- SUMMARY ---");
for (i = 0; i < seriesCount; i++) {
    if (hasTdTomato) {
        print(imageNames[i] + ": " + sox9Counts[i] + " ROI marker+ cells, " + overlapCounts[i] + " TdTomato+ (" + d2s(overlapPercentages[i], 1) + "%)");
    } else {
        print(imageNames[i] + ": " + sox9Counts[i] + " ROI marker+ cells");
    }
}

showMessage("Analysis Complete!", 
    "Processed " + seriesCount + " images.\n \n" +
    "Results saved to:\n" + outputDir + "\n \n" +
    "- Cell counts and TdTomato overlap\n" +
    "- Sholl analysis by TdTomato status\n \n" +
    "Check the Log window for detailed results.");
