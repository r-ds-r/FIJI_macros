
// Batch processing macro for .lif files
// Analyzes cells from Ch2 and checks for Ch3 colocalization
// Two-phase processing: 1) Collect thresholds, 2) Analyze with stored thresholds

// Function to parse filename and extract experimental conditions
function parseConditions(imageName) {
    // Expected format: [day]_[genotype][sex]-[#]_[marker]_10x
    // Example: d7_GNM-1_s100b_10x
    // But Bio-Formats adds project name prefix: "project.lif - d7_GNM-1_s100b_10x"
    
    conditionsResult = newArray(5); // day, genotype, sex, replicate, marker
    
    // Remove the .lif project name prefix if present
    cleanName = imageName;
    if (indexOf(imageName, " - ") > 0) {
        splitParts = split(imageName, " - ");
        cleanName = splitParts[1]; // Take the part after " - "
    }
    
    // Split by underscore
    nameParts = split(cleanName, "_");
    
    if (nameParts.length >= 3) {
        // Extract day (e.g., "d7")
        conditionsResult[0] = nameParts[0];
        
        // Extract genotype, sex, and replicate from second part (e.g., "GNM-1")
        genoSexRep = nameParts[1];
        
        // Find the dash to separate genotype+sex from replicate
        dashIdx = indexOf(genoSexRep, "-");
        if (dashIdx > 0) {
            genoSex = substring(genoSexRep, 0, dashIdx);
            conditionsResult[3] = substring(genoSexRep, dashIdx + 1); // replicate
            
            // Last character before dash is sex (M or F)
            sexIdx = lengthOf(genoSex) - 1;
            conditionsResult[2] = substring(genoSex, sexIdx, sexIdx + 1); // sex
            
            // Everything before sex is genotype
            conditionsResult[1] = substring(genoSex, 0, sexIdx); // genotype
        }
        
        // Extract marker (e.g., "s100b")
        if (nameParts.length >= 3) {
            conditionsResult[4] = nameParts[2];
        }
    } else {
        // If parsing fails, fill with "unknown"
        for (ii = 0; ii < 5; ii++) {
            conditionsResult[ii] = "unknown";
        }
    }
    
    return conditionsResult;
}

// Close any open images and clear ROI manager
close("*");
run("Clear Results");
roiManager("reset");

// Ask user if they have pre-existing thresholds
Dialog.create("Threshold Mode");
Dialog.addMessage("Do you have pre-existing threshold values from a previous run?");
Dialog.addChoice("Mode:", newArray("Collect new thresholds (Phase 1 + 2)", "Use existing thresholds (Phase 2 only)"));
Dialog.show();
thresholdMode = Dialog.getChoice();

useExistingThresholds = (thresholdMode == "Use existing thresholds (Phase 2 only)");

// Prompt user to select the .lif file
lifFile = File.openDialog("Select .lif file");
if (lifFile == "") {
    exit("No file selected");
}

// Get directory to save results
dir = File.getDirectory(lifFile);
resultsDir = dir + "Results" + File.separator;
File.makeDirectory(resultsDir);

// Get the number of series (images) in the .lif file
run("Bio-Formats Macro Extensions");
Ext.setId(lifFile);
Ext.getSeriesCount(seriesCount);

print("\\Clear");
print("Processing " + seriesCount + " images from:");
print(lifFile);
print("---");

// Initialize threshold storage arrays
sliceInterval = 5; // Threshold interval
maxThresholds = 100; // Maximum number of threshold measurements
ch2Thresholds = newArray(seriesCount * maxThresholds);
ch3Thresholds = newArray(seriesCount * maxThresholds);
thresholdSlices = newArray(seriesCount * maxThresholds);
seriesSliceCounts = newArray(seriesCount);

// Create randomized order for blinding
seriesOrder = newArray(seriesCount);
for (i = 0; i < seriesCount; i++) {
    seriesOrder[i] = i;
}

if (useExistingThresholds) {
    // Load existing thresholds from CSV
    print("===");
    print("LOADING EXISTING THRESHOLDS");
    print("---");
    
    thresholdFile = File.openDialog("Select Threshold_Values.csv file");
    if (thresholdFile == "") {
        exit("No threshold file selected");
    }
    
    print("Loading thresholds from: " + thresholdFile);
    
    // Read the CSV file
    fileString = File.openAsString(thresholdFile);
    lines = split(fileString, "\n");
    
    // Skip header line, process data lines
    for (lineIdx = 1; lineIdx < lines.length; lineIdx++) {
        if (lengthOf(lines[lineIdx]) > 0) {
            columns = split(lines[lineIdx], ",");
            if (columns.length >= 5) {
                seriesIdx = parseInt(columns[0]) - 1; // Convert to 0-based
                sliceNum = parseInt(columns[2]);
                ch2Val = parseFloat(columns[3]);
                ch3Val = parseFloat(columns[4]);
                
                // Calculate storage index
                thresholdIdx = (sliceNum / sliceInterval) - 1;
                storageIdx = seriesIdx * maxThresholds + thresholdIdx;
                
                ch2Thresholds[storageIdx] = ch2Val;
                ch3Thresholds[storageIdx] = ch3Val;
            }
        }
    }
    
    // Get slice counts by opening each series briefly
    for (s = 0; s < seriesCount; s++) {
        run("Bio-Formats Importer", 
            "open=[" + lifFile + "] " +
            "autoscale " +
            "color_mode=Composite " +
            "view=Hyperstack " +
            "stack_order=XYCZT " +
            "series_" + (s+1));
        getDimensions(width, height, channels, tempSlices, frames);
        seriesSliceCounts[s] = tempSlices;
        close();
    }
    
    print("Thresholds loaded successfully.");
    print("Proceeding directly to Phase 2 analysis...");
    print("");
    
} else {
    // Original Phase 1: Collect thresholds

// Shuffle array (Fisher-Yates algorithm)
for (i = seriesCount - 1; i > 0; i--) {
    j = floor(random * (i + 1));
    temp = seriesOrder[i];
    seriesOrder[i] = seriesOrder[j];
    seriesOrder[j] = temp;
}

print("Randomized processing order: " + String.join(seriesOrder, ", "));
print("===");
print("");

//========================================
// PHASE 1: COLLECT THRESHOLDS
//========================================
print("PHASE 1: Threshold Collection");
print("---");

sliceInterval = 5; // Threshold interval
maxThresholds = 100; // Maximum number of threshold measurements

// Storage arrays for thresholds
// Format: thresholdData[seriesIndex][channelIndex][sliceIndex] = threshold value
ch2Thresholds = newArray(seriesCount * maxThresholds);
ch3Thresholds = newArray(seriesCount * maxThresholds);
thresholdSlices = newArray(seriesCount * maxThresholds);
seriesSliceCounts = newArray(seriesCount);

// Loop through each series in random order for threshold collection
for (idx = 0; idx < seriesCount; idx++) {
    s = seriesOrder[idx];
    
    print("Collecting thresholds for series " + (s+1) + " (" + (idx+1) + " of " + seriesCount + ")");
    
    // Open the image with Bio-Formats
    run("Bio-Formats Importer", 
        "open=[" + lifFile + "] " +
        "autoscale " +
        "color_mode=Composite " +
        "view=Hyperstack " +
        "stack_order=XYCZT " +
        "series_" + (s+1));
    
    originalTitle = getTitle();
    seriesName = originalTitle;
    
    // Parse experimental conditions
    conditions = parseConditions(seriesName);
    day = conditions[0];
    genotype = conditions[1];
    sex = conditions[2];
    replicate = conditions[3];
    marker = conditions[4];
    
    getDimensions(width, height, channels, slices, frames);
    seriesSliceCounts[s] = slices;
    
    print("  Processing blinded image " + (idx+1) + " of " + seriesCount);
    print("  " + slices + " slices - collecting thresholds every " + sliceInterval + " slices");
    
    if (channels < 3) {
        print("  WARNING: Image has less than 3 channels. Skipping...");
        close();
        continue;
    }
    
    // Extract Ch2 and Ch3 stacks
    selectWindow(originalTitle);
    run("Duplicate...", "title=ch2_temp duplicate channels=2");
    selectWindow(originalTitle);
    run("Duplicate...", "title=ch3_temp duplicate channels=3");
    close(originalTitle);
    
    // Collect Ch2 thresholds - present each key slice individually
    print("  Collecting Ch2 (S100B) thresholds...");
    selectWindow("ch2_temp");
    
    numThresholdSlices = floor(slices / sliceInterval);
    
    for (i = 0; i < numThresholdSlices; i++) {
        z = (i + 1) * sliceInterval;
        if (z > slices) break;
        
        selectWindow("ch2_temp");
        setSlice(z);
        run("Duplicate...", "title=ch2_slice_" + z);
        
        // Apply auto threshold
        setAutoThreshold("Otsu dark");
        
        // Show threshold dialog for manual adjustment
        run("Threshold...");
        waitForUser("Ch2 Threshold - Slice " + z + " of " + slices, 
            "Adjust threshold for Ch2 at slice " + z + ".\\n" +
            "This threshold will be applied to slices " + ((i * sliceInterval) + 1) + "-" + z + ".\\n\\n" +
            "Adjust using the Threshold window if needed.\\n" +
            "Click OK when satisfied.");
        
        // Get the threshold value
        getThreshold(lower, upper);
        storageIndex = s * maxThresholds + i;
        ch2Thresholds[storageIndex] = lower;
        thresholdSlices[storageIndex] = z;
        
        // Close threshold window and slice
        selectWindow("Threshold");
        run("Close");
        selectWindow("ch2_slice_" + z);
        close();
    }
    
    close("ch2_temp");
    
    // Collect Ch3 thresholds - present each key slice individually
    print("  Collecting Ch3 (GFAP) thresholds...");
    selectWindow("ch3_temp");
    
    for (i = 0; i < numThresholdSlices; i++) {
        z = (i + 1) * sliceInterval;
        if (z > slices) break;
        
        selectWindow("ch3_temp");
        setSlice(z);
        run("Duplicate...", "title=ch3_slice_" + z);
        
        // Apply auto threshold
        setAutoThreshold("Triangle dark");
        
        // Show threshold dialog for manual adjustment
        run("Threshold...");
        waitForUser("Ch3 Threshold - Slice " + z + " of " + slices, 
            "Adjust threshold for Ch3 (GFAP) at slice " + z + ".\\n" +
            "This threshold will be applied to slices " + ((i * sliceInterval) + 1) + "-" + z + ".\\n\\n" +
            "Adjust using the Threshold window if needed.\\n" +
            "Click OK when satisfied.");
        
        // Get the threshold value
        getThreshold(lower, upper);
        storageIndex = s * maxThresholds + i;
        ch3Thresholds[storageIndex] = lower;
        
        // Close threshold window and slice
        selectWindow("Threshold");
        run("Close");
        selectWindow("ch3_slice_" + z);
        close();
    }
    
    close("ch3_temp");
    
    // Force garbage collection to free memory
    run("Collect Garbage");
    
    print("  Thresholds collected.");
    print("---");
}

} // End of Phase 1 (if not using existing thresholds)

print("");
print("===");
print("PHASE 2: Analysis with Stored Thresholds");
print("---");
print("");

// Save threshold values to file for reproducibility (only if new thresholds were collected)
if (!useExistingThresholds) {
    thresholdLogPath = resultsDir + "Threshold_Values.csv";
    thresholdLog = File.open(thresholdLogPath);
    print(thresholdLog, "Series,Image_Name,Slice,Ch2_Threshold,Ch3_Threshold");

    for (idx = 0; idx < seriesCount; idx++) {
        s = seriesOrder[idx];
        
        // Get image name (will be populated during processing, use placeholder for now)
        imageName = "Series_" + (s+1);
        
        // Write threshold values for this series
        numSlices = seriesSliceCounts[s];
        for (z = sliceInterval; z <= numSlices; z += sliceInterval) {
            thresholdIndex = (z / sliceInterval - 1);
            storageIndex = s * maxThresholds + thresholdIndex;
            ch2Val = ch2Thresholds[storageIndex];
            ch3Val = ch3Thresholds[storageIndex];
            print(thresholdLog, (s+1) + "," + imageName + "," + z + "," + ch2Val + "," + ch3Val);
        }
    }
    File.close(thresholdLog);
    print("Threshold values saved to: Threshold_Values.csv");
    print("");
}

// Create summary results table (expanded for condition data)
summaryResults = newArray(seriesCount * 10); // day, genotype, sex, replicate, marker, image, cells, coloc, percentage
row = 0;

// Loop through each series in order (randomized if new thresholds, sequential if existing)
for (idx = 0; idx < seriesCount; idx++) {
    if (useExistingThresholds) {
        s = idx; // Process in sequential order
    } else {
        s = seriesOrder[idx]; // Use randomized order from Phase 1
    }
    
    print("Analyzing series " + (s+1) + " (" + (idx+1) + " of " + seriesCount + ")");
    
    // Open the image with Bio-Formats
    run("Bio-Formats Importer", 
        "open=[" + lifFile + "] " +
        "autoscale " +
        "color_mode=Composite " +
        "view=Hyperstack " +
        "stack_order=XYCZT " +
        "series_" + (s+1));
    
    originalTitle = getTitle();
    seriesName = originalTitle;
    
    // Parse experimental conditions from filename
    conditions = parseConditions(seriesName);
    day = conditions[0];
    genotype = conditions[1];
    sex = conditions[2];
    replicate = conditions[3];
    marker = conditions[4];
    
    // Get image properties
    getDimensions(width, height, channels, slices, frames);
    print("  Image: " + seriesName);
    print("  Conditions: Day=" + day + ", Genotype=" + genotype + ", Sex=" + sex + ", Rep=" + replicate + ", Marker=" + marker);
    print("  Dimensions: " + width + "x" + height + ", " + channels + " channels, " + slices + " slices");
    
    // Check if we have at least 3 channels
    if (channels < 3) {
        print("  WARNING: Image has less than 3 channels. Skipping...");
        close();
        continue;
    }
    
    // Remove Channel 1 by creating a new hyperstack without it
    run("Duplicate...", "title=working duplicate");
    close(originalTitle);
    selectWindow("working");
    
    // Extract Channel 2 (S100B - nuclear marker) for 3D segmentation
    run("Duplicate...", "title=ch2_stack duplicate channels=2");
    
    // Create 3D mask from Channel 2 using stored slice-specific thresholds
    selectWindow("ch2_stack");
    run("Duplicate...", "title=ch2_mask duplicate");
    
    // Apply slice-specific thresholds by processing each slice
    selectWindow("ch2_mask");
    for (z = 1; z <= slices; z++) {
        setSlice(z);
        
        // Determine which threshold to use for this slice
        thresholdIndex = floor((z - 1) / sliceInterval);
        storageIndex = s * maxThresholds + thresholdIndex;
        thresholdValue = ch2Thresholds[storageIndex];
        
        // Apply threshold to this slice using pixel manipulation
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                pixelValue = getPixel(x, y);
                if (pixelValue >= thresholdValue) {
                    setPixel(x, y, 255);
                } else {
                    setPixel(x, y, 0);
                }
            }
        }
    }
    
    // Convert to 8-bit if needed
    selectWindow("ch2_mask");
    run("8-bit");
    
    // Clean up mask in 3D (optional - requires 3D ImageJ Suite plugin)
    // Uncomment if you have the plugin installed:
    // run("3D Fill Holes");
    
    // 3D segmentation using 3D Objects Counter
    // Requires: 3D ImageJ Suite plugin (built into FIJI)
    run("3D Objects Counter", 
        "threshold=128 " +
        "slice=1 " +
        "min.=500 " +  // Adjust minimum nucleus volume (voxels)
        "max.=100000 " +  // Adjust maximum nucleus volume
        "objects surfaces");
    
    // Get cell count from 3D Objects Counter results
    cellCount = 0;
    if (nResults > 0) {
        cellCount = nResults;
    }
    print("  Detected cells (Ch2 3D): " + cellCount);
    
    // Save the Ch2 binary mask for reproducibility
    selectWindow("ch2_mask");
    saveAs("Tiff", resultsDir + seriesName + "_Ch2_mask.tif");
    close(); // Close after saving
    
    // Close unnecessary windows
    if (isOpen("ch2_stack")) {
        selectWindow("ch2_stack");
        close();
    }
    
    // Save the 3D objects map
    selectWindow("Objects map of ch2_mask");
    rename("objects_map");
    
    // Close unneeded windows
    if (isOpen("Surface map of ch2_mask")) {
        selectWindow("Surface map of ch2_mask");
        close();
    }
    
    // Now check 3D colocalization with Channel 3 (GFAP)
    if (cellCount > 0) {
        
        // Extract Channel 3 (GFAP)
        selectWindow("working");
        run("Duplicate...", "title=ch3_stack duplicate channels=3");
        
        // Threshold Channel 3 to detect GFAP fibrils in 3D using stored slice-specific thresholds
        run("Duplicate...", "title=ch3_thresh duplicate");
        
        // Apply slice-specific thresholds by processing each slice
        selectWindow("ch3_thresh");
        for (z = 1; z <= slices; z++) {
            setSlice(z);
            
            // Determine which threshold to use for this slice
            thresholdIndex = floor((z - 1) / sliceInterval);
            storageIndex = s * maxThresholds + thresholdIndex;
            thresholdValue = ch3Thresholds[storageIndex];
            
            // Apply threshold to this slice using pixel manipulation
            for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                    pixelValue = getPixel(x, y);
                    if (pixelValue >= thresholdValue) {
                        setPixel(x, y, 255);
                    } else {
                        setPixel(x, y, 0);
                    }
                }
            }
        }
        
        // Convert to 8-bit if needed
        selectWindow("ch3_thresh");
        run("8-bit");
        
        // Count colocalization by checking overlap between objects and GFAP
        colocCount = 0;
        
        selectWindow("objects_map");
        
        // For each detected cell (3D object)
        for (i = 1; i <= cellCount; i++) {
            
            // Create mask for this specific object
            selectWindow("objects_map");
            run("Duplicate...", "title=single_object duplicate");
            setThreshold(i, i);
            run("Convert to Mask", "method=Default background=Dark black");
            
            // Dilate the object mask slightly to capture peri-nuclear GFAP
            run("Options...", "iterations=2 count=1 black do=Dilate stack");
            
            // Check overlap with GFAP signal
            imageCalculator("AND create stack", "single_object", "ch3_thresh");
            selectWindow("Result of single_object");
            
            // Count non-zero voxels in the overlap
            getRawStatistics(nPixels, mean);
            overlapVoxels = nPixels * mean / 255;
            
            // If there's overlap, count as positive
            if (overlapVoxels > 10) {  // Adjust minimum overlap threshold
                colocCount++;
            }
            
            // Clean up temporary images
            close("single_object");
            close("Result of single_object");
        }
        
        print("  Cells with Ch3 signal (3D): " + colocCount);
        print("  Percentage: " + d2s((colocCount/cellCount)*100, 2) + "%");
        
        // Save the Ch3 binary mask for reproducibility
        selectWindow("ch3_thresh");
        saveAs("Tiff", resultsDir + seriesName + "_Ch3_mask.tif");
        close();
        
        // Close ch3_stack to free memory
        if (isOpen("ch3_stack")) {
            selectWindow("ch3_stack");
            close();
        }
        
        // Save the objects map
        selectWindow("objects_map");
        saveAs("Tiff", resultsDir + seriesName + "_3D_objects.tif");
        close();
        
        // Save the 3D Objects Counter results table
        if (nResults > 0) {
            saveAs("Results", resultsDir + seriesName + "_3D_objects_data.csv");
        }
        
        // Create and save overlay composite
        selectWindow("working");
        saveAs("Tiff", resultsDir + seriesName + "_channels.tif");
        close();
        
    } else {
        print("  No cells detected");
        colocCount = 0;
    }
    
    // Store results in summary with conditions
    summaryResults[row * 10] = day;
    summaryResults[row * 10 + 1] = genotype;
    summaryResults[row * 10 + 2] = sex;
    summaryResults[row * 10 + 3] = replicate;
    summaryResults[row * 10 + 4] = marker;
    summaryResults[row * 10 + 5] = seriesName;
    summaryResults[row * 10 + 6] = cellCount;
    summaryResults[row * 10 + 7] = colocCount;
    if (cellCount > 0) {
        summaryResults[row * 10 + 8] = d2s((colocCount/cellCount)*100, 2);
    } else {
        summaryResults[row * 10 + 8] = "0";
    }
    summaryResults[row * 10 + 9] = "";
    row++;
    
    // Clean up
    close("*");
    run("Clear Results");
    roiManager("reset");
    
    // Force garbage collection to free memory
    run("Collect Garbage");
    
    print("---");
}

// Create summary table with experimental conditions
title = "Summary Results";
Table.create(title);
for (i = 0; i < row; i++) {
    Table.set("Day", i, summaryResults[i * 10]);
    Table.set("Genotype", i, summaryResults[i * 10 + 1]);
    Table.set("Sex", i, summaryResults[i * 10 + 2]);
    Table.set("Replicate", i, summaryResults[i * 10 + 3]);
    Table.set("Marker", i, summaryResults[i * 10 + 4]);
    Table.set("Image", i, summaryResults[i * 10 + 5]);
    Table.set("Total Cells (Ch2)", i, summaryResults[i * 10 + 6]);
    Table.set("Cells with Ch3", i, summaryResults[i * 10 + 7]);
    Table.set("Percentage", i, summaryResults[i * 10 + 8]);
}
Table.update;

// Save summary results
Table.save(resultsDir + "Summary_Results.csv");

print("===");
print("Processing complete!");
print("Results saved to: " + resultsDir);

// Show results folder
if (File.exists(resultsDir)) {
    exec("open", resultsDir);
}
