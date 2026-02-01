// Export training images from .lif file for Cellpose
// Maintains blinding by using randomized names
// Supports 2D slices or 3D z-stacks with multi-channel export

// Close any open images
close("*");

// User input dialog
Dialog.create("Cellpose Training Image Exporter");
Dialog.addMessage("=== Export Mode ===");
Dialog.addChoice("Export mode:", newArray("Individual 2D slices", "Full 3D z-stacks"));
Dialog.addMessage("\n=== Channel Selection ===");
Dialog.addCheckbox("Channel 1", false);
Dialog.addCheckbox("Channel 2", true);
Dialog.addCheckbox("Channel 3", true);
Dialog.addCheckbox("Channel 4", false);
Dialog.addMessage("\n=== Cropping Options ===");
Dialog.addCheckbox("Apply random cropping", true);
Dialog.addNumber("Crop width (pixels):", 512);
Dialog.addNumber("Crop height (pixels):", 512);
Dialog.addNumber("Number of crops per image:", 1);
Dialog.addMessage("\n=== 2D Slice Options (ignored for 3D mode) ===");
Dialog.addNumber("Slice interval (export every Nth slice):", 5);
Dialog.addMessage("Or specify custom slices (comma-separated, e.g., 5,10,15,20,25):");
Dialog.addString("Custom slices (leave empty for interval):", "", 30);
Dialog.show();

exportMode = Dialog.getChoice();
ch1Export = Dialog.getCheckbox();
ch2Export = Dialog.getCheckbox();
ch3Export = Dialog.getCheckbox();
ch4Export = Dialog.getCheckbox();
applyCropping = Dialog.getCheckbox();
cropWidth = Dialog.getNumber();
cropHeight = Dialog.getNumber();
numCropsPerImage = Dialog.getNumber();
sliceInterval = Dialog.getNumber();
customSlices = Dialog.getString();

is3DMode = (exportMode == "Full 3D z-stacks");

// Count selected channels
selectedChannels = newArray(0);
if (ch1Export) selectedChannels = Array.concat(selectedChannels, 1);
if (ch2Export) selectedChannels = Array.concat(selectedChannels, 2);
if (ch3Export) selectedChannels = Array.concat(selectedChannels, 3);
if (ch4Export) selectedChannels = Array.concat(selectedChannels, 4);

if (selectedChannels.length == 0) {
    exit("No channels selected!");
}

// Prompt user to select the .lif file
lifFile = File.openDialog("Select .lif file");
if (lifFile == "") {
    exit("No file selected");
}

// Prompt user to select output directory
exportDir = getDirectory("Choose directory to save training images");
if (exportDir == "") {
    exit("No output directory selected");
}

if (is3DMode) {
    // For 3D mode, create a single output folder
    if (selectedChannels.length == 1) {
        outputDir = exportDir + "3D_Ch" + selectedChannels[0] + File.separator;
    } else {
        channelString = "";
        for (c = 0; c < selectedChannels.length; c++) {
            channelString += selectedChannels[c];
            if (c < selectedChannels.length - 1) channelString += "-";
        }
        outputDir = exportDir + "3D_Ch" + channelString + File.separator;
    }
    File.makeDirectory(outputDir);
} else {
    // For 2D mode, create channel-specific folders
    if (ch1Export) {
        ch1Dir = exportDir + "Channel_1" + File.separator;
        File.makeDirectory(ch1Dir);
    }
    if (ch2Export) {
        ch2Dir = exportDir + "Channel_2" + File.separator;
        File.makeDirectory(ch2Dir);
    }
    if (ch3Export) {
        ch3Dir = exportDir + "Channel_3" + File.separator;
        File.makeDirectory(ch3Dir);
    }
    if (ch4Export) {
        ch4Dir = exportDir + "Channel_4" + File.separator;
        File.makeDirectory(ch4Dir);
    }
}

// Get the number of series (images) in the .lif file
run("Bio-Formats Macro Extensions");
Ext.setId(lifFile);
Ext.getSeriesCount(seriesCount);

print("\\Clear");
print("Exporting training images from: " + lifFile);
print("Export mode: " + exportMode);
print("Total series: " + seriesCount);
print("Selected channels: " + String.join(selectedChannels, ", "));
if (applyCropping) {
    print("Random cropping: " + cropWidth + "x" + cropHeight + " pixels");
    print("Crops per image: " + numCropsPerImage);
}
print("---");

// Parse custom slices if provided (for 2D mode)
if (!is3DMode) {
    useCustomSlices = (lengthOf(customSlices) > 0);
    if (useCustomSlices) {
        customSliceArray = split(customSlices, ",");
        for (i = 0; i < customSliceArray.length; i++) {
            customSliceArray[i] = parseInt(trim(customSliceArray[i]));
        }
        print("Using custom slices: " + customSlices);
    } else {
        print("Using slice interval: " + sliceInterval);
    }
    print("---");
}

// Generate random image IDs for blinding
random("seed", getTime());
imageIDs = newArray(seriesCount);
for (i = 0; i < seriesCount; i++) {
    imageIDs[i] = round(random * 100000);
}

// Global slice counter for unique IDs across all series
globalSliceCounter = 0;

// Loop through each series
for (s = 0; s < seriesCount; s++) {
    
    print("Processing series " + (s+1) + " of " + seriesCount);
    
    // Open the image with Bio-Formats
    run("Bio-Formats Importer", 
        "open=[" + lifFile + "] " +
        "autoscale " +
        "color_mode=Grayscale " +
        "view=Hyperstack " +
        "stack_order=XYCZT " +
        "series_" + (s+1));
    
    originalTitle = getTitle();
    getDimensions(width, height, channels, slices, frames);
    
    print("  Dimensions: " + width + "x" + height + ", " + channels + " channels, " + slices + " slices");
    
    // Determine if cropping is possible for this series
    canCropThisSeries = true;
    if (applyCropping && is3DMode) {
        maxX = width - cropWidth;
        maxY = height - cropHeight;
        if (maxX < 0 || maxY < 0) {
            print("  WARNING: Image too small for cropping (" + width + "x" + height + "), skipping crop");
            canCropThisSeries = false;
        }
    }
    
    if (is3DMode) {
        // Export entire z-stack with selected channels (with multiple crops if requested)
        
        // Loop for multiple crops
        for (cropNum = 0; cropNum < numCropsPerImage; cropNum++) {
            
            // Generate random crop coordinates for this crop iteration
            if (applyCropping && canCropThisSeries) {
                maxX = width - cropWidth;
                maxY = height - cropHeight;
                cropX = round(random * maxX);
                cropY = round(random * maxY);
                applyCroppingThisCrop = true;
                if (numCropsPerImage > 1) {
                    print("  3D crop " + (cropNum+1) + "/" + numCropsPerImage + " region: " + cropX + "," + cropY + " to " + (cropX+cropWidth) + "," + (cropY+cropHeight));
                } else {
                    print("  3D crop region: " + cropX + "," + cropY + " to " + (cropX+cropWidth) + "," + (cropY+cropHeight));
                }
            } else {
                applyCroppingThisCrop = false;
            }
            
            if (selectedChannels.length == 1) {
                // Single channel z-stack
                ch = selectedChannels[0];
                selectWindow(originalTitle);
                run("Duplicate...", "title=export_stack duplicate channels=" + ch);
                
                // Apply crop if enabled
                if (applyCroppingThisCrop) {
                    makeRectangle(cropX, cropY, cropWidth, cropHeight);
                    run("Crop");
                }
                
                // Save with randomized name (include crop number if multiple crops)
                if (numCropsPerImage > 1) {
                    filename = "train_" + IJ.pad(imageIDs[s], 6) + "_crop" + (cropNum+1) + ".tif";
                } else {
                    filename = "train_" + IJ.pad(imageIDs[s], 6) + ".tif";
                }
                saveAs("Tiff", outputDir + filename);
                close();
                
            } else {
                // Multi-channel z-stack
                // Need to merge channels in correct order
                
                // Extract each selected channel
                channelNames = newArray(selectedChannels.length);
                for (c = 0; c < selectedChannels.length; c++) {
                    selectWindow(originalTitle);
                    run("Duplicate...", "title=ch" + selectedChannels[c] + "_temp duplicate channels=" + selectedChannels[c]);
                    channelNames[c] = "ch" + selectedChannels[c] + "_temp";
                }
                
                // Merge channels - build merge command
                mergeCmd = "c1=[" + channelNames[0] + "]";
                for (c = 1; c < selectedChannels.length; c++) {
                    mergeCmd += " c" + (c+1) + "=[" + channelNames[c] + "]";
                }
                mergeCmd += " create";
                
                run("Merge Channels...", mergeCmd);
                
                // Rename and ensure correct dimension order (Z x C x Y x X)
                rename("export_stack");
                
                // Apply crop if enabled (to entire merged stack)
                if (applyCroppingThisCrop) {
                    makeRectangle(cropX, cropY, cropWidth, cropHeight);
                    run("Crop");
                }
                
                // Save with randomized name (include crop number if multiple crops)
                if (numCropsPerImage > 1) {
                    filename = "train_" + IJ.pad(imageIDs[s], 6) + "_crop" + (cropNum+1) + ".tif";
                } else {
                    filename = "train_" + IJ.pad(imageIDs[s], 6) + ".tif";
                }
                saveAs("Tiff", outputDir + filename);
                close();
            }
            
            if (numCropsPerImage > 1) {
                print("  Exported 3D stack crop " + (cropNum+1) + ": " + filename);
            } else {
                print("  Exported 3D stack: " + filename);
            }
        }
        
    } else {
        // 2D mode - export individual slices
        
        // Determine which slices to export
        if (useCustomSlices) {
            slicesToExport = newArray(0);
            for (i = 0; i < customSliceArray.length; i++) {
                if (customSliceArray[i] <= slices) {
                    slicesToExport = Array.concat(slicesToExport, customSliceArray[i]);
                }
            }
        } else {
            // Use interval
            numExportSlices = floor(slices / sliceInterval);
            slicesToExport = newArray(numExportSlices);
            for (i = 0; i < numExportSlices; i++) {
                slicesToExport[i] = (i + 1) * sliceInterval;
                if (slicesToExport[i] > slices) {
                    slicesToExport[i] = slices;
                }
            }
        }
        
        print("  Exporting " + slicesToExport.length + " slices");
        
        // Export each selected slice for each selected channel (with multiple crops if requested)
        for (i = 0; i < slicesToExport.length; i++) {
            z = slicesToExport[i];
            
            // Loop for multiple crops per slice
            for (cropNum = 0; cropNum < numCropsPerImage; cropNum++) {
                sliceID = globalSliceCounter;
                globalSliceCounter++;
                
                // Generate random crop coordinates for this slice and crop iteration (2D mode)
                if (applyCropping) {
                    maxX = width - cropWidth;
                    maxY = height - cropHeight;
                    if (maxX < 0 || maxY < 0) {
                        sliceCropX = 0;
                        sliceCropY = 0;
                        applyCroppingThisSlice = false;
                    } else {
                        sliceCropX = round(random * maxX);
                        sliceCropY = round(random * maxY);
                        applyCroppingThisSlice = true;
                    }
                } else {
                    applyCroppingThisSlice = false;
                }
                
                // Channel 1
                if (ch1Export && channels >= 1) {
                    selectWindow(originalTitle);
                    Stack.setPosition(1, z, 1);
                    run("Duplicate...", "title=temp_export");
                    if (applyCroppingThisSlice) {
                        makeRectangle(sliceCropX, sliceCropY, cropWidth, cropHeight);
                        run("Crop");
                    }
                    filename = "train_" + IJ.pad(sliceID, 6) + ".tif";
                    saveAs("Tiff", ch1Dir + filename);
                    close();
                }
                
                // Channel 2
                if (ch2Export && channels >= 2) {
                    selectWindow(originalTitle);
                    Stack.setPosition(2, z, 1);
                    run("Duplicate...", "title=temp_export");
                    if (applyCroppingThisSlice) {
                        makeRectangle(sliceCropX, sliceCropY, cropWidth, cropHeight);
                        run("Crop");
                    }
                    filename = "train_" + IJ.pad(sliceID, 6) + ".tif";
                    saveAs("Tiff", ch2Dir + filename);
                    close();
                }
                
                // Channel 3
                if (ch3Export && channels >= 3) {
                    selectWindow(originalTitle);
                    Stack.setPosition(3, z, 1);
                    run("Duplicate...", "title=temp_export");
                    if (applyCroppingThisSlice) {
                        makeRectangle(sliceCropX, sliceCropY, cropWidth, cropHeight);
                        run("Crop");
                    }
                    filename = "train_" + IJ.pad(sliceID, 6) + ".tif";
                    saveAs("Tiff", ch3Dir + filename);
                    close();
                }
                
                // Channel 4
                if (ch4Export && channels >= 4) {
                    selectWindow(originalTitle);
                    Stack.setPosition(4, z, 1);
                    run("Duplicate...", "title=temp_export");
                    if (applyCroppingThisSlice) {
                        makeRectangle(sliceCropX, sliceCropY, cropWidth, cropHeight);
                        run("Crop");
                    }
                    filename = "train_" + IJ.pad(sliceID, 6) + ".tif";
                    saveAs("Tiff", ch4Dir + filename);
                    close();
                }
            }
        }
        
        if (numCropsPerImage > 1) {
            print("  Exported " + slicesToExport.length + " slices x " + numCropsPerImage + " crops = " + (slicesToExport.length * numCropsPerImage) + " images per channel");
        } else {
            print("  Exported " + slicesToExport.length + " slices per channel");
        }
    }
    
    close(originalTitle);
}

print("---");
print("Export complete!");
print("Images saved to: " + exportDir);

// Open export directory
if (File.exists(exportDir)) {
    exec("open", exportDir);
}