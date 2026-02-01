/*
 * Spinal Cord Grey Matter ROI Creator
 * 
 * This macro processes TIF images of spinal dorsal horn sections.
 * For each image, the user draws a freehand line separating white matter 
 * (top) from grey matter (bottom). The macro then:
 * 1. Saves the drawn line as an ROI file (for verification/adjustment)
 * 2. Creates and saves a binary mask (1=grey matter to keep, 0=white matter to remove)
 * 
 * The binary masks are 32-bit with 0-1 values for proper multiplication.
 * This ensures that multiplying with original images preserves pixel values.
 * 
 * Expected structure:
 * Input_Folder/
 *   - image1.tif
 *   - image2.tif
 * Output will create:
 * Input_Folder/ROIs/
 *   - image1_line.roi
 *   - image1_mask.tif
 *   - image2_line.roi
 *   - image2_mask.tif
 */

// Clean up before starting
run("Close All");
run("Clear Results");
roiManager("reset");
print("\\Clear");

// Get input directory
inputDir = getDirectory("Choose folder containing TIF images");
outputROIDir = inputDir + "ROIs" + File.separator;
File.makeDirectory(outputROIDir);

// Get list of TIF files
fileList = getFileList(inputDir);
tifFiles = newArray();
for (i = 0; i < fileList.length; i++) {
    if (endsWith(fileList[i], ".tif") || endsWith(fileList[i], ".tiff")) {
        tifFiles = Array.concat(tifFiles, fileList[i]);
    }
}

if (tifFiles.length == 0) {
    exit("No TIF files found in the selected directory!");
}

print("Found " + tifFiles.length + " TIF images to process.");
print("Output will be saved in: " + outputROIDir);
print("---");

// Set tool to freehand line
setTool("freeline");

// Process each image
for (i = 0; i < tifFiles.length; i++) {
    // Open image
    imagePath = inputDir + tifFiles[i];
    open(imagePath);
    originalID = getImageID();
    originalTitle = getTitle();
    baseName = File.nameWithoutExtension;
    
    print("Processing: " + originalTitle + " (" + (i+1) + "/" + tifFiles.length + ")");
    
    // Check if ROI already exists (skip if user wants)
    lineROIPath = outputROIDir + baseName + "_line.roi";
    maskPath = outputROIDir + baseName + "_mask.tif";
    
    if (File.exists(lineROIPath) && File.exists(maskPath)) {
        Dialog.create("ROI Already Exists");
        Dialog.addMessage("ROI and mask already exist for:\n" + originalTitle);
        Dialog.addChoice("Action:", newArray("Skip this image", "Redraw ROI", "Cancel macro"));
        Dialog.show();
        action = Dialog.getChoice();
        
        if (action == "Skip this image") {
            print("  -> Skipped (ROI already exists)");
            close();
            continue;
        } else if (action == "Cancel macro") {
            print("Macro cancelled by user.");
            close();
            exit();
        }
    }
    
    // Instructions for user
    if (i == 0) {
        waitForUser("Draw White-Grey Matter Boundary", 
            "For each image:\n \n" +
            "1. Draw a FREEHAND LINE across the image\n" +
            "   separating white matter (TOP) from grey matter (BOTTOM)\n \n" +
            "2. The line should span the entire width of the image\n \n" +
            "3. Click OK when you're done drawing\n \n" +
            "The macro will:\n" +
            "- Save your line as an ROI file\n" +
            "- Create a binary mask (white=keep, black=remove)\n" +
            "- Save the mask for analysis");
    }
    
    // Wait for user to draw line
    waitForUser("Draw Boundary Line", 
        "Draw the freehand line separating white (TOP) from grey matter (BOTTOM)\n \n" +
        "Image: " + originalTitle + "\n" +
        "Progress: " + (i+1) + "/" + tifFiles.length + "\n \n" +
        "Click OK when ready.");
    
    // Check if a line was drawn
    if (selectionType() != 6) {  // 6 = freehand line
        showMessage("Error", "Please draw a freehand line selection!\nSkipping this image.");
        print("  -> ERROR: No freehand line drawn. Skipped.");
        close();
        continue;
    }
    
    // Save the line ROI
    roiManager("Add");
    roiManager("Select", 0);
    roiManager("Save", lineROIPath);
    print("  -> Saved line ROI: " + baseName + "_line.roi");
    
    // Create binary mask
    // Get image dimensions
    getDimensions(width, height, channels, slices, frames);
    
    // Get line coordinates
    getSelectionCoordinates(xpoints, ypoints);
    
    // Create a binary mask image (32-bit for 0-1 values)
    newImage("Mask", "32-bit black", width, height, 1);
    maskID = getImageID();
    
    // Fill everything below the line with 1 (grey matter = keep)
    // Strategy: for each x position, find the y coordinate on the line
    // and fill from that y coordinate to the bottom
    
    setForegroundColor(1, 1, 1);
    
    // For each column in the image
    for (x = 0; x < width; x++) {
        // Find the y value at this x position on the line
        // Use linear interpolation between line points
        yAtX = height; // default to bottom if x is beyond line range
        
        for (j = 0; j < xpoints.length - 1; j++) {
            x1 = xpoints[j];
            x2 = xpoints[j + 1];
            y1 = ypoints[j];
            y2 = ypoints[j + 1];
            
            // Check if x is between these two points
            if ((x >= x1 && x <= x2) || (x >= x2 && x <= x1)) {
                // Linear interpolation
                if (x2 != x1) {
                    yAtX = y1 + (y2 - y1) * (x - x1) / (x2 - x1);
                } else {
                    yAtX = y1;
                }
                break;
            }
        }
        
        // Draw a line from yAtX to bottom of image for this x position
        if (yAtX < height) {
            makeLine(x, yAtX, x, height-1);
            run("Draw", "slice");
        }
    }
    
    // Clean up the mask - fill any gaps
    run("Select None");
    run("Fill Holes");
    
    // Save the mask
    saveAs("Tiff", maskPath);
    print("  -> Saved mask: " + baseName + "_mask.tif");
    close(); // close mask
    
    // Clear ROI manager
    roiManager("reset");
    
    // Close original image
    selectImage(originalID);
    close();
    
    print("  -> Complete!");
}

print("---");
print("Processing complete!");
print("Total images processed: " + tifFiles.length);
print("Output location: " + outputROIDir);
print("");
print("Next steps for analysis:");
print("1. Use the *_mask.tif files to mask your channels");
print("2. Method: Use 'Image Calculator' to multiply your channel by the mask");
print("3. Then measure IntDen and Area on the masked channel");
