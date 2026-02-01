// General-purpose macro to export .lif images as TIF and PNG
// - TIF: original images
// - PNG: max-projection composite with auto-adjusted contrast and scale bar

// Clean up any open images
run("Close All");
run("Clear Results");

// Prompt user to select .lif file
lifPath = File.openDialog("Select .lif file");
lifDir = File.getParent(lifPath);
lifName = File.getName(lifPath);
lifBaseName = replace(lifName, ".lif", "");

// Create output directories
tifDir = lifDir + File.separator + lifBaseName + "_TIF" + File.separator;
pngDir = lifDir + File.separator + lifBaseName + "_PNG" + File.separator;
File.makeDirectory(tifDir);
File.makeDirectory(pngDir);

print("=== Processing " + lifName + " ===");
print("Output directories:");
print("  TIF: " + tifDir);
print("  PNG: " + pngDir);
print("");

// Get number of series (images) in the .lif file
run("Bio-Formats Macro Extensions");
Ext.setId(lifPath);
Ext.getSeriesCount(seriesCount);

print("Found " + seriesCount + " image(s) in .lif file");
print("");

// Process each series
for (i = 0; i < seriesCount; i++) {
    
    // Open the current series
    run("Bio-Formats Importer", "open=[" + lifPath + "] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT series_" + (i+1));
    
    // Get the current image title
    fullTitle = getTitle();
    
    // Extract image name (after " - ")
    // Format is typically "project-name.lif - ImageName"
    if (indexOf(fullTitle, " - ") > 0) {
        imageName = substring(fullTitle, indexOf(fullTitle, " - ") + 3);
    } else {
        imageName = fullTitle;
        // Remove .lif extension if present
        imageName = replace(imageName, ".lif", "");
    }
    
    // Handle nested folders in .lif structure (e.g., "folder/subfolder/imagename")
    // Extract only the final image name after the last "/"
    lastSlashIndex = lastIndexOf(imageName, "/");
    if (lastSlashIndex > 0) {
        imageName = substring(imageName, lastSlashIndex + 1);
    }
    
    // Replace any remaining problematic characters in filename
    imageName = replace(imageName, "/", "_");
    imageName = replace(imageName, "\\", "_");
    
    print("Processing [" + (i+1) + "/" + seriesCount + "]: " + imageName);
    
    // Get image properties
    getDimensions(width, height, channels, slices, frames);
    print("  Dimensions: " + width + "x" + height + ", " + channels + " channel(s), " + slices + " slice(s)");
    
    // === EXPORT 1: Save as TIF (original) ===
    selectWindow(fullTitle);
    saveAs("Tiff", tifDir + imageName + ".tif");
    tifImage = getTitle();
    
    // === EXPORT 2: Create max projection and save as PNG ===
    
    // Create max projection if Z-stack
    selectWindow(tifImage);
    if (slices > 1) {
        run("Z Project...", "projection=[Max Intensity]");
        maxTitle = getTitle();
        selectWindow(tifImage);
        close();
        selectWindow(maxTitle);
    } else {
        // Already a single slice, just duplicate for processing
        run("Duplicate...", "duplicate");
    }
    
    pngImage = getTitle();
    
    // Convert to composite if not already
    Stack.setDisplayMode("composite");
    
    // Auto-adjust brightness and contrast for each channel
    for (ch = 1; ch <= channels; ch++) {
        Stack.setChannel(ch);
        resetMinAndMax();
        run("Enhance Contrast", "saturated=0.35");
    }
    
    // Make all channels visible - build channel string
    channelString = "";
    for (ch = 1; ch <= channels; ch++) {
        channelString = channelString + "1";
    }
    Stack.setDisplayMode("composite");
    Stack.setActiveChannels(channelString);
    
    // Add scale bar
    // Scale bar parameters: width=50 micrometers, height=4, color=white, background=none, location=lower right
    run("Scale Bar...", "width=50 height=4 font=14 color=White background=None location=[Lower Right] bold");
    
    // Flatten the image to burn in the composite and scale bar
    run("Flatten");
    flatImage = getTitle();
    
    // Save as PNG
    saveAs("PNG", pngDir + imageName + ".png");
    
    print("  Saved TIF: " + imageName + ".tif");
    print("  Saved PNG: " + imageName + ".png");
    print("");
    
    // Close all images for this series
    close("*");
}

// Close Bio-Formats extensions
Ext.close();

print("=== Processing Complete ===");
print("Processed " + seriesCount + " image(s)");
print("TIF files saved to: " + tifDir);
print("PNG files saved to: " + pngDir);
