// Macro to export channels from .lif images
// - Individual channels as separate TIF files (user-selectable)
// - Composite TIF with selected channels only (user-selectable)

// Clean up any open images
run("Close All");
run("Clear Results");

// Prompt user to select .lif file
lifPath = File.openDialog("Select .lif file");
lifDir = File.getParent(lifPath);
lifName = File.getName(lifPath);
lifBaseName = replace(lifName, ".lif", "");

// Define log file path
logFilePath = lifDir + File.separator + lifBaseName + "_export_log.txt";

print("=== Analyzing " + lifName + " ===");

// Get number of series and channel information
run("Bio-Formats Macro Extensions");
Ext.setId(lifPath);
Ext.getSeriesCount(seriesCount);

// Open first series to determine number of channels and channel names
run("Bio-Formats Importer", "open=[" + lifPath + "] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT series_1");
getDimensions(width, height, channels, slices, frames);

// Get channel names
channelNames = newArray(channels);
for (ch = 1; ch <= channels; ch++) {
    Stack.setChannel(ch);
    getLut(reds, greens, blues);
    // Try to get channel name from metadata
    Property.set("Info", getMetadata("Info"));
    channelNames[ch-1] = "Channel " + ch;
}

close("*");

print("Found " + seriesCount + " image(s) with " + channels + " channel(s)");
print("");

// === USER DIALOG: Select export options ===
Dialog.create("Channel Export Options");
Dialog.addMessage("Select channels to export:");
Dialog.addMessage("──────────────────────────────────");

// Individual channel exports
Dialog.addMessage("Export Individual Channels:");
exportIndividual = newArray(channels);
for (ch = 0; ch < channels; ch++) {
    Dialog.addCheckbox(channelNames[ch], true);
}

Dialog.addMessage("──────────────────────────────────");

// Composite channel selection
Dialog.addMessage("Include in Composite TIF:");
includeComposite = newArray(channels);
for (ch = 0; ch < channels; ch++) {
    Dialog.addCheckbox(channelNames[ch], true);
}

Dialog.addMessage("──────────────────────────────────");
Dialog.addCheckbox("Create Composite TIF", true);
Dialog.addMessage("");
Dialog.addCheckbox("Export as Max Projection (instead of Z-stack)", false);

Dialog.show();

// Get user selections for individual exports
for (ch = 0; ch < channels; ch++) {
    exportIndividual[ch] = Dialog.getCheckbox();
}

// Get user selections for composite
for (ch = 0; ch < channels; ch++) {
    includeComposite[ch] = Dialog.getCheckbox();
}

createComposite = Dialog.getCheckbox();
exportMaxProjection = Dialog.getCheckbox();

// Validate selections
anyIndividual = false;
for (ch = 0; ch < channels; ch++) {
    if (exportIndividual[ch]) {
        anyIndividual = true;
    }
}

anyComposite = false;
for (ch = 0; ch < channels; ch++) {
    if (includeComposite[ch]) {
        anyComposite = true;
    }
}

if (!anyIndividual && !createComposite) {
    print("ERROR: No export options selected. Exiting.");
    exit("No export options selected. Please run the macro again.");
}

if (createComposite && !anyComposite) {
    print("ERROR: Composite selected but no channels chosen. Exiting.");
    exit("Composite TIF selected but no channels chosen for composite.");
}

// Create output directories
print("Creating output directories...");

// Individual channel folders
channelDirs = newArray(channels);
for (ch = 0; ch < channels; ch++) {
    if (exportIndividual[ch]) {
        channelDirs[ch] = lifDir + File.separator + lifBaseName + "_Ch" + (ch+1) + File.separator;
        File.makeDirectory(channelDirs[ch]);
        print("  Channel " + (ch+1) + ": " + channelDirs[ch]);
    }
}

// Composite folder
if (createComposite && anyComposite) {
    compositeDir = lifDir + File.separator + lifBaseName + "_Composite" + File.separator;
    File.makeDirectory(compositeDir);
    print("  Composite: " + compositeDir);
}

print("");
print("=== Processing Images ===");
if (exportMaxProjection) {
    print("Export mode: Max Projection");
} else {
    print("Export mode: Full Z-stack");
}
print("");

// Track processing stats
processedCount = 0;
skippedCount = 0;

// Process each series
for (i = 0; i < seriesCount; i++) {
    
    // Open the current series
    run("Bio-Formats Importer", "open=[" + lifPath + "] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT series_" + (i+1));
    
    // Get the current image title
    fullTitle = getTitle();
    
    // Extract image name (after " - ")
    if (indexOf(fullTitle, " - ") > 0) {
        imageName = substring(fullTitle, indexOf(fullTitle, " - ") + 3);
    } else {
        imageName = fullTitle;
        imageName = replace(imageName, ".lif", "");
    }
    
    // Handle nested folders in .lif structure
    lastSlashIndex = lastIndexOf(imageName, "/");
    if (lastSlashIndex > 0) {
        imageName = substring(imageName, lastSlashIndex + 1);
    }
    
    // Replace problematic characters in filename
    imageName = replace(imageName, "/", "_");
    imageName = replace(imageName, "\\", "_");
    imageName = replace(imageName, " ", "_");
    
    print("Processing [" + (i+1) + "/" + seriesCount + "]: " + imageName);
    
    // Get image properties
    getDimensions(width, height, currentChannels, slices, frames);
    print("  Dimensions: " + width + "x" + height + ", " + currentChannels + " channel(s), " + slices + " slice(s)");
    
    // Validate image has enough channels for selected operations
    if (currentChannels < channels) {
        print("  SKIPPED: Image has fewer channels (" + currentChannels + ") than expected (" + channels + ")");
        print("");
        close("*");
        skippedCount++;
        continue;
    }
    
    originalImage = getTitle();
    
    // === EXPORT INDIVIDUAL CHANNELS ===
    for (ch = 0; ch < channels; ch++) {
        if (exportIndividual[ch]) {
            selectWindow(originalImage);
            
            // Duplicate the specific channel
            run("Duplicate...", "duplicate channels=" + (ch+1));
            channelImage = getTitle();
            
            // Apply max projection if selected
            if (exportMaxProjection && slices > 1) {
                run("Z Project...", "projection=[Max Intensity] all");
                close(channelImage);
                channelImage = getTitle();
            }
            
            // Save as TIF
            saveAs("Tiff", channelDirs[ch] + imageName + "_Ch" + (ch+1) + ".tif");
            print("  Saved Channel " + (ch+1) + ": " + imageName + "_Ch" + (ch+1) + ".tif");
            
            close();
        }
    }
    
    // === EXPORT COMPOSITE WITH SELECTED CHANNELS ===
    if (createComposite && anyComposite) {
        selectWindow(originalImage);
        
        // Identify which channels to include
        channelsToKeep = newArray();
        for (ch = 0; ch < channels; ch++) {
            if (includeComposite[ch]) {
                channelsToKeep = Array.concat(channelsToKeep, ch+1);
            }
        }
        
        print("  Creating composite with channels: " + String.join(channelsToKeep, ","));
        
        // Duplicate the full image first
        selectWindow(originalImage);
        run("Duplicate...", "duplicate");
        tempImage = getTitle();
        
        // Apply max projection if needed
        if (exportMaxProjection && slices > 1) {
            run("Z Project...", "projection=[Max Intensity] all");
            close(tempImage);
            tempImage = getTitle();
        }
        
        // Split channels
        selectWindow(tempImage);
        run("Split Channels");
        
        // Close unwanted channels and keep track of channels to merge
        channelImagesToMerge = newArray();
        for (ch = 1; ch <= channels; ch++) {
            channelImageName = "C" + ch + "-" + tempImage;
            shouldKeep = false;
            for (k = 0; k < channelsToKeep.length; k++) {
                if (channelsToKeep[k] == ch) {
                    shouldKeep = true;
                }
            }
            if (shouldKeep) {
                channelImagesToMerge = Array.concat(channelImagesToMerge, channelImageName);
            } else {
                if (isOpen(channelImageName)) {
                    selectWindow(channelImageName);
                    close();
                }
            }
        }
        
        // Merge the selected channels back together if more than 1
        if (channelImagesToMerge.length > 1) {
            // Build merge command
            mergeCmd = "c1=[" + channelImagesToMerge[0] + "]";
            for (m = 1; m < channelImagesToMerge.length; m++) {
                mergeCmd = mergeCmd + " c" + (m+1) + "=[" + channelImagesToMerge[m] + "]";
            }
            mergeCmd = mergeCmd + " create";
            run("Merge Channels...", mergeCmd);
            compositeImage = getTitle();
            
            // Set composite mode and activate all channels
            Stack.setDisplayMode("composite");
            allChannelsString = "";
            for (chn = 1; chn <= channelImagesToMerge.length; chn++) {
                allChannelsString = allChannelsString + "1";
            }
            Stack.setActiveChannels(allChannelsString);
        } else {
            // Only one channel, just use it directly
            compositeImage = channelImagesToMerge[0];
        }
        
        // Save as TIF
        selectWindow(compositeImage);
        saveAs("Tiff", compositeDir + imageName + "_Composite.tif");
        print("  Saved Composite: " + imageName + "_Composite.tif (" + channelImagesToMerge.length + " channel(s))");
        
        close();
    }
    
    // Close original image
    selectWindow(originalImage);
    close();
    
    processedCount++;
    print("");
}

// Close Bio-Formats extensions
Ext.close();

print("=== Processing Complete ===");
print("Processed: " + processedCount + " image(s)");
if (skippedCount > 0) {
    print("Skipped: " + skippedCount + " image(s) (incompatible channel count)");
}
print("");

// Print summary
if (anyIndividual) {
    print("Individual channels saved to:");
    for (ch = 0; ch < channels; ch++) {
        if (exportIndividual[ch]) {
            print("  Channel " + (ch+1) + ": " + channelDirs[ch]);
        }
    }
}

if (createComposite && anyComposite) {
    print("Composite images saved to: " + compositeDir);
}

print("");
print("Log saved to: " + logFilePath);
print("Export complete!");

// Save log to text file
logContent = getInfo("log");
File.saveString(logContent, logFilePath);
