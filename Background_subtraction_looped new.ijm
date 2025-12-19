// Loop for subracting background from images for BrainJ
// Processes all images in a user-selected folder,
// pausing for manual ROI drawing and rotation.

// 1. Ask user to choose a directory
dir = getDirectory("Choose a Directory");
list = getFileList(dir);

// Set the processing parameters
setBackgroundColor(0, 0, 0);

// 2. Loop through all files in the directory
for (i=0; i<list.length; i++) {
    // 3. Open the image
    open(dir + list[i]);
    
    // 4. Pause for manual user input
    // A dialog will pop up and wait for you to click "OK"
    waitForUser("Action Required", "Draw your ROI and rotate the image as needed, then click OK.");

    // 5. Check if you drew an ROI
    getSelectionBounds(x, y, w, h);
    if (w==0 && h==0) {
        // No ROI was found, so it skips this image
        print("Skipping " + getTitle() + ": No ROI found.");
    } else {
        // An ROI was found, so it processes the image
        print("Processing " + getTitle());
        run("Clear Outside", "stack");
        run("Save");
    }
    
    // 6. Close the image to move to the next one
    // (This runs whether you processed it or skipped it)
    close();
}

showMessage("Processing Complete!");
