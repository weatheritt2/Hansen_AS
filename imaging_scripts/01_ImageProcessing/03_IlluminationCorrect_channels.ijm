// Script loads the identified focus slices for all channels and use all focus images across the fields-of-view to calculate illumination correction profiles using the BaSIC plugin. Both flat-field and dark-field profiles are calculated. The calculated profiles are then applied to the focus images to generate illumination corrected images.
// The illuminated corrected images are saved as .tiffs.
// User is prompted to specify input and output directories upon running the script.

run("Fresh Start"); //Fresh start closes all open images, clears ROI manager and results, it also sets backgrounds to black
output = getDirectory("Select the destination directory")

run("Bio-Formats Macro Extensions");
userChosenDirectory = getDirectory("Select the intial folder"); //SELECT A FOLDER ABOVE THE CHANNEL IMAGES SO I
processBioFormatsFiles (userChosenDirectory);

//Not sure what this wrapped function does
function processBioFormatsFiles(currentDirectory) {
	function illumination_correct(channel_name){
		files_BOTTOM = currentDirectory + "BOTTOM/";
		//files_LOWWER = currentDirectory + "LOWWER\\";
		//files_MIDDLE = currentDirectory + "MIDDLE\\";
		iteration_count = 0;
		channels = newArray(channel_name); //newArray("CHANNEL1"); //"CHANNEL2","CHANNEL3","CHANNEL4",
		for(ch = 0; ch <= channels.length - 1; ch++){
		files_path = files_BOTTOM + channels[ch];
		files = getFileList(files_path);
			for(i = 0; i <= files.length - 1; i++){
				run("Bio-Formats Importer", "open=[" + files_path +"/"+ files[i] + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"); //IMPORTAN, WILL CRASH IF INDEX FILE IS IN DIRECTORY
			}
			run("Images to Stack", "name="+channels[ch]+" title=_ch"+substring(channels[ch],7)+" use");
			wait(100); // Add a 100ms pause before saving
			run("Collect Garbage"); //SHOULD RE-CLAIM UNUSED RAM BUT MAY NOT WORK, WE SHALL SEE
			run("BaSiC ", "processing_stack="+channels[ch]+" flat-field=None dark-field=None shading_estimation=[Estimate shading profiles] shading_model=[Estimate both flat-field and dark-field] setting_regularisationparametes=Automatic temporal_drift=Ignore correction_options=[Compute shading and correct images] lambda_flat=0.50 lambda_dark=0.50");
			wait(100); // Add a 100ms pause before saving
			close(channels[ch]);
			wait(100); // Add a 100ms pause before saving
			selectWindow("Flat-field:"+channels[ch]);
			wait(100); // Add a 100ms pause before saving
			saveAs("Tiff", files_BOTTOM + "/CORRECTION_PROFILES/FLAT_FIELD_PROFILE_"+channels[ch]+".tiff");
			close();
			wait(100); // Add a 100ms pause before saving
			selectWindow("Dark-field:"+channels[ch]);
			wait(100); // Add a 100ms pause before saving
			saveAs("Tiff", files_BOTTOM + "/CORRECTION_PROFILES/DARK_FIELD_PROFILE_"+channels[ch]+".tiff");
			wait(100); // Add a 100ms pause before saving
			close();
			selectWindow("Corrected:"+channels[ch]);
			wait(100); // Add a 100ms pause before saving
			Stack.setXUnit("um"); //Dimensions are then set for all images per stack (pixel size is 0.148148148um x 0.148148148um)
			run("Properties...", "channels=1 slices="+(files.length)+" frames=1 pixel_width=0.148148148 pixel_height=0.148148148 voxel_depth=10000.0000"); //Of note, this function also requires no. of images in stack
			run("Stack to Images");
			corrected_images = getList("image.titles");
			Array.sort(corrected_images);
			// Replace the corrected images saving section with this code (CODE BELOW BY CLAUDE)
			corrected_images = getList("image.titles");
			Array.sort(corrected_images);
			print("Number of images to process: " + corrected_images.length); // Debug info
			
			for(i = 0; i <= corrected_images.length - 1; i++){
			    image_name = corrected_images[i];
			    
			    // Check if window exists before selecting
			    if (isOpen(image_name)) {
			        selectWindow(image_name);
			        wait(200); // Increased wait time
			        
			        // Store the intended filename before saving
			        save_name = files_BOTTOM + "/CORRECTED_"+channels[ch]+"/" + substring(image_name,0,16)+"_corr";
			        print("Saving image " + (i+1) + " of " + corrected_images.length + ": " + save_name); // Debug info
			        
			        // Save with error checking
			        if (File.exists(save_name)) {
			            File.delete(save_name); // Remove any existing file to prevent issues
			        }
			        saveAs("Tiff", save_name);
			        wait(300); // Increased wait time after saving
			        
			        // Verify save was successful
			        if (File.exists(save_name)) {
			            close(); // Only close if save was successful
			            wait(100);
			        } else {
			            print("Warning: Failed to save " + save_name);
			        }
			    } else {
			        print("Warning: Window " + image_name + " not found");
			    }
			    run("Collect Garbage");
			}
			
			// Add verification step
			run("Close All");
			wait(200);
			}}
illumination_correct("CHANNEL1")
illumination_correct("CHANNEL2")
illumination_correct("CHANNEL3")
illumination_correct("CHANNEL4")
}

