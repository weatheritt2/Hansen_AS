// Script to load the maxZ .tiffs of the nucleus and ester channel and  and perform illumination correction on them using the BaSIC plugin, using the pre-calculated flat-field and dark-field profiles that were obtained from the focus images from the respective channels. The corrected maxZ images are then outputted as .tifs.

run("Fresh Start"); //Fresh start closes all open images, clears ROI manager and results, it also sets backgrounds to black
output = getDirectory("Select the destination directory")

run("Bio-Formats Macro Extensions");
userChosenDirectory = getDirectory("Select the intial folder"); //SELECT A FOLDER ABOVE THE CHANNEL IMAGES SO (Select PROCESSED_IMAGES)
processBioFormatsFiles (userChosenDirectory);

function processBioFormatsFiles(currentDirectory) {
	function flat_field_max_z(image_folder){
		channel = substring(image_folder, 6, 14);
		positions = newArray("BOTTOM"); //"LOWWER", "MIDDLE"
		correction_profile_folder = "CORRECTION_PROFILES/";
		iteration_count = 0;
		for(i = 0; i <= positions.length - 1; i++){
			files_path = currentDirectory + positions[i] + "/" + image_folder + "/";
			files = getFileList(files_path);
			for(j = 0; j <= files.length - 1; j++){
				run("Bio-Formats Importer", "open=[" + files_path + files[j] + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"); //IMPORTAN, WILL CRASH IF INDEX FILE IS IN DIRECTORY
				iteration_count = iteration_count + 1;
			}	
			stack_name = positions[i]+"_"+channel+"_MAX";
			run("Images to Stack", "name="+stack_name+" title=max use");
			correction_folder_path = currentDirectory + positions[i] + "/" + correction_profile_folder;
			run("Bio-Formats Importer", "open=[" + correction_folder_path + "DARK_FIELD_PROFILE_"+channel+".tiff] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_");
			run("Bio-Formats Importer", "open=[" + correction_folder_path + "FLAT_FIELD_PROFILE_"+channel+".tiff] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_");
			run("BaSiC ", "processing_stack="+stack_name+" flat-field=FLAT_FIELD_PROFILE_"+channel+".tiff dark-field=DARK_FIELD_PROFILE_"+channel+".tiff shading_estimation=[Skip estimation and use predefined shading profiles] shading_model=[Estimate both flat-field and dark-field] setting_regularisationparametes=Automatic temporal_drift=Ignore correction_options=[Compute shading and correct images] lambda_flat=0.50 lambda_dark=0.50");
			close(stack_name);
			close("DARK_FIELD_PROFILE_"+channel+".tiff");
			close("FLAT_FIELD_PROFILE_"+channel+".tiff");
			close("Flat-field:"+stack_name);
			close("Dark-field:"+stack_name);
			selectWindow("Corrected:"+stack_name);
			Stack.setXUnit("um"); //Dimensions are then set for all images per stack (pixel size is 0.148148148um x 0.148148148um)
			run("Properties...", "channels=1 slices="+iteration_count+" frames=1 pixel_width=0.148148148 pixel_height=0.148148148 voxel_depth=10000.0000"); //Of note, this function also requires no. of images in stack
			run("Stack to Images");
			// Replace the corrected images saving section with this code (CODE BELOW BY CLAUDE)
			corrected_images = getList("image.titles");
			Array.sort(corrected_images);
			print("Number of images to process: " + corrected_images.length); // Debug info
			
			for(k = 0; k <= corrected_images.length - 1; k++){
			    image_name = corrected_images[k];
			    
			    // Check if window exists before selecting
			    if (isOpen(image_name)) {
			        selectWindow(image_name);
			        wait(200); // Increased wait time
			        
			        // Store the intended filename before saving
			        save_name = currentDirectory + positions[i] + "/CORRECTED_MAX_Z_" + channel +"/" + image_name + "_corr";
			        print("Saving image " + (k+1) + " of " + corrected_images.length + ": " + save_name); // Debug info
			        
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
flat_field_max_z("MAX_Z_CHANNEL1");
run("Fresh Start"); //Fresh start closes all open images, clears ROI manager and results, it also sets backgrounds to black
flat_field_max_z("MAX_Z_CHANNEL2");
}	

