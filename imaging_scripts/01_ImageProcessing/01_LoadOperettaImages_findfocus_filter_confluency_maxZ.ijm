 //This Macro will iteratively load the single slice images of each FOV that is in a given directory. For each FOV, the macro will create a Z-stack of the image slices, and use the intensity variance of the specified channel (the ESTER one) to identify which slice in the stack is in focus. It will then record this slice name and use it as the focus slice for all channels for the given FOV. It will save write this to a .txt file.
 // The macro will also create a subset max-Z projection of specified channels, which is used to create a max-Z of the nucleus and cell body channels for segmentation purposes. These maxZs are saved as single channel .tiffs.
 // The number of slices in a FOV has to be specified. Also, the channel used for variance-based focus selection has to be specified (this operation requires the find_focus macro). The variance-threshold used has to be manually specified before the script is run.
 // The script also removes FOVs without cells and which are confluent (based on how many % of the image has pixel values that correspond to foreground) - these values also need to be manually specified.
// After setting the variables listed above, you can run the script. When running it, you will be prompted to specify input and output directories.

//Define the output folder
run("Fresh Start"); //Fresh start closes all open images, clears ROI manager and results, it also sets backgrounds to black
output = getDirectory("Select the destination directory")

//Ask the user to select directory containing input images
run("Bio-Formats Macro Extensions");
userChosenDirectory = getDirectory("Select the intial folder");
processBioFormatsFiles (userChosenDirectory);

//
function processBioFormatsFiles(currentDirectory) {
	start_time = getTime();
	//Generates a list (array) of the files in the input directory
	fileList = getFileList(currentDirectory);
	//Initial section with IF conditions generates five different sub-lists of the directory file list depending on what channel the file corresponds to (or if the file is not labelled according to the PerkinElmer Operetta standard)
	channel1 = newArray(); //Creates an empty array (what would be a vector in R) that is used for iteration
	channel2 = newArray();
	channel3 = newArray();
	channel4 = newArray();
	other_file = newArray();
	for(file = 0; file <= fileList.length-1; file++){ //THIS FOR LOOP IS TO CREATE DIFFERENT ARRAYS (LISTS) OF FILE NAMES ACCORDING TO CHANNEL NAMING IN ORDER TO BE ABLE ONLY TO IMPORT IMAGES FROM A GIVEN CHANNEL
		x = newArray(fileList[file]);
		if(x[0].length == 30){ //Only does the following if the file name of x is 30 characters long
		ID_SUBSTRING = substring(x[0],13,16); // Selectes the characters of position 13-16 and creates a string. (Note, 0-indexing is used, so it corresponds to character 14 to 16 (the last character in the substring is EXCLUDED)
		} else {
			ID_SUBSTRING = "OTHER_FILE_TYPE";} //If the file name is NOT 30 characters long, the file is assigned to OTHER FILE TYPE LIST
		if(ID_SUBSTRING == "ch1"){ //Creates the list of channel 1
			// NOTE, I have change ch1 images to channel2 array AS THE ORIGINAL SCRIPT had ESTER in channel 2, but in this image series they are ch1
			channel2 = Array.concat(channel2,x[0]); //Array.concat binds together two arrays. This is why we need an empty array for starting the iteration using the for loop.
		} 
		if(ID_SUBSTRING == "ch2"){ //Creates the list of channel 2
			channel3 = Array.concat(channel3, x[0]);
		}
		if(ID_SUBSTRING == "ch3"){ //Creates the list of channel 3
			// NOTE, I have change ch1 images to channel2 array AS THE ORIGINAL SCRIPT had ESTER in channel 2, but in this image series they are ch1
			channel4 = Array.concat(channel4, x[0]);
		}
		if(ID_SUBSTRING == "ch4"){  //Creates the list of channel 4
			// NOTE, I have change ch1 images to channel2 array AS THE ORIGINAL SCRIPT had ESTER in channel 2, but in this image series they are ch1
			channel1 = Array.concat(channel1, x[0]);
		}
		if(ID_SUBSTRING == "OTHER_FILE_TYPE"){
			other_file = Array.concat(other_file,x[0]);
		}
		else{
			continue;
		}}
			

//For loop uses the sub-list generated for channel 2 (which corresponds to my NHS ester images)
slices_per_stack = 8; //Various parts of the pipeline depends on the no. of slices per stack. THIS CAN BE MODIFIED WITH THIS VARIABLE.
channel2_BOT_FOVs = newArray();
channel2_LOW_FOVs = newArray();
channel2_MID_FOVs = newArray();
CONFLUENT_channel2_BOT_FOVs = newArray();
CONFLUENT_channel2_LOW_FOVs = newArray();
CONFLUENT_channel2_MID_FOVs = newArray();
EMPTY_FOVs = newArray();
for (image_ch2 = 0; image_ch2 <= channel2.length-1; image_ch2++) { //Note, the for loop is asked to run for the length of the fileList -1 (as 0-indexing is used when accessing the files).
		iteration_count = image_ch2*1 + 1; //This variable is to count what iteration of the for loop we're on.
		iteration_stack = iteration_count % slices_per_stack; //This variable is used to determine when we have opened all images that make up a stack. 14 is used as the determinant as each of my stacks are 14 images tall.
		position_number = floor(iteration_count/slices_per_stack); //This variable counts the position of the FOV within the and also depends on no. of images in a stack. NOTE UNCOMPLETE, it has to reset when reaching something (corresponding to transition to new well)
		print("Iteration: " + iteration_count);  //Just to see what iteration count we're on
		print("Iteration_stack: " + iteration_stack); //And to see what stack/FOV in a well we're in. Again, I doesn't reset when the next well is reached and 
		//The following function opens the files in imageJ. Note, FIJI will CRASH if file labelled index is in input directory
		run("Bio-Formats Importer", "open=[" + currentDirectory + channel2[image_ch2] + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"); //IMPORTAN, WILL CRASH IF INDEX FILE IS IN DIRECTORY
		if(iteration_stack == 0){  //Iteration stack == 0 every time 14 images have been opened, which corresponds to all images from a stack
			print("Position: ", position_number," END OF STACK"); //Sanity check to see that the condition wroks
			well_pos = substring(channel2[image_ch2],0,9); //Takes the row, col, and FOV data from the file to use as subsequent name
			run("Images to Stack", "name="+well_pos+" title="+well_pos+" use"); //Stacks all open images and gives the stack the ROW/COL/FOV name (this generates the stack from a given FOV).
			//Stack.setXUnit("um"); //Dimensions are then set for all images per stack (pixel size is 0.148148148um x 0.148148148um)
			//run("Properties...", "channels=1 slices="+slices_per_stack+" frames=1 pixel_width=0.148148148 pixel_height=0.148148148 voxel_depth=10000.0000"); //Of note, this function also requires no. of images in stack
			slice_order = newArray(); //empty variable used for iteration below (used to verify that images are stacked in correct order. sanity check).
			//Tje for loop and associated IF/ELSE condition below is to verify that the images making up a stack are put together in correct order, otherwise an error is raised and the script stops
			for(i = 1; i <= slices_per_stack; i++){ //LENGTH IS THE NO. OF SLICES IN STACK
				setSlice(i); //Selects slice i in stack
				slice_name = getInfo("slice.label"); //Gets label from the given slice in the stack and stores as slice_name
				slice_number = newArray(substring(slice_name,10,12)); //the information about what no. in the stack that an image had upon image aqcuistion is held in index 10-12. This is extracted and stored as slice_number
				slice_order = Array.concat(slice_order,slice_number); //The slice numbers are collected in slice_order array through iterations over the labels from the stack. 
				//If stacked correctly, the slice number will be ordered from: 01,02,03 etc. to 14.
			}//If condition below tests if all 14 slices are correctly stacked
			if(slice_order[0] == "01" && slice_order[1] == "02" && slice_order[2] == "03" && slice_order[3] == "04" && slice_order[4] == "05" && slice_order[5] == "06" && slice_order[6] == "07" && slice_order[7] == "08"){ 
				print("IMAGES STACKED CORRECTLY");
			}else{ //If the images are not correctly stacked, the macro is exited and the following message is generated.
				exit("IMAGE STACKING ORDER ERROR OCCURED FOR LAST OPEN STACK: \n LOOK INTO ORDERING OF FILES IN DIRECTORY IMAGE NAMES OR STACKING MACRO");
			}
//If the script doesn't terminate, the stack with the correctly stacked images are then selected
			selectWindow(well_pos);
			focused_slice_variance_THRESHOLD = 40;
			run("Find focused slices", "select="+focused_slice_variance_THRESHOLD+" variance=0.000 select_only verbose"); //The function Find focused slices is then run on the stack. The variance is set at 30%, emperically determined from looking at 6 different wells in SET 1, REP1. Note, this may be a good thing to verify for each well plate (so each replicate, that this no. makes sense).
			selectWindow(well_pos); //The original stacked image is then selected and closed.
			close();
			selectWindow("Focused slices of "+well_pos+"_"+focused_slice_variance_THRESHOLD+".0%"); //The stack of images in focus is then selected. This function doesn't like the padded values in terms of slice no. however and the slices are therefore not properly ordered, which is why the following is done as it is.
			run("Stack to Images"); //The stack of focussed images in converted to individual images
			focus_images = getList("image.titles"); //An array (vector of strings) is then generated that has the title of all open images.
			Array.sort(focus_images); //The array with titles is then sorted according to characters/numerics which gives us the slices in the correct order
			bottom_image = focus_images[0]; //The first slice in focus (index 0) is selected and labelled bottom_image. This corresponds to the first image that is deemed to be in focus and is the plane where the cells adhere to the surface.
			//lower_image = focus_images[1]; //The second slice in focus (index 1) is selected and labelled lower_image. This corresponds to the second image that is deemed to be in focus and is 0.5um above the plane where the cells adhere to the surface.
			//middle_image = focus_images[2]; //An additional image labelled middle image is also selected, which corresponds to two positions above the bottom slice (two index positions higer, so index 2). This corresponds to 1um when step size is 0.5um.
			selectWindow(bottom_image); //Selects bottom image
			run("Duplicate...", " "); //Duplicates selected image
			run("Threshold..."); //Opens the threshold tool
			//NOTE, THIS THRESHOLD IS ABRITRARY AND DETERMINED EMPERICALLY. THE LOWER VALUE WILL DETERMINE WHAT THE REQUIRED THRESHOLD IS. Currently, it is set to 285.
			minimum_pixel_THRESHOLD = 300;
			relative_coverage_lower_THRESHOLD = 1;
			setThreshold(minimum_pixel_THRESHOLD, 65535, "raw"); //Sets a threshhold on the duplicated image where everything with pixel values above 285 are set as foreground and everything below is background
			run("Convert to Mask"); //The thresholded image is converted to mask so all values above 285 are set to 255 and everything below is set to 0 (the mask is 8-bit)
			run("Create Selection"); //Selects foreground (pixels with value of 255) in selected image, which is the mask
			run("Set Measurements...", "area stack limit redirect=None decimal=3"); // TEST
			run("Measure"); //Measures the selection: we need the area measurement as that will tell us how much area is determined to be foreground
			close(); //The mask is then closed
			coverage_area = getResult("Area", 0); //From the measurement, the area measurement is extracted and stored as variable coverage area (note, the 0 in getResult specifies that it is the first row)
			relative_coverage = (coverage_area/194122498.484 * 100); //The relative coverage is then calculated, which is defined as the coverage_area divided by total image area (102400) * 100.
			print("The relative coverage is: "+relative_coverage); //Prints the coverage area (not needed for anything but seeing the process
			close("Results"); //Closes result window so can be run again for next iteration
			close("Threshold"); //Closes threshold window so it can be run again for next iteration
			
			if(relative_coverage >= relative_coverage_lower_THRESHOLD && relative_coverage <= 100){

				//lower_image = focus_images[1]; //The second slice in focus (index 1) is selected and labelled lower_image. This corresponds to the second image that is deemed to be in focus and is 0.5um above the plane where the cells adhere to the surface.
				//middle_image = focus_images[2]; //An additional image labelled middle image is also selected, which corresponds to two positions above the bottom slice (two index positions higer, so index 2). This corresponds to 1um when step size is 0.5um.
				bottom_stack = substring(bottom_image,0,12)+"_ch2_STACK";
				//lower_stack = substring(lower_image,0,12)+"_ch2_STACK";
				//middle_stack = substring(middle_image,0,12)+"_ch2_STACK";
				for(i = 0; i <= focus_images.length-1; i++){
						image = focus_images[i];
						if(i <= 1){
							print("Part of max Z-stack");
						}
						if(i > 1){
							selectWindow(image);
							close();
						}}
			//The following if conditions allows saving images above/below confluence threshold in different directories as images with too high confluency will be hard to segment 
			
			relative_coverage_THRESHOLD = 95; 	//Currently, the threshhold is set to 95 <- AGAIN, ARBITRARY
			if(relative_coverage >= relative_coverage_lower_THRESHOLD && relative_coverage < relative_coverage_THRESHOLD){
				channel2_BOT_FOVs = Array.concat(channel2_BOT_FOVs, substring(bottom_image,0,12)+"_ch2.tiff");
				//channel2_LOW_FOVs = Array.concat(channel2_LOW_FOVs, substring(lower_image,0,12)+"_ch2.tiff");
				//channel2_MID_FOVs = Array.concat(channel2_MID_FOVs, substring(middle_image,0,12)+"_ch2.tiff");
				print("BELOW "+relative_coverage_THRESHOLD+"% COVERAGE");
				//
				run("Images to Stack", "name="+bottom_stack+" use keep");		
				selectWindow(focus_images[0]);
				close();
				//run("Images to Stack", "name="+lower_stack+" use keep");
				//selectWindow(focus_images[1]);
				//close();
				//run("Images to Stack", "name="+middle_stack+" use keep");
				//
				selectWindow(bottom_stack);
				run("Z Project...", "projection=[Max Intensity]");
				selectWindow("MAX_"+bottom_stack);
				saveAs("Tiff", output + "/BOTTOM//MAX_Z_CHANNEL2/"+substring(bottom_stack,0,16)+"_max");
				//
				//selectWindow(lower_stack);
				//run("Z Project...", "projection=[Max Intensity]");
				//selectWindow("MAX_"+lower_stack);
				//saveAs("Tiff", output + "\\LOWWER\\MAX_Z_CHANNEL2\\"+substring(lower_stack,0,16)+"_max");
				//
				//selectWindow(middle_stack);
				//run("Z Project...", "projection=[Max Intensity]");
				//selectWindow("MAX_"+middle_stack);
				//saveAs("Tiff", output + "\\MIDDLE\\MAX_Z_CHANNEL2\\"+substring(middle_stack,0,16)+"_max");
					
			}if(relative_coverage >= relative_coverage_THRESHOLD && relative_coverage <= 100){
				CONFLUENT_channel2_BOT_FOVs = Array.concat(CONFLUENT_channel2_BOT_FOVs, substring(bottom_image,0,12)+"_ch2.tiff");
				//CONFLUENT_channel2_LOW_FOVs = Array.concat(CONFLUENT_channel2_LOW_FOVs, substring(lower_image,0,12)+"_ch2.tiff");
				//CONFLUENT_channel2_MID_FOVs = Array.concat(CONFLUENT_channel2_MID_FOVs, substring(middle_image,0,12)+"_ch2.tiff");
				print("ABOVE "+relative_coverage_THRESHOLD+"% COVERAGE");
				//
				run("Images to Stack", "name="+bottom_stack+" use keep");		
				selectWindow(focus_images[0]);
				close();
				//run("Images to Stack", "name="+lower_stack+" use keep");
				//selectWindow(focus_images[1]);
				//close();
				//run("Images to Stack", "name="+middle_stack+" use keep");
				//
				selectWindow(bottom_stack);
				run("Z Project...", "projection=[Max Intensity]");
				selectWindow("MAX_"+bottom_stack);
				saveAs("Tiff", output + "/BOTTOM/MAX_Z_CONFLUENT_CHANNEL2/"+substring(bottom_stack,0,16)+"_max");
				//
				//selectWindow(lower_stack);
				//run("Z Project...", "projection=[Max Intensity]");
				//selectWindow("MAX_"+lower_stack);
				//saveAs("Tiff", output + "\\LOWWER\\MAX_Z_CONFLUENT_CHANNEL2\\"+substring(lower_stack,0,16)+"_max");
				//
				//selectWindow(middle_stack);
				//run("Z Project...", "projection=[Max Intensity]");
				//selectWindow("MAX_"+middle_stack);
				//saveAs("Tiff", output + "\\MIDDLE\\MAX_Z_CONFLUENT_CHANNEL2\\"+substring(middle_stack,0,16)+"_max");
			}
			if(relative_coverage < 0 || relative_coverage > 100){ //This part of the if conditions is used if coverage is negative or above 100, in which case something has gone wrong. It terminates the script and raises an error
				exit("ERROR IN CALCULTING RELATIVE COVERAGE");
			}
			run("Close All"); //CLOSES ALL WINDOWS
			//
			for(i = 0; i <= 1; i++){
				focus_image_channel_1 = substring(focus_images[i],0,12)+"-ch4sk1fk1fl1.tiff";
				run("Bio-Formats Importer", "open=[" + currentDirectory + focus_image_channel_1 + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"); //IMPORTAN, WILL CRASH IF INDEX FILE IS IN DIRECTORY
			}
			focus_images_channel_1 = getList("image.titles");
			Array.sort(focus_images_channel_1);
			bottom_stack_channel_1 = substring(focus_images_channel_1[0],0,12) + "_ch1_STACK";
			//lower_stack_channel_1 = substring(focus_images_channel_1[1],0,12)+"_ch1_STACK";
			//middle_stack_channel_1 = substring(focus_images_channel_1[2],0,12)+"_ch1_STACK";
			//
			run("Images to Stack", "name="+bottom_stack_channel_1+" use keep");		
			selectWindow(focus_images_channel_1[0]);
			close();
			//run("Images to Stack", "name="+lower_stack_channel_1+" use keep");
			//selectWindow(focus_images_channel_1[1]);
			//close();
			//run("Images to Stack", "name="+middle_stack_channel_1+" use keep");
			//
			selectWindow(bottom_stack_channel_1);
			run("Z Project...", "projection=[Max Intensity]");
			//selectWindow(lower_stack_channel_1);
			//run("Z Project...", "projection=[Max Intensity]");
			//selectWindow(middle_stack_channel_1);
			//run("Z Project...", "projection=[Max Intensity]");
			
			if(relative_coverage >= relative_coverage_lower_THRESHOLD && relative_coverage < relative_coverage_THRESHOLD){
				selectWindow("MAX_"+bottom_stack_channel_1);
				saveAs("Tiff", output + "/BOTTOM/MAX_Z_CHANNEL1/"+substring(bottom_stack_channel_1,0,12)+"_ch1_max");
			//	selectWindow("MAX_"+lower_stack_channel_1);
			//	saveAs("Tiff", output + "\\LOWWER\\MAX_Z_CHANNEL1\\"+substring(lower_stack_channel_1,0,12)+"_ch1_max");
			//	selectWindow("MAX_"+middle_stack_channel_1);
			//	saveAs("Tiff", output + "\\MIDDLE\\MAX_Z_CHANNEL1\\"+substring(middle_stack_channel_1,0,12)+"_ch1_max");
			}if(relative_coverage >= relative_coverage_THRESHOLD && relative_coverage <= 100){
				selectWindow("MAX_"+bottom_stack_channel_1);
				saveAs("Tiff", output + "/BOTTOM/MAX_Z_CONFLUENT_CHANNEL1/"+substring(bottom_stack_channel_1,0,12)+"_ch1_max");
			//	selectWindow("MAX_"+lower_stack_channel_1);
			//	saveAs("Tiff", output + "\\LOWWER\\MAX_Z_CONFLUENT_CHANNEL1\\"+substring(lower_stack_channel_1,0,12)+"_ch1_max");
			//	selectWindow("MAX_"+middle_stack_channel_1);
			//	saveAs("Tiff", output + "\\MIDDLE\\MAX_Z_CONFLUENT_CHANNEL1\\"+substring(middle_stack_channel_1,0,12)+"_ch1_max");
			}
			}if(relative_coverage > 0 && relative_coverage <= relative_coverage_lower_THRESHOLD){
				EMPTY_FOVs = Array.concat(EMPTY_FOVs, substring(bottom_image,0,12)+"_ch2.tiff");
				print("EMPTY FOV: "+relative_coverage_THRESHOLD+"% COVERAGE");
			}
			run("Close All"); //CLOSES ALL WINDOWS
			Array.print(focus_images);
			run("Collect Garbage"); //SHOULD RE-CLAIM UNUSED RAM BUT MAY NOT WORK, WE SHALL SEE
		} //END OF IF CONDITION THAT IS TRIGGERED WHEN A NUMBER OF IMAGES CORRESPONDING TO NO. OF IMAGES_IN_STACK variable has been loaded in.
}



channel2_BOT_FOVs_STRING = String.join(channel2_BOT_FOVs, ",")
symbolic = output + "/BOTTOM/BOTTOM_FOVs.txt"
File.saveString(channel2_BOT_FOVs_STRING, symbolic)

//channel2_LOW_FOVs_STRING = String.join(channel2_LOW_FOVs, ",")
//symbolic = output + "\\LOWWER\\LOWWER_FOVs.txt"
//File.saveString(channel2_LOW_FOVs_STRING, symbolic)

//channel2_MID_FOVs_STRING = String.join(channel2_MID_FOVs, ",")
//symbolic = output + "\\MIDDLE\\MIDDLE_FOVs.txt"
//File.saveString(channel2_MID_FOVs_STRING, symbolic)


CONFLUENT_channel2_BOT_FOVs_STRING = String.join(CONFLUENT_channel2_BOT_FOVs, ",")
symbolic = output + "/BOTTOM/CONFLUENT_BOTTOM_FOVs.txt"
File.saveString(CONFLUENT_channel2_BOT_FOVs_STRING, symbolic)

//CONFLUENT_channel2_LOW_FOVs_STRING = String.join(CONFLUENT_channel2_LOW_FOVs, ",")
//symbolic = output + "\\LOWWER\\CONFLUENT_LOWWER_FOVs.txt"
//File.saveString(CONFLUENT_channel2_LOW_FOVs_STRING, symbolic)

//CONFLUENT_channel2_MID_FOVs_STRING = String.join(CONFLUENT_channel2_MID_FOVs, ",")
//symbolic = output + "\\MIDDLE\\CONFLUENT_MIDDLE_FOVs.txt"
//File.saveString(CONFLUENT_channel2_MID_FOVs_STRING, symbolic)

EMPTY_FOVs_STRING = String.join(EMPTY_FOVs, ",")
symbolic = output + "/EMPTY_FOVs.txt"
File.saveString(EMPTY_FOVs_STRING, symbolic)

end_time = getTime();
run_time = (end_time - start_time)/60000;
print("RUN TIME (MINUTES) WAS: "+run_time);
}


