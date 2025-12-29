// Script to load the illumination corrected maxZ images of the cell body and nucleus channels, and then convert them from from 16-bit .tiffs into 8-bit RGBs. The RGB will have the nuclear signal in as the red (R) channel, and the cell body signal as the green (G) channel.
// The conversion also performs some scaling of the channels. Here, this scaling is based on manually inputted values that are applied to all images in the given experimental repeat. These values are determined by inspecting the images.
// The 8-bit RGBs are saved as .tifs and can be used for CellPose-based segmentation using the provied model (created by transfer learning) and for Ilastik-based cell classification (using the the provided model, which uses the 8-bit RGB in conjuction with the CellPose-based segmentation outputs in the form of .pngs)

run("Fresh Start"); //Fresh start closes all open images, clears ROI manager and results, it also sets backgrounds to black
output = getDirectory("Select the destination directory")

run("Bio-Formats Macro Extensions");
userChosenDirectory = getDirectory("Select the intial folder"); //SELECT A FOLDER ABOVE THE CHANNEL IMAGES SO (Select PROCESSED_IMAGES)
processBioFormatsFiles (userChosenDirectory);

function processBioFormatsFiles(currentDirectory) {
	function composite_rgb_maxz(POSITION){
		image_folder_channel1 = currentDirectory + POSITION + "/CORRECTED_MAX_Z_CHANNEL1/";
		image_folder_channel2 = currentDirectory + POSITION + "/CORRECTED_MAX_Z_CHANNEL2/";
		files_channel1 = getFileList(image_folder_channel1);
		files_channel2 = getFileList(image_folder_channel2);
		for(i = 0; i <= files_channel1.length - 1; i++){
			run("Bio-Formats Importer", "open=[" + image_folder_channel1 + files_channel1[i] + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"); //IMPORTAN, WILL CRASH IF INDEX FILE IS IN DIRECTORY
			run("Bio-Formats Importer", "open=[" + image_folder_channel2 + files_channel2[i] + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"); //IMPORTAN, WILL CRASH IF INDEX FILE IS IN DIRECTORY
			run("Merge Channels...", "c1="+files_channel1[i]+" c2="+files_channel2[i]+" create");
			Stack.setChannel(1);
			setMinAndMax(100, 5500);
			Stack.setChannel(2);
			//run("Enhance Contrast", "saturated=0.20");
			setMinAndMax(65, 1750);
			run("RGB Color");
			close("Composite");
			selectWindow("Composite (RGB)");
			wait(100); // Add a 100ms pause before saving
			image_name = substring(files_channel1[i],0,12)+"_ch"+substring(files_channel1[i],15,16)+"_ch"+substring(files_channel2[i],15,16)+"_RGB";
			rename(image_name);
			saveAs("Tiff", currentDirectory + POSITION + "/" + "RGB_CHANNEL1_2_DAPI_FIXED/" + image_name);
			wait(100); // Add a 100ms pause before saving
			run("Close All"); //CLOSES ALL WINDOWS
			run("Collect Garbage"); //SHOULD RE-CLAIM UNUSED RAM BUT MAY NOT WORK, WE SHALL SEE
		}
	}
	composite_rgb_maxz("BOTTOM");
}