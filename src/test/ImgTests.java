package test;

import image.tools.ITools;
import image.tools.IViewer;
import img.ImageManipulation;

public class ImgTests {
	
	public static void main(String[] args){
		String path = "C:/Users/Public/Pictures/Sample Pictures";
		path = "C:/Users/Harold/Pictures";
		byte[][][] img = ImageManipulation.loadImage(path+"/"+"desert.jpg");
		new IViewer("orig", ImageManipulation.getBufferedImage(img));
		double[][][] img2 = ITools.getNormColors(ITools.byte2double(img));
		new IViewer("norm", ImageManipulation.getBufferedImage(ITools.normilizeColor(new double[][][]{img2[0], img2[1], img2[2]})));
		new IViewer("gray", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[3])));
		new IViewer("red", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[0])));
		new IViewer("green", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[1])));
		new IViewer("blue", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[2])));
		
	}

}
