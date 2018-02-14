package test;

import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferShort;
import java.awt.image.DataBufferUShort;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import image.tools.ITools;
import image.tools.IViewer;
import img.ImageManipulation;

public class ImgTests {
	
	public static void main(String[] args){
//		String path = "C:/Users/Public/Pictures/Sample Pictures";
//		path = "C:/Users/Harold/Pictures";
//		byte[][][] img = ImageManipulation.loadImage(path+"/"+"desert.jpg");
//		new IViewer("orig", ImageManipulation.getBufferedImage(img));
//		double[][][] img2 = ITools.getNormColors(ITools.byte2double(img));
//		new IViewer("norm", ImageManipulation.getBufferedImage(ITools.normilizeColor(new double[][][]{img2[0], img2[1], img2[2]})));
//		new IViewer("gray", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[3])));
//		new IViewer("red", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[0])));
//		new IViewer("green", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[1])));
//		new IViewer("blue", ImageManipulation.getGrayBufferedImage(ITools.normilize(img2[2])));
		
//		for(int i=0; i<180; i+=20){
//			new IViewer(i+"",ImageManipulation.getGrayBufferedImage(ITools.normilize(ITools.getEdgeKernel(i))));
//		}
		
		try {
			BufferedImage image = ImageIO.read(new File("/home/harold/14NOV27004452-M2AS-054191978040_01_P001.TIF"));
			short[] pixels = ((DataBufferUShort)image.getRaster().getDataBuffer()).getData();
			System.out.println(image.getRaster().getNumBands());
			System.out.println("Img loaded: "+pixels);
			
			final int width = image.getWidth();
			final int height = image.getHeight();

			final int pixelLength = image.getRaster().getNumBands();
			short[][][] result = new short[pixelLength][height][width];
			int pixel = 0;
			for (int row = 0; row < height; row++) {
				for (int col = 0; col < width; col++) {
					for (int band = 0; band < pixelLength; band++) {
						result[band][row][col] = pixels[pixel++];
					}
				}
			}

			short max = result[0][0][0], min = max;
			for (int band = 0; band < pixelLength; band++) {
			    for (int row = 0; row < height; row++) {
					for (int col = 0; col < width; col++) {
						if(result[band][row][col]>max)
							max = result[band][row][col];
						if(result[band][row][col]<min)
							min = result[band][row][col];
					}
				}
			}
			System.out.println("max: "+max);
			System.out.println("min: "+min);
			
			new IViewer(image);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

}
