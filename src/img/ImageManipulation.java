package img;

import static java.awt.image.BufferedImage.*;

import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferUShort;
import java.awt.image.PixelInterleavedSampleModel;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class ImageManipulation {
	
	/**
	 * Loads an image into an array of bytes. The first index is the color band,
	 * the second is the vertical coordinate, and the third is the horizontal one.
	 * @param img Path to the image
	 * @return Array containing the image
	 */
	public static byte[][][] loadImage(String img){
		try {
			return ImageManipulation.loadImage(
					ImageIO.read(new File(img)));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static short[][][] loadTiffImage(String img){
		try {
			return ImageManipulation.loadTiffImage(
					ImageIO.read(new File(img)));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static short[][][] loadTiffImage(BufferedImage image){

		switch (image.getType()) {
		case BufferedImage.TYPE_CUSTOM:

		}

		WritableRaster raster = image.getRaster();
		final short[] pixels = ((DataBufferUShort) raster.getDataBuffer()).getData();
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
		return result;
	}
	
	/**
	 * Loads an image into an array of bytes. The first index is the color band,
	 * the second is the vertical coordinate, and the third is the horizontal one.
	 * @param image bufferd image
	 * @return Array containing the image
	 */
	public static byte[][][] loadImage(BufferedImage image){
		  final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
	      final int width = image.getWidth();
	      final int height = image.getHeight();
	      final boolean hasAlphaChannel = image.getAlphaRaster() != null;
	      	      	      
	      switch(image.getType()){
	      case TYPE_3BYTE_BGR:{
//	    	  System.out.println("3BYTE_BGR");
	    	  byte[][][] result = new byte[4][height][width];
	    	  final int pixelLength = 3;
		         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
		            result[3][row][col] = (byte) 255; // alpha
		            result[2][row][col] = pixels[pixel]; // blue
		            result[1][row][col] = pixels[pixel + 1]; // green
		            result[0][row][col] = pixels[pixel + 2]; // red
		            col++;
		            if (col == width) {
		               col = 0;
		               row++;
		            }
		         }
		      return result;
	      }
	      case TYPE_4BYTE_ABGR:
	    	  System.out.println("4BYTE_ABGR");
	    	  break;
	      case TYPE_4BYTE_ABGR_PRE:
	    	  System.out.println("4BYTE_ABGR_PRE");
	    	  break; 
	      case TYPE_INT_RGB:
	    	  System.out.println("INT_RGB");
	    	  break;
	      case TYPE_INT_BGR:
	    	  System.out.println("INT_BGR");
	    	  break;
	      case TYPE_INT_ARGB:
	    	  System.out.println("INT_ARGB");
	    	  break;
	      case TYPE_INT_ARGB_PRE:
	    	  System.out.println("INT_ARGB_PRE");
	    	  break;
	      case 	TYPE_BYTE_BINARY:{
			System.out.println("TYPE_BYTE_BINARY");		
			byte[][][] result = new byte[1][height][width];
			for (int pixel = 0, row = 0, col = 0; pixel < pixels.length ; pixel++) {
				for(int i=0; i<8 && row<height; i++){
					byte c = (byte) (((pixels[pixel]>>(7-i)) & 0b01)*255);
					result[0][row][col] = c;
					col++;
					if (col == width) {
						col = 0;
						row++;
						i=0;
						pixel++;
					}
				}
			}
			return result;
	      }
	      case TYPE_BYTE_GRAY:{
//	    	  System.out.println("TYPE_BYTE_GRAY");
	    	  byte[][][] result = new byte[1][height][width];
	    	  final int pixelLength = 1;
	          for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
	             result[0][row][col] = pixels[pixel]; // gray
	             col++;
	             if (col == width) {
	                col = 0;
	                row++;
	             }
	          }
	    	  return result;
	      }
	    	  
	      }
	      

    	  System.out.println(image.getType());
	      
	      byte[][][] result = new byte[4][height][width];
	      if (hasAlphaChannel) {
	         final int pixelLength = 4;
	         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
	            result[3][row][col] = pixels[pixel]; // alpha
	            result[2][row][col] = pixels[pixel + 1]; // blue
	            result[1][row][col] = pixels[pixel + 2]; // green
	            result[0][row][col] = pixels[pixel + 3]; // red
	            col++;
	            if (col == width) {
	               col = 0;
	               row++;
	            }
	         }
	      } else {
	         final int pixelLength = 3;
	         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
	            result[3][row][col] = (byte) 255; // alpha
	            result[2][row][col] = pixels[pixel]; // blue
	            result[1][row][col] = pixels[pixel + 1]; // green
	            result[0][row][col] = pixels[pixel + 2]; // red
	            col++;
	            if (col == width) {
	               col = 0;
	               row++;
	            }
	         }
	      }

	      return result;
	}

	public static BufferedImage getBufferedImage(byte[][][] pixels){
		int height = pixels[0].length;
		int width = pixels[0][0].length;
		byte[] imgData = new byte[height*width*3];
		int index = 0;
		for(int h = 0; h < height; h++){
			for(int w = 0; w < width; w++){
//				imgData[index++] = pixels[3][h][w]; //alpha
				imgData[index++] = pixels[0][h][w]; //red
				imgData[index++] = pixels[1][h][w]; //green
				imgData[index++] = pixels[2][h][w]; //blue
			}
		}
	    BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
	    img2.setData(Raster.createRaster(
	    		new PixelInterleavedSampleModel(DataBuffer.TYPE_BYTE,width,height, 3, 3 * width, new int[]{0,1,2})
	    		, new DataBufferByte(imgData, imgData.length)
	    		, null ) );   
	    return img2;				
	}
	
	public static byte[][][] getRGBAArray(String img){
		try {
			return ImageManipulation.getRGBAArray(
					ImageIO.read(new File(img)));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static byte[][][] getRGBAArray(BufferedImage image){
		  final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
	      final int width = image.getWidth();
	      final int height = image.getHeight();
	      final boolean hasAlphaChannel = image.getAlphaRaster() != null;
	      
	      switch(image.getType()){
	      case TYPE_3BYTE_BGR:
	    	  System.out.println("3BYTE_BGR");
	    	  break;
	      case TYPE_4BYTE_ABGR:
	    	  System.out.println("4BYTE_ABGR");
	    	  break;
	      case TYPE_4BYTE_ABGR_PRE:
	    	  System.out.println("4BYTE_ABGR_PRE");
	    	  break; 
	      case TYPE_INT_RGB:
	    	  System.out.println("INT_RGB");
	    	  break;
	      case TYPE_INT_BGR:
	    	  System.out.println("INT_BGR");
	    	  break;
	      case TYPE_INT_ARGB:
	    	  System.out.println("INT_ARGB");
	    	  break;
	      case TYPE_INT_ARGB_PRE:
	    	  System.out.println("INT_ARGB_PRE");
	    	  break;
	     default:
	    	  System.out.println("INT_TYPE: "+image.getType());
	      }
	      
	      byte[][][] result = new byte[height][width][4];
	      if (hasAlphaChannel) {
	         final int pixelLength = 4;
	         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
	            result[row][col][3] = pixels[pixel]; // alpha
	            result[row][col][2] = pixels[pixel + 1]; // blue
	            result[row][col][1] = pixels[pixel + 2]; // green
	            result[row][col][0] = pixels[pixel + 3]; // red
	            col++;
	            if (col == width) {
	               col = 0;
	               row++;
	            }
	         }
	      } else {
	         final int pixelLength = 3;
	         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
	            result[row][col][3] = (byte) 255; // alpha
	            result[row][col][2] = pixels[pixel]; // blue
	            result[row][col][1] = pixels[pixel + 1]; // green
	            result[row][col][0] = pixels[pixel + 2]; // red
	            col++;
	            if (col == width) {
	               col = 0;
	               row++;
	            }
	         }
	      }

	      return result;
	}
	
	// red- band = 1; green- band = 2; blue- band = 3;
	public static byte[][] getRGBAArrayBand(BufferedImage image, int band){
//		for(String s: ImageIO.getReaderFileSuffixes())
//			System.out.println(s);
//		if(image==null)System.out.println("Cant read image");
		  final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
	      final int width = image.getWidth();
	      final int height = image.getHeight();
	      final boolean hasAlphaChannel = image.getAlphaRaster() != null;
	      final int pixelLength = image.getSampleModel().getNumBands();

	      byte[][] result = new byte[height][width];
	      if (hasAlphaChannel) {
	         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
	            result[row][col] = pixels[pixel + band];
	            col++;
	            if (col == width) {
	               col = 0;
	               row++;
	            }
	         }
	      } else {
	         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
	            result[row][col] = pixels[pixel + band - 1];
	            col++;
	            if (col == width) {
	               col = 0;
	               row++;
	            }
	         }
	      }

	      return result;
	}
	
	public static BufferedImage getRGBArrayBufferedImage(byte[][][] pixels){
		int height = pixels.length;
		int width = pixels[0].length;
		byte[] imgData = new byte[height*width*3];
		int index = 0;
		for(int h = 0; h < height; h++){
			for(int w = 0; w < width; w++){
//				imgData[index++] = pixels[h][w][3]; //alpha
				imgData[index++] = pixels[h][w][0]; //red
				imgData[index++] = pixels[h][w][1]; //green
				imgData[index++] = pixels[h][w][2]; //blue
			}
		}
	    BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
	    img2.setData(Raster.createRaster(
	    		new PixelInterleavedSampleModel(DataBuffer.TYPE_BYTE,width,height, 3, 3 * width, new int[]{0,1,2})
	    		, new DataBufferByte(imgData, imgData.length)
	    		, null ) );
	    
	    return img2;				
	}
		
	public static void saveRGBArrayAsImage(byte[][][] pixels, String path, String imgType){
		int height = pixels.length;
		int width = pixels[0].length;
		byte[] imgData = new byte[height*width*3];
		int index = 0;
		for(int h = 0; h < height; h++){
			for(int w = 0; w < width; w++){
//				imgData[index++] = pixels[h][w][3]; //alpha
				imgData[index++] = pixels[h][w][2]; //blue
				imgData[index++] = pixels[h][w][1]; //green
				imgData[index++] = pixels[h][w][0]; //red
			}
		}
	    BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
	    img2.setData(Raster.createRaster(
	    		new PixelInterleavedSampleModel(DataBuffer.TYPE_BYTE,width,height, 3, 3 * width, new int[]{0,1,2})
	    		, new DataBufferByte(imgData, imgData.length)
	    		, null ) );
	    
	    try {
			ImageIO.write(img2, imgType, new File(path));
		} catch (IOException e) {
			e.printStackTrace();
		}
				
	}
	
	public static boolean writeImage(byte[][][] pixels, String path){
		int height = pixels[0].length;
		int width = pixels[0][0].length;
		byte[] imgData = new byte[height*width*3];
		int index = 0;
		for(int h = 0; h < height; h++){
			for(int w = 0; w < width; w++){
				imgData[index++] = pixels[0][h][w]; //blue
				imgData[index++] = pixels[1][h][w]; //green
				imgData[index++] = pixels[2][h][w]; //red
			}
		}
	    BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
	    img2.setData(Raster.createRaster(
	    		new PixelInterleavedSampleModel(DataBuffer.TYPE_BYTE,width,height, 3, 3 * width, new int[]{0,1,2})
	    		, new DataBufferByte(imgData, imgData.length)
	    		, null ) );
	    String type = path.substring(path.lastIndexOf(".")+1, path.length());
	    try {
			return ImageIO.write(img2, type, new File(path));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return false;	
	}
	
	public static boolean writeImage(byte[][] pixels, String path){
		int height = pixels.length;
		int width = pixels[0].length;
		byte[] imgData = new byte[height*width*3];
		int index = 0;
		for(int h = 0; h < height; h++){
			for(int w = 0; w < width; w++){
				imgData[index++] = pixels[h][w]; //blue
			}
		}
	    BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
	    img2.setData(Raster.createRaster(
	    		new PixelInterleavedSampleModel(DataBuffer.TYPE_BYTE,width,height, 1, 1 * width, new int[]{0})
	    		, new DataBufferByte(imgData, imgData.length)
	    		, null ) );
	    String type = path.substring(path.lastIndexOf(".")+1, path.length());
	    try {
			return ImageIO.write(img2, type, new File(path));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return false;	
	}
	
	public static BufferedImage getGrayBufferedImage(byte[][] pixels){
		int height = pixels.length;
		int width = pixels[0].length;
		byte[] imgData = new byte[height*width];
		int index = 0;
		for(int h = 0; h < height; h++){
			for(int w = 0; w < width; w++){
				imgData[index++] = pixels[h][w];
			}
		}
	    BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
	    img2.setData(Raster.createRaster(
	    		new PixelInterleavedSampleModel(DataBuffer.TYPE_BYTE,width,height, 1, 1 * width, new int[]{0})
	    		, new DataBufferByte(imgData, imgData.length)
	    		, null ) );
	    
	    return img2;				
	}
	
/*	public static void saveGrayArrayAsImage(byte[][] pixels, String path, String imgType){
		int height = pixels.length;
		int width = pixels[0].length;
		byte[] imgData = new byte[height*width];
		int index = 0;
		for(int h = 0; h < height; h++){
			for(int w = 0; w < width; w++){
				imgData[index++] = pixels[h][w];
			}
		}
	    BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
	    img2.setData(Raster.createRaster(
	    		new PixelInterleavedSampleModel(DataBuffer.TYPE_BYTE,width,height, 1, 1 * width, new int[]{0})
	    		, new DataBufferByte(imgData, imgData.length)
	    		, null ) );
	    
	    try {
	    //Beging new Code
	    	JPEGImageWriteParam jpegParams = new JPEGImageWriteParam(null);
	    	jpegParams.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
	    	jpegParams.setCompressionQuality(1f);
	    	// use IIORegistry to get the available services
	    	IIORegistry registry = IIORegistry.getDefaultInstance();
	    	// return an iterator for the available ImageWriterSpi for jpeg images
	    	Iterator<ImageWriterSpi> services = registry.getServiceProviders(ImageWriterSpi.class,
	    	                                                 new ServiceRegistry.Filter() {   
	    	        @Override
	    	        public boolean filter(Object provider) {
	    	            if (!(provider instanceof ImageWriterSpi)) return false;

	    	            ImageWriterSpi writerSPI = (ImageWriterSpi) provider;
	    	            String[] formatNames = writerSPI.getFormatNames();
	    	            for (int i = 0; i < formatNames.length; i++) {
	    	                if (formatNames[i].equalsIgnoreCase("JPEG")) {
	    	                    return true;
	    	                }
	    	            }

	    	            return false;
	    	        }
	    	    },
	    	   true);
	    	//...assuming that servies.hasNext() == true, I get the first available service.
	    	ImageWriterSpi writerSpi = services.next();
	    	ImageWriter writer = writerSpi.createWriterInstance();

	    	// specifies where the jpg image has to be written
	    	writer.setOutput(new FileImageOutputStream(
	    	  new File(path)));

	    	// writes the file with given compression level 
	    	// from your JPEGImageWriteParam instance
	    	writer.write(null, new IIOImage(img2, null, null), jpegParams);
	    //End New Code
	    	
//			ImageIO.write(img2, imgType, new File(path));
		} catch (IOException e) {
			e.printStackTrace();
		}
				
	}*/
	
	public static byte[][] getGrayImage(BufferedImage image){
		final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
		final int width = image.getWidth();
		final int height = image.getHeight();

		switch (image.getType()) {
		case TYPE_3BYTE_BGR:
			System.out.println("3BYTE_BGR");
			break;
		case TYPE_4BYTE_ABGR:
			System.out.println("4BYTE_ABGR");
			break;
		case TYPE_4BYTE_ABGR_PRE:
			System.out.println("4BYTE_ABGR_PRE");
			break;
		case TYPE_INT_RGB:
			System.out.println("INT_RGB");
			break;
		case TYPE_INT_BGR:
			System.out.println("INT_BGR");
			break;
		case TYPE_INT_ARGB:
			System.out.println("INT_ARGB");
			break;
		case TYPE_INT_ARGB_PRE:
			System.out.println("INT_ARGB_PRE");
			break;

		}

		byte[][] result = new byte[height][width];
		final int pixelLength = 3;
		for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
			result[row][col] = (byte) ((pixels[pixel] & 0x0FF) / 3f
					+ (pixels[pixel + 1] & 0x0FF) / 3f + (pixels[pixel + 2] & 0x0FF) / 3f); // red
			col++;
			if (col == width) {
				col = 0;
				row++;
			}
		}

		return result;
	}
	
	// [x][y]
	public static byte[][] getGrayImg(BufferedImage image){
		final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
	    final int width = image.getWidth();
	    final int height = image.getHeight();
	    byte[][] result = new byte[width][height];
	    	    
	    final int pixelLength = 1;
        for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
           result[col][row] = pixels[pixel]; // gray
           col++;
           if (col == width) {
              col = 0;
              row++;
           }
        }
        return result;
	}

}
