package img;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import image.tools.ITools;
import image.tools.IViewer;

/**
 * TERRAIN.BIL - Height field
 * TERRAIN.HDR - Height field header
 * TERRAIN.NRM - Normal's
 * TERRAIN.SBI - Optimized terrain layout info
 * TERRAIN.TID - Texture IDs
 * 
 * @author Harold
 *
 */
public class BIL {
	
	public static byte[][][] loadImage(String img){ //*.BIL
		byte[][][] out = new byte[1][1024][1024];
		try {
			InputStreamReader in = new InputStreamReader(new FileInputStream(img));
			int col = 0;
			int row = 0;
			while(in.ready()){
				int b = in.read() | in.read()<<8;
				out[0][row][col] = (byte) (b>>8);
				col++;
				if(col==1024){
					col=0;
					row++;
				}	
			}
			in.close();
			return out;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}	
		return null;
	}
	
	public static byte[][][] loadElevationImage(String img){ //*.BIL
		byte[][][] out = new byte[3][1024][1024];
		try {
			InputStreamReader in = new InputStreamReader(new FileInputStream(img));
			int col = 0;
			int row = 0;
			while(in.ready()){
				int b = in.read()|in.read()<<8;
//				java.awt.Color c = getElevationColor(b); // 65535
				java.awt.Color c = ITools.getNormalColorGrad(b/65535.0, 0.15);
				out[0][row][col] = (byte) c.getRed();
				out[1][row][col] = (byte) c.getGreen();
				out[2][row][col] = (byte) c.getBlue();
				col++;
				if(col==1024){
					col=0;
					row++;
				}	
			}
			in.close();
			return out;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}	
		return null;
	}
	
	public static byte[][][] loadNormalImage(String img){ //*.NRM
		byte[][][] out = new byte[2][1024][1024];
		try {
			InputStreamReader in = new InputStreamReader(new FileInputStream(img));
			int col = 0;
			int row = 0;
			while(in.ready()){
				out[0][row][col] = (byte) (in.read());
				out[1][row][col] = (byte) (in.read());
				col++;
				if(col==1024){
					col=0;
					row++;
				}	
			}
			in.close();
			return out;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}	
		return null;
	}
	
	public static int[][] loadTextureIDImage(String img){ //*.TID
		int[][] out = new int[1024][1024];
		try {
			InputStreamReader in = new InputStreamReader(new FileInputStream(img));
			int col = 0;
			int row = 0;
			while(in.ready()){
				int b = in.read()|in.read()<<8;
				out[row][col] = b;
				col++;
				if(col==1024){
					col=0;
					row++;
				}	
			}
			in.close();
			return out;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}	
		return null;
	}
	
	public static void main(String[] args){
		byte[][][] img = loadImage("F:/Games/Falcon BMS 4.33/Data/Terrdata/korea/terrain/TERRAIN_POLAK.BIL");
		byte[][][] nImg = loadNormalImage("F:/Games/Falcon BMS 4.33/Data/Terrdata/korea/terrain/TERRAIN_POLAK.NRM");
		byte[][][] elev = loadElevationImage("F:/Games/Falcon BMS 4.33/Data/Terrdata/korea/terrain/TERRAIN_POLAK.BIL");
//		new IViewer(ImageManipulation.getGrayBufferedImage(img[0]));
//		new IViewer("NRM", ImageManipulation.getGrayBufferedImage(nImg[0]));
//		new IViewer("Elev", ImageManipulation.getBufferedImage(elev));
		
//		ImageManipulation.saveGrayArrayAsImage(img[0], "F:/Games/Falcon BMS 4.33/Data/Terrdata/korea/terrain/TERRAIN_POLAK.png", "png");
//		ImageManipulation.saveGrayArrayAsImage(nImg[0], "F:/Games/Falcon BMS 4.33/Data/Terrdata/korea/terrain/TERRAIN_POLAK_NRM.png", "png");
//		ImageManipulation.writeImage(elev, "F:/Games/Falcon BMS 4.33/Data/Terrdata/korea/terrain/TERRAIN_POLAK_ELEV_MAP.png");
		
		byte[][][] cBar = new byte[3][800][800]; 
		for(int i=0; i<800; i++){
			java.awt.Color c = ITools.getNormalColorGrad(i/800.0, 0.15);
			for(int x=0; x<800; x++){
				cBar[0][i][x] = (byte) c.getRed();
				cBar[1][i][x] = (byte) c.getGreen();
				cBar[2][i][x] = (byte) c.getBlue();
			}
		}
		new IViewer("Colors", ImageManipulation.getBufferedImage(cBar));
	}
	
	public static java.awt.Color getElevationColor(int e){
		e /= 10;
		if(e<500){
			return new java.awt.Color( 169, 185, 210);
		}
		if(e<1000){
			return new java.awt.Color( 177, 188, 156);
		}
		if(e<2000){
			return new java.awt.Color( 142, 160, 120);
		}
		if(e<5000){
			return new java.awt.Color( 197, 188, 181);
		}
		if(e<7000){
			return new java.awt.Color( 212, 195, 141);
		}
		if(e<9000){
			return new java.awt.Color( 223, 223, 223);
		}
		System.out.println("Elevation too high: "+e);
		return new java.awt.Color( 0, 0, 0);
	}

}
