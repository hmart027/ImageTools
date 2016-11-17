package img;

import image.tools.ITools;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.TreeMap;

public class JPEGEncoding {
	
	public static final int MAX_BUFFER_SIZE = 0xFFFFFF;
	public static byte[] buffer = new byte[MAX_BUFFER_SIZE];
	public static byte bitPointer = 0;
	public static int bytePointer = 0;

	private static final int[][] QUANTIZATION_MATRIX = new int[][]{
		new int[]{16,11,10,16,24,40,51,61},
		new int[]{12,12,14,19,26,58,60,55},
		new int[]{14,13,16,24,40,57,69,56},
		new int[]{14,17,22,29,51,87,80,62},
		new int[]{18,22,37,56,68,109,103,77},
		new int[]{24,35,55,64,81,104,113,92},
		new int[]{49,64,78,87,103,121,120,101},
		new int[]{72,92,95,98,112,100,103,99}
	};

	private static final int DC_CODES_COUNT 	= 16;
	private static final int[] DC_CODES 		= new int[DC_CODES_COUNT];
	private static final int[] DC_CODES_LENGTHS = new int[DC_CODES_COUNT];

	private static final int AC_CODES_COUNT 	= 256;
	private static final int[] AC_CODES 		= new int[AC_CODES_COUNT];
	private static final int[] AC_CODES_LENGTHS = new int[AC_CODES_COUNT];
	
	private static TreeMap<Integer, ArrayList<DCCode>> dcCodesByLenght = new TreeMap<>();
	private static TreeMap<Integer, ArrayList<ACCode>> acCodesByLenght = new TreeMap<>();
	
	private static final int[] CAT_MASKS = new int[]{
		0b000000000000000,
		0b000000000000001,
		0b000000000000011,
		0b000000000000111,
		0b000000000001111,
		0b000000000011111,
		0b000000000111111,
		0b000000001111111,
		0b000000011111111,
		0b000000111111111,
		0b000001111111111,
		0b000011111111111,
		0b000111111111111,
		0b001111111111111,
		0b011111111111111,
		0b111111111111111
	};
	
	static{
		loadDCCodes();
		loadACCodes();
//		for(int a: dcCodesByLenght.keySet()){
//			for(DCCode d: dcCodesByLenght.get(a))
//			System.out.println(a+"-"+d.cat+", "+d.val+" ,"+d.len);
//		}
	}
	
	public static boolean encodeImage(byte[][] img){
		byte[][] l_img = img;
		int h = l_img.length;
		int w = l_img[0].length;
		if(h==0 || w==0) return false;
		if(h!=w) l_img = ITools.zeroPadding(l_img);
		h = l_img.length;
		w = l_img[0].length;
		if((int)(Math.log(h)/Math.log(2))%1 != 0) return false;
		
		System.out.println(h+","+w);
		System.out.println(Math.log(h)/Math.log(2));
		
		int itterations = h/8;
		int yOffset = 0, xOffset = 0;
		byte[][] t_img = new byte[8][8];
		double[][] dct;
		int prevDC = 0;

		System.out.println("Starting to encode");
		for(int yO=0; yO<itterations; yO++){
			yOffset = yO*8;
			for (int xO=0; xO<itterations; xO++){
				xOffset = xO*8;
				for(int y=0; y<8; y++){
					for(int x=0; x<8; x++){
						t_img[y][x] = l_img[y+yOffset][x+xOffset];
					}
				}
//				dct = ImageTools.getFastCousine(t_img);
				dct = ITools.getDiscreteCousineTransform(t_img);
				prevDC = encodeBlock(dct, prevDC);
			}
		}		
		System.out.println("Last Byte Offset: "+bytePointer);
		return true;
	}
	
	public static int encodeBlock(double[][] dct, int prevDC){
		int[][] qdct = quantizeDCT(dct);
		int coef = 0, zeroRun = 0;
		int cDC = qdct[zigzagY[0]][zigzagX[0]]-prevDC;
		int[] encodedVal = encodeDC(cDC);
		bufferAppend(encodedVal[0], encodedVal[1]);		// Encode DC value
		for(int i = 1; i<64; i++){
			coef = qdct[zigzagY[i]][zigzagX[i]];
			if(coef == 0){
				zeroRun++;
			}else{
				encodedVal = encodeAC(zeroRun, coef);
				bufferAppend(encodedVal[0], encodedVal[1]);
				zeroRun = 0;
			}
			if(zeroRun>=16)break;
		}
		bufferAppend(0b01010, 4); 	// End of Block
		return cDC;
	}
	
	private static void bufferAppend(int val, int len){
		int rlen = len, needed = 0, shift = 0, nval =0;
		while(rlen > 0){
			needed = 8 - bitPointer;
			shift = rlen - needed;
			if(shift < 0){
				nval = val << -shift;
				buffer[bytePointer] |= nval;
				bitPointer += rlen;
				if(bitPointer >= 8){
					bitPointer = 0;
					bytePointer++;
				}
				rlen = 0;
			}else{
				nval = val >> shift;
				buffer[bytePointer] |= nval;
				bitPointer = 0;
				bytePointer++;
				rlen -= needed;
			}
		}
	}
	
	private static int getCategory(int i){
		if(i==0) return 0;
		return (int) (Math.log(Math.abs(i))/Math.log(2.0))+1;
	}
	
	private static int[] encodeDC(int dc){
		int cat = getCategory(dc);
		int val = dc;
		if(dc<0) val = (-val) ^ CAT_MASKS[cat];
		return new int[]{(DC_CODES[cat] << cat) | (val & CAT_MASKS[cat]), DC_CODES_LENGTHS[cat]};
	}

	private static int[] encodeAC(int zeroRun, int ac){
		int cat = getCategory(ac);
		int index = zeroRun * 16 + cat;
		int val = ac;
		if(ac<0) val = (-val) ^ CAT_MASKS[cat];
		return new int[]{(AC_CODES[index] << cat) | (val & CAT_MASKS[cat]), AC_CODES_LENGTHS[index]};
	}
	
	private static void loadACCodes(){
		try {
			java.io.BufferedReader in = new java.io.BufferedReader(new java.io.InputStreamReader(Class.forName("img.JPEGEncoding").getResourceAsStream("AC_Code.txt")));
			int index = 0;
			while(in.ready()){
				String line = in.readLine();
				String runCat = line.substring(0, line.indexOf(","));
				line = line.substring(line.indexOf(",")+1, line.length());
				String code = line.substring(0, line.indexOf(","));
				line = line.substring(line.indexOf(",")+1, line.length());
				index = getACIndex(runCat);
				AC_CODES[index]			= parseBitArray(code);
				AC_CODES_LENGTHS[index] = Integer.parseInt(line);
				ACCode ac = new ACCode((index>>4)&0x0F, index&0x0F, AC_CODES[index], AC_CODES_LENGTHS[index]);
				ArrayList<ACCode> acs = acCodesByLenght.get(ac.len - ac.cat);
				if(acs==null) acs = new ArrayList<>();
				acs.add(ac);
				acCodesByLenght.put(ac.len - ac.cat, acs);
			}
			in .close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	private static void loadDCCodes(){
		try {
			java.io.BufferedReader in = new java.io.BufferedReader(new java.io.InputStreamReader(Class.forName("img.JPEGEncoding").getResourceAsStream("DC_Code.txt")));
			int index = 0;
			while(in.ready()){
				String line = in.readLine();
				if(line.length()<=0) continue;
				String cat = line.substring(0, line.indexOf(","));
				line = line.substring(line.indexOf(",")+1, line.length());
				String code = line.substring(0, line.indexOf(","));
				line = line.substring(line.indexOf(",")+1, line.length());
				index = new BigInteger(cat, 16).intValue();
				DC_CODES[index]			= parseBitArray(code);
				DC_CODES_LENGTHS[index] = Integer.parseInt(line);
				DCCode dc = new DCCode(index, DC_CODES[index], DC_CODES_LENGTHS[index]);
				ArrayList<DCCode> dcs = dcCodesByLenght.get(dc.len - dc.cat);
				if(dcs==null) dcs = new ArrayList<>();
				dcs.add(dc);
				dcCodesByLenght.put(dc.len - dc.cat, dcs);
			}
			in .close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	private static int getACIndex(String runCat){
		int out = 0;
		out =  new BigInteger(runCat.substring(0, runCat.indexOf('/')), 16).intValue() * 0x10;
		out += new BigInteger(runCat.substring(runCat.indexOf('/')+1, runCat.length()), 16).intValue();
		return out;
	}
	
	private static int parseBitArray(String bitArray){
		int out = 0;
		for(char c: bitArray.toCharArray()){
			out <<= 1;
			if(c == '1') out |= 1;
		}			
		return out;
	}

	private static int[][] quantizeDCT(double[][] dct){
		int[][] out = new int[8][8];
		int t = 0, r = 0;
		for(int y = 0; y<8; y++){
			for(int x = 0; x<8; x++){
				t = ((int)dct[y][x]/QUANTIZATION_MATRIX[y][x]);
				r = (((int)dct[y][x]%QUANTIZATION_MATRIX[y][x])*10);
				if(r>=  5*QUANTIZATION_MATRIX[y][x]) t++;
				if(r<= -5*QUANTIZATION_MATRIX[y][x]) t--;
				out[y][x] = t;
			}
		}
		return out;
	}

	private static final int[] zigzagX = new int[]{
		0,1,0,
		0,1,2,3,2,1,0,
		0,1,2,3,4,5,4,3,2,1,0,
		0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,
		1,2,3,4,5,6,7,7,6,5,4,3,2,
		3,4,5,6,7,7,6,5,4,
		5,6,7,7,6,7
	};
	private static final int[] zigzagY = new int[]{
		0,0,1,
		2,1,0,0,1,2,3,
		4,3,2,1,0,0,1,2,3,4,5,
		6,5,4,3,2,1,0,0,1,2,3,4,5,6,7,
		7,6,5,4,3,2,1,2,3,4,5,6,7,
		7,6,5,4,3,4,5,6,7,
		7,6,5,6,7,
		7
	};

	private static double[][] img;
	
	public static double[][] decodeImage(byte[] img, int start, int bitLenght, int N){
		int byteIntex = 0, bitIndex = 0, totalBitIndex = 0;
		int xO = 0, yO = 0;
		double[][] out = new double[N][N];
		while(totalBitIndex < bitLenght){
			int cInt = loadInt(img, byteIntex, bitIndex);
			DCCode dcCode = findDCCode(cInt);
			
			ACCode lAC = null;
			do{
				
			}while(lAC.ac != 0b1010);//End of Block found
		}
		return out;
	}
	
	public static int[] decodeBlock(byte[] img, int byteI, int bit, int xO, int yO){
		DCCode dcCode;
		ACCode acCode;
		ArrayList<ACCode> acCodes = new ArrayList<>();
		int cInt = loadInt(img, byteI, bit);
		dcCode = findDCCode(cInt);
		int l = dcCode.len;
		if(l>(8-bit)){
			int remB = l-(8-bit);
			byteI += 1+remB/8;
			bit = remB%8;
		}else{
			bit += l;
		}
		do{
			cInt = loadInt(img, byteI, bit);
			acCode = findACCode(cInt);
			acCodes.add(acCode);
			l = acCode.len;
			if(l>(8-bit)){
				int remB = l-(8-bit);
				byteI += 1+remB/8;
				bit = remB%8;
			}else{
				bit += l;
			}
		}while(acCode.ac != 0b1010);//End of Block found
		
		System.out.println("DC: "+dcCode.dc);
		for(ACCode ac: acCodes)
			System.out.println("AC: "+ac.zero+","+ac.ac);
		
		return new int[]{byteI, bit};
	}
	
	public static DCCode findDCCode(int b){
		DCCode code=null;
		int shift = 31;
		int mask = 0;
		int maskL = 1;
		int val = 0;
		while(shift>0){
			maskL = 32-shift;
			mask = (b>>shift) & getInclusiveMask(maskL);
			ArrayList<DCCode> pDC = dcCodesByLenght.get(maskL);
			if(pDC != null && pDC.size()>0){
				for(DCCode dc: pDC){
					if(mask == dc.val){
						code = new DCCode(dc.cat, dc.val, dc.len);
						val = (b>>(shift-dc.cat)) & getInclusiveMask(dc.cat); 
						System.out.println(mask);
						code.dc = val;
						return code;
					}
				}
			}
			shift--;
		}
		return code;
	}
	
	public static ACCode findACCode(int b){
		ACCode code=null;
		int shift = 31;
		int mask = 0;
		int maskL = 1;
		int val = 0;
		while(shift>0){
			maskL = 32-shift;
			mask = (b>>shift) & getInclusiveMask(maskL);
			ArrayList<ACCode> pAC = acCodesByLenght.get(maskL);
			if(pAC != null && pAC.size()>0){
				for(ACCode ac: pAC){
					if(mask == ac.val){
						code = new ACCode(ac.zero, ac.cat, ac.val, ac.len);
						val = (b>>(shift-ac.cat)) & getInclusiveMask(ac.cat); 
						System.out.println(mask);
						code.ac = val;
						return code;
					}
				}
			}
			shift--;
		}
		return code;
	}
	
	public static int loadInt(byte[] array, int byteIndex, int bitIndexinByte){
		int out = 0;
		int index = 0;
		int shift = 31;
		while(shift>0 && index<array.length){
			if(shift<=8){
				out |= (array[index]>>(8-shift)) & getInclusiveMask(shift);
				break;
			}
			if(bitIndexinByte==0){
				out |= (array[index] & 0xFF)<<(shift-7);
				shift-=8;
			}
			if(bitIndexinByte>0){
				int l = 8-bitIndexinByte;
				out |= (array[index]&getInclusiveMask(l))<<(shift-l+1);
				shift-=l;
				bitIndexinByte = 0;
			}
			index++;
		}
		return out;
	}
	
	private static int getInclusiveMask(int n){
		return (int) (Math.pow(2, n)-1);
	}
	
	public static class DCCode{
		public int cat, val, len, dc;
		public DCCode(int cat, int val, int len){
			this.cat = cat;
			this.val = val;
			this.len = len;
		}
		public boolean equals(Object o){
			if(o==null)return false;
			if(!o.getClass().equals(this.getClass())) return false;
			return cat==((DCCode)o).cat;
		}
	}
	
	public static class ACCode{
		public int zero, cat, val, len, ac;
		public ACCode(int zero, int cat, int val, int len){
			this.zero = zero;
			this.cat = cat;
			this.val = val;
			this.len = len;
		}
		public boolean equals(Object o){
			if(o==null)return false;
			if(!o.getClass().equals(this.getClass())) return false;
			return (zero==((ACCode)o).zero) && (cat==((ACCode)o).cat);
		}
	}
}
