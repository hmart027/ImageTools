package image.tools;

import math2.Complex;

public class ITools {
	
//	public static byte[][] crop(byte[][] img, int w, int h){
//		int h0 = img.length;
//		int w0 = img[0].length;
//		
//		return null;
//	}
	
	public static boolean CENTERED_FFT = true;
	
	public static byte[][] setAspectOneToOne(byte[][] img){
		int h0 = img.length;
		int w0 = img[0].length;
		int ht = h0;	// Target Height
		int wt = w0;	// Target Width
		
		if(wt>ht) ht = wt;
		
		return null;
	}

	public static byte[][] shrinkByTwo(byte[][] img){
		int h0 = img.length;
		int h1 = h0/2;
		
		byte[][] out = new byte[h1][h1];
		
		for(int y = 0; y<h1; y++){
			for(int x = 0; x<h1; x++){
				out[y][x] = (byte)((img[y*2][x*2] & 0x0FF)/4f + (img[y*2][x*2+1] & 0x0FF)/4f + (img[y*2+1][x*2] & 0x0FF)/4f + (img[y*2+1][x*2+1] & 0x0FF)/4f);
			}
		}
		
		return out;
	}

	public static Complex[][] getDiscreteFourier(byte[][] img){
		int h0 = img.length;
		
		Complex[][] out = new Complex[h0][h0];
		
		for(int v = 0; v<h0; v++){
			for(int u = 0; u<h0; u++){
				out[v][u] = new Complex(0, 0);
				for(int y = 0; y<h0; y++){
					for(int x = 0; x<h0; x++){
						Complex r = new Complex((img[y][x]&0x0FF)*Math.pow(-1, x+y), 0);
						Complex i = new Complex(0, -2*Math.PI*(u*x+v*y)/(double)h0);
						out[v][u] = out[v][u].plus(r.times(i.exp()));
					}
				}
				out[v][u] = out[v][u].div(new Complex(Math.pow(h0, 2),0));
			}
		}
		
		return out;
	}
	
	public static double[][] getInverseDiscreteFourier(Complex[][] img){
		int h0 = img.length;
		
		double[][] out = new double[h0][h0];
		Complex t = null;
		
		for(int y = 0; y<h0; y++){
			for(int x = 0; x<h0; x++){
				out[y][x] = 0;
				t = new Complex(0,0);
				for(int v = 0; v<h0; v++){
					for(int u = 0; u<h0; u++){
						Complex r = new Complex(Math.pow(-1, x+y), 0);
						Complex i = new Complex(0, 2*Math.PI*(u*x+v*y)/(double)h0);
						t=t.plus(img[v][u].times(r).times(i.exp()));
					}
				}
				out[y][x] = t.real();
			}
		}
		
		return out;
	}
	
	public static double[][] getDiscreteFourierMod(Complex[][] fourier){
		int h0 = fourier.length;		
		double[][] out = new double[h0][h0];	
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < h0; x++) {
				out[y][x] = Math.pow(Math.pow(fourier[y][x].real(),2) + Math.pow(fourier[y][x].imag(),2), 0.5);
			}
		}		
		return out;
	}
	
	public static double[][] getReal(Complex[][] fourier){
		int h0 = fourier.length;	
		double[][] out = new double[h0][h0];	
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < h0; x++) {
				out[y][x] = fourier[y][x].real();
			}
		}
		return out;
	}
	
	public static double[][] getImaginary(Complex[][] fourier){
		int h0 = fourier.length;	
		double[][] out = new double[h0][h0];	
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < h0; x++) {
				out[y][x] = fourier[y][x].imag();
			}
		}
		return out;
	}
	
	public static double[][] getDiscreteCousineTransform(byte[][] img){
		int h0 = img.length;
		
		double[][] out = new double[h0][h0];
		double t1 = Math.pow(2d/(double)h0, 0.5);
		double t2 = 1d/Math.pow(2.0, 0.5);
		double c = Math.PI/(2*(double)h0);
		
		for(int v = 0; v<h0; v++){
			for(int u = 0; u<h0; u++){
				out[v][u] = 0;
				for(int y = 0; y<h0; y++){
					for(int x = 0; x<h0; x++){
						out[v][u] += t1*(img[y][x]&0x0FF)*Math.cos((2*x+1)*u*c)*Math.cos((2*y+1)*v*c);
					}
				}
				if(v==0) out[v][u] *= t2;
				if(u==0) out[v][u] *= t2;
			}
		}
		return out;
	}
	
	public static double[][] getInverseDiscreteCousineTransform(double[][] DCT){
		int h0 = DCT.length;
		
		double[][] out = new double[h0][h0];
		double t1 = Math.pow(2d/(double)h0, 0.5);
		double t2 = 1d/Math.pow(2.0, 0.5);
		double c = Math.PI/(2*(double)h0);
		double temp = 0;
		
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < h0; x++) {
				out[y][x] = 0;
				temp = 0;
				for (int v = 0; v < h0; v++) {
					for (int u = 0; u < h0; u++) {
						temp = t1*DCT[v][u]*Math.cos((2*x+1)*u*c)*Math.cos((2*y+1)*v*c);
						if(u==0) temp *= t2;
						if(v==0) temp *= t2;
						out[y][x] += temp;
					}
				}
			}
		}		
		return out;
	}
	
	public static byte[][] normilize(double[][] img){
		int h0 = img.length;
		int w0 = img[0].length;
		
		byte[][] out = new byte[h0][w0];
		double max = 0, min = 0, range = 0;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				if(img[y][x] > max) max = img[y][x];
				if(img[y][x] < min) min = img[y][x];  
			}
		}
		range = max-min;
//		System.out.println("Max: "+max+", Min: "+min);
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				out[y][x] = (byte) ((img[y][x]-min)/range*255); 
			}
		}
		return out;
	}
	
	public static byte[][] normilize(float[][] img){
		int h0 = img.length;
		int w0 = img[0].length;
		
		byte[][] out = new byte[h0][w0];
		double max = 0, min = 0, range = 0;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				if(img[y][x] > max) max = img[y][x];
				if(img[y][x] < min) min = img[y][x];  
			}
		}
		range = max-min;
//		System.out.println("Max: "+max+", Min: "+min);
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				out[y][x] = (byte) ((img[y][x]-min)/range*255); 
			}
		}
		return out;
	}
	
	public static byte[][][] normilize(float[][][] img){
		int c0 = img.length;
		int h0 = img[0].length;
		int w0 = img[0][0].length;
		
		byte[][][] out = new byte[c0][h0][w0];
		double max = 0, min = 0, range = 0;

		for (int c = 0; c < c0; c++) {
			for (int y = 0; y < h0; y++) {
				for (int x = 0; x < w0; x++) {
					if(img[c][y][x] > max) max = img[c][y][x];
					if(img[c][y][x] < min) min = img[c][y][x];  
				}
			}
		}
		range = max-min;
//		System.out.println("Max: "+max+", Min: "+min);
		for (int c = 0; c < c0; c++) {
			for (int y = 0; y < h0; y++) {
				for (int x = 0; x < w0; x++) {
					out[c][y][x] = (byte) ((img[c][y][x]-min)/range*255); 
				}
			}
		}
		return out;
	}
	
	/**
	 * Normalizes a color image. The color bvand index is expected to be the first one.
	 * @param img
	 * @return
	 */
	public static byte[][][] normilizeColor(double[][][] img){
		int c0 = img.length;
		int h0 = img[0].length;
		int w0 = img[0][0].length;
		
		byte[][][] out = new byte[c0][h0][w0];
		double max = 0, min = 0, range = 0;
		for (int c = 0; c < c0; c++) {
			for (int y = 0; y < h0; y++) {
				for (int x = 0; x < w0; x++) {
					if(img[c][y][x] > max) max = img[c][y][x];
					if(img[c][y][x] < min) min = img[c][y][x];  
				}
			}
			range = max-min;
//			System.out.println("Max "+c+": "+max[c]+", Min "+c+": "+min[c]);
		}
		for (int c = 0; c < c0; c++) {
			for (int y = 0; y < h0; y++) {
				for (int x = 0; x < w0; x++) {
					out[c][y][x] = (byte) ((img[c][y][x]-min)/range*255); 
				}
			}
		}
		return out;
	}
	
	public static byte[][] normilizeKeepZero(double[][] img){
		int h0 = img.length;
		int w0 = img[0].length;
		
		byte[][] out = new byte[h0][w0];
		double max = 0, min = 0, range = 0;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				if(img[y][x] > max) max = img[y][x];
				if(img[y][x] < min) min = img[y][x];  
			}
		}
		range = Math.abs(max);
//		System.out.println("Max: "+max+", Min: "+min);
		if(Math.abs(min)>range) range = Math.abs(min);
		range*=2;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				out[y][x] = (byte) ((img[y][x]-min)/range*255); 
			}
		}
		return out;
	}
	
	public static byte[][] normilize(byte[][] img1){
		int[] s = {img1.length, img1[0].length};
		int tval;
		
		byte[][] out = new byte[s[0]][s[1]];
		double max = 0, min = 0, range = 0;
		for (int y = 0; y < s[0]; y++) {
			for (int x = 0; x < s[1]; x++) {
				tval = (img1[y][x] & 0x0FF);
				if(tval > max) max = tval;
				if(tval < min) min = tval;  
			}
		}
//		System.out.println(min+", "+max);
		range = max-min;
		for (int y = 0; y < s[0]; y++) {
			for (int x = 0; x < s[1]; x++) {
				tval = (img1[y][x] & 0x0FF);
				out[y][x] = (byte) ((tval-min)/range*255); 
			}
		}
		return out;
	}
	
	public static byte[][] logScale(double[][] fourier){
		int h0 = fourier.length;
		byte[][] out = new byte[h0][h0];
		double max = 0;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < h0; x++) {
				double t = Math.abs(fourier[y][x]);
				if(t > max) max = t;
			}
		}
		double scale1 = 255.0/Math.log10(256.0);
		double scale2 = 255.0/max;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < h0; x++) {
				double t = Math.abs(fourier[y][x]);
				out[y][x] = (byte) (scale1*Math.log10(1.0 + scale2*t)); 
			}
		}
		return out;
	}
	
	public static double[][] getWalshTransform(byte[][] img){
		int N = img.length;
		int N2 = (int) (Math.log(N)/Math.log(2));
		if(N != img[0].length) return null;
		
		double[][] out    = new double[N][N];
		for(int v = 0; v<N; v++){
			for(int u = 0; u<N; u++){
				out[v][u] = 0;
				for(int y = 0; y<N; y++){
					for(int x = 0; x<N; x++){
						double p = 1;
						for(int i = 0; i<N2; i++){
							float e = (float)((x>>i) & 0x01) * (float)((u>>(N2-1-i)) & 0x01) + 
									(float)((y>>i) & 0x01) * (float)((v>>(N2-1-i)) & 0x01);
							p *= Math.pow(-1.0, e);
						}
						out[v][u] += (double)(img[y][x] & 0x0FF)*p;
					}
				}
				out[v][u] *= Math.pow(N,-2);
			}
		}
		
		return out;
	}
	
	public static double[][] getReverseWalshTransform(double[][] img){
		int N = img.length;
		int N2 = (int) (Math.log(N)/Math.log(2));
		if(N != img[0].length) return null;
		
		double[][] out    = new double[N][N];
		for(int v = 0; v<N; v++){
			for(int u = 0; u<N; u++){
				out[v][u] = 0;
				for(int y = 0; y<N; y++){
					for(int x = 0; x<N; x++){
						double p = 1;
						for(int i = 0; i<N2; i++){
							float e = (float)((x>>i) & 0x01) * (float)((u>>(N2-1-i)) & 0x01) + 
									(float)((y>>i) & 0x01) * (float)((v>>(N2-1-i)) & 0x01);
							p *= Math.pow(-1.0, e);
						}
						out[v][u] += img[y][x]*p;
					}
				}
			}
		}
		
		return out;
	}
	
	/**
	 * Fast Foutier Transform (FFT)
	 */
	public static Complex[][] getFastFourier(byte[][] img){
		int N = img.length;
		int nBits = (int) (Math.log(N)/Math.log(2));
		//Generate W array
		Complex[] w = new Complex[N];
		for(int v = 0; v<N; v++){
			w[v] = new Complex(0, -2.0 * Math.PI / (double) N * (double)v).exp();
		}
		return FFT(1,N,nBits,img,w);
	}

	private static Complex[][] FFT(int level, int oSize, int nBits, byte[][] img, Complex[] w){
		int levelSize = oSize/level;
		//At lowest possible level
		if(levelSize<=2){
			int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
			double c = Math.pow(oSize, -2);
			double exp = -2.0*Math.PI/(double)levelSize;
			Complex[][] fft = new Complex[oSize][oSize];
			Complex F00, F01, F10, F11;
			for(int v = 0; v<fft.length; v+=2){
				for(int u = 0; u<fft.length; u+=2){
					x0 = bitReversal(u, nBits);
					x1 = bitReversal(u+1, nBits);
					y0 = bitReversal(v, nBits);
					y1 = bitReversal(v+1, nBits);
					if(CENTERED_FFT){
						F00 = new Complex((img[y0][x0] & 0x0FF)*Math.pow(-1, y0+x0), 0);
						F01 = new Complex((img[y0][x1] & 0x0FF)*Math.pow(-1, y0+x1), 0);
						F10 = new Complex((img[y1][x0] & 0x0FF)*Math.pow(-1, y1+x0), 0);
						F11 = new Complex((img[y1][x1] & 0x0FF)*Math.pow(-1, y1+x1), 0);
					}else{
						F00 = new Complex((img[y0][x0] & 0x0FF), 0);
						F01 = new Complex((img[y0][x1] & 0x0FF), 0);
						F10 = new Complex((img[y1][x0] & 0x0FF), 0);
						F11 = new Complex((img[y1][x1] & 0x0FF), 0);
					}
					x0=0;x1=1;y0=0;y1=1;
					// u=0;v=0
					fft[v][u] = (new Complex(c, 0)).times(F00.plus(F01).plus(F10).plus(F11));
					// u=1;v=0
					fft[v][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F01.times(new Complex(0,exp*(double)(x1)).exp())).plus //x1
							(F10.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F11.times(new Complex(0,exp*(double)(x1)).exp())) //x1
							);
					// u=0;v=1
					fft[v+1][u] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F01.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F10.times(new Complex(0,exp*(double)(y1)).exp())).plus //y1
							(F11.times(new Complex(0,exp*(double)(y1)).exp()))
							);
					// u=1;v=1
					fft[v+1][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0+y0)).exp())).plus
							(F01.times(new Complex(0,exp*(double)(x1+y0)).exp())).plus
							(F10.times(new Complex(0,exp*(double)(x0+y1)).exp())).plus
							(F11.times(new Complex(0,exp*(double)(x1+y1)).exp()))
							);
				}
			}
			return fft;
		}
		Complex[][] fftLower = FFT(level*2, oSize, nBits, img,w);//Lower Level FFT
		Complex[][] fft = new Complex[oSize][oSize];//FFT to return
		int rad = levelSize/2;
		for(int vOffset = 0; vOffset<oSize; vOffset+=levelSize){
			for(int uOffset = 0; uOffset<oSize; uOffset+=levelSize){
				for(int v = 0; v<rad; v++){
					for(int u = 0; u<rad; u++){
						int cv = vOffset+ v;
						int cu = uOffset+ u;
						fft[cv][cu] 			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv][cu+rad] 		= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu]			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu+rad] 	= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));	
					}
				}
			}
		}
		return fft;
	}

	/**
	 * Color Fast Foutier Transform (FFT)
	 */
	public static Complex[][][] getColorFastFourier(byte[][][] img){
		int N = img.length;
		int nBits = (int) (Math.log(N)/Math.log(2));
		//Generate W array
		Complex[] w = new Complex[N];
		for(int v = 0; v<N; v++){
			w[v] = new Complex(0, -2.0 * Math.PI / (double) N * (double)v).exp();
		}
		return colorFFT(1,N,nBits,img,w);
	}

	private static Complex[][][] colorFFT(int level, int oSize, int nBits, byte[][][] img, Complex[] w){
		int levelSize = oSize/level;
		//At lowest possible level
		if(levelSize<=2){
			int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
			double c = Math.pow(oSize, -2);
			double exp = -2.0*Math.PI/(double)levelSize;
			Complex[][][] fft = new Complex[oSize][oSize][3];
			Complex F00, F01, F10, F11;
			for(int v = 0; v<fft.length; v+=2){
				for(int u = 0; u<fft.length; u+=2){
					x0 = bitReversal(u, nBits);
					x1 = bitReversal(u+1, nBits);
					y0 = bitReversal(v, nBits);
					y1 = bitReversal(v+1, nBits);
					for(int co=0; co<3; co++){
						if(CENTERED_FFT){
							F00 = new Complex((img[y0][x0][co] & 0x0FF)*Math.pow(-1, y0+x0), 0);
							F01 = new Complex((img[y0][x1][co] & 0x0FF)*Math.pow(-1, y0+x1), 0);
							F10 = new Complex((img[y1][x0][co] & 0x0FF)*Math.pow(-1, y1+x0), 0);
							F11 = new Complex((img[y1][x1][co] & 0x0FF)*Math.pow(-1, y1+x1), 0);
						}else{
							F00 = new Complex((img[y0][x0][co] & 0x0FF), 0);
							F01 = new Complex((img[y0][x1][co] & 0x0FF), 0);
							F10 = new Complex((img[y1][x0][co] & 0x0FF), 0);
							F11 = new Complex((img[y1][x1][co] & 0x0FF), 0);
						}
						x0=0;x1=1;y0=0;y1=1;
						// u=0;v=0
						fft[v][u][co] = (new Complex(c, 0)).times(F00.plus(F01).plus(F10).plus(F11));
						// u=1;v=0
						fft[v][u+1][co] = (new Complex(c, 0)).times(
								(F00.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
								(F01.times(new Complex(0,exp*(double)(x1)).exp())).plus //x1
								(F10.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
								(F11.times(new Complex(0,exp*(double)(x1)).exp())) //x1
								);
						// u=0;v=1
						fft[v+1][u][co] = (new Complex(c, 0)).times(
								(F00.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
								(F01.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
								(F10.times(new Complex(0,exp*(double)(y1)).exp())).plus //y1
								(F11.times(new Complex(0,exp*(double)(y1)).exp()))
								);
						// u=1;v=1
						fft[v+1][u+1][co] = (new Complex(c, 0)).times(
								(F00.times(new Complex(0,exp*(double)(x0+y0)).exp())).plus
								(F01.times(new Complex(0,exp*(double)(x1+y0)).exp())).plus
								(F10.times(new Complex(0,exp*(double)(x0+y1)).exp())).plus
								(F11.times(new Complex(0,exp*(double)(x1+y1)).exp()))
								);
					}
				}
			}
			return fft;
		}
		Complex[][][] fftLower = colorFFT(level*2, oSize, nBits, img,w);//Lower Level FFT
		Complex[][][] fft = new Complex[oSize][oSize][3];//FFT to return
		int rad = levelSize/2;
		for(int vOffset = 0; vOffset<oSize; vOffset+=levelSize){
			for(int uOffset = 0; uOffset<oSize; uOffset+=levelSize){
				for(int v = 0; v<rad; v++){
					for(int u = 0; u<rad; u++){
						int cv = vOffset+ v;
						int cu = uOffset+ u;
						for(int c=0; c<3; c++){
							fft[cv][cu][c] 			= fftLower[cv][cu][c].plus(fftLower[cv][cu+rad][c].times(w[u*level])).plus(fftLower[cv+rad][cu][c].times(w[v*level])).plus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
							fft[cv][cu+rad][c] 		= fftLower[cv][cu][c].minus(fftLower[cv][cu+rad][c].times(w[u*level])).plus(fftLower[cv+rad][cu][c].times(w[v*level])).minus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
							fft[cv+rad][cu][c]		= fftLower[cv][cu][c].plus(fftLower[cv][cu+rad][c].times(w[u*level])).minus(fftLower[cv+rad][cu][c].times(w[v*level])).minus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
							fft[cv+rad][cu+rad][c] 	= fftLower[cv][cu][c].minus(fftLower[cv][cu+rad][c].times(w[u*level])).minus(fftLower[cv+rad][cu][c].times(w[v*level])).plus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
						}
					}
				}
			}
		}
		return fft;
	}
	
	
	public static Complex[][] getFastFourier(double[][] img){
		int N = img.length;
		int nBits = (int) (Math.log(N)/Math.log(2));
		//Generate W array
		Complex[] w = new Complex[N];
		for(int v = 0; v<N; v++){
			w[v] = new Complex(0, -2.0 * Math.PI / (double) N * (double)v).exp();
		}
		return FFT(1,N,nBits,img,w);
	}

	private static Complex[][] FFT(int level, int oSize, int nBits, double[][] img, Complex[] w){
		int levelSize = oSize/level;
		//At lowest possible level
		if(levelSize<=2){
			int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
			double c = Math.pow(oSize, -2);
			double exp = -2.0*Math.PI/(double)levelSize;
			Complex[][] fft = new Complex[oSize][oSize];
			Complex F00, F01, F10, F11;
			for(int v = 0; v<fft.length; v+=2){
				for(int u = 0; u<fft.length; u+=2){
					x0 = bitReversal(u, nBits);
					x1 = bitReversal(u+1, nBits);
					y0 = bitReversal(v, nBits);
					y1 = bitReversal(v+1, nBits);
					if(CENTERED_FFT){
						F00 = new Complex(img[y0][x0]*Math.pow(-1, y0+x0), 0);
						F01 = new Complex(img[y0][x1]*Math.pow(-1, y0+x1), 0);
						F10 = new Complex(img[y1][x0]*Math.pow(-1, y1+x0), 0);
						F11 = new Complex(img[y1][x1]*Math.pow(-1, y1+x1), 0);
					}else{
						F00 = new Complex(img[y0][x0], 0);
						F01 = new Complex(img[y0][x1], 0);
						F10 = new Complex(img[y1][x0], 0);
						F11 = new Complex(img[y1][x1], 0);
					}
					x0=0;x1=1;y0=0;y1=1;
					// u=0;v=0
					fft[v][u] = (new Complex(c, 0)).times(F00.plus(F01).plus(F10).plus(F11));
					// u=1;v=0
					fft[v][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F01.times(new Complex(0,exp*(double)(x1)).exp())).plus //x1
							(F10.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F11.times(new Complex(0,exp*(double)(x1)).exp())) //x1
							);
					// u=0;v=1
					fft[v+1][u] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F01.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F10.times(new Complex(0,exp*(double)(y1)).exp())).plus //y1
							(F11.times(new Complex(0,exp*(double)(y1)).exp()))
							);
					// u=1;v=1
					fft[v+1][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0+y0)).exp())).plus
							(F01.times(new Complex(0,exp*(double)(x1+y0)).exp())).plus
							(F10.times(new Complex(0,exp*(double)(x0+y1)).exp())).plus
							(F11.times(new Complex(0,exp*(double)(x1+y1)).exp()))
							);
				}
			}
			return fft;
		}
		Complex[][] fftLower = FFT(level*2, oSize, nBits, img,w);//Lower Level FFT
		Complex[][] fft = new Complex[oSize][oSize];//FFT to return
		int rad = levelSize/2;
		for(int vOffset = 0; vOffset<oSize; vOffset+=levelSize){
			for(int uOffset = 0; uOffset<oSize; uOffset+=levelSize){
				for(int v = 0; v<rad; v++){
					for(int u = 0; u<rad; u++){
						int cv = vOffset+ v;
						int cu = uOffset+ u;
						fft[cv][cu] 			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv][cu+rad] 		= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu]			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu+rad] 	= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));	
					}
				}
			}
		}
		return fft;
	}
	
	public static double[][] getReverseFastFourier(Complex[][] fft){
		int N = fft.length;
		int nBits = (int) (Math.log(N)/Math.log(2));
		//Generate W array
		Complex[] w = new Complex[N];
		for(int v = 0; v<N; v++){
			w[v] = new Complex(0, 2.0 * Math.PI / (double) N * (double)v).exp();
		}
		Complex[][] rFFT = rFFT(1,N,nBits,fft,w);
		double[][] out = new double[N][N];
		for(int y=0; y<N; y++)
			for(int x=0; x<N; x++){
				out[y][x] = rFFT[y][x].real(); 
				if(CENTERED_FFT){
					out[y][x] *= Math.pow(-1, y+x);
				}
			}
		return out;
	}
	
	private static Complex[][] rFFT(int level, int oSize, int nBits, Complex[][] img, Complex[] w){
		int levelSize = oSize/level;
		//At lowest possible level
		if(levelSize<=2){
			int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
			double c = 1;
			double exp = 2.0*Math.PI/(double)levelSize;
			Complex[][] fft = new Complex[oSize][oSize];
			Complex F00, F01, F10, F11;
			for(int v = 0; v<fft.length; v+=2){
				for(int u = 0; u<fft.length; u+=2){
					x0 = bitReversal(u, 	nBits);
					x1 = bitReversal(u+1, 	nBits);
					y0 = bitReversal(v, 	nBits);
					y1 = bitReversal(v+1, 	nBits);
					F00 = img[y0][x0];
					F01 = img[y0][x1];
					F10 = img[y1][x0];
					F11 = img[y1][x1];
					x0=0;x1=1;y0=0;y1=1;
					// u=0;v=0
					fft[v][u] = (new Complex(c, 0)).times(F00.plus(F01).plus(F10).plus(F11));
					// u=1;v=0
					fft[v][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F01.times(new Complex(0,exp*(double)(x1)).exp())).plus //x1
							(F10.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F11.times(new Complex(0,exp*(double)(x1)).exp())) //x1
							);
					// u=0;v=1
					fft[v+1][u] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F01.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F10.times(new Complex(0,exp*(double)(y1)).exp())).plus //y1
							(F11.times(new Complex(0,exp*(double)(y1)).exp()))
							);
					// u=1;v=1
					fft[v+1][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0+y0)).exp())).plus
							(F01.times(new Complex(0,exp*(double)(x1+y0)).exp())).plus
							(F10.times(new Complex(0,exp*(double)(x0+y1)).exp())).plus
							(F11.times(new Complex(0,exp*(double)(x1+y1)).exp()))
							);
				}
			}
			return fft;
		}
		Complex[][] fftLower = rFFT(level*2, oSize, nBits, img,w);//Lower Level FFT
		Complex[][] fft = new Complex[oSize][oSize];//FFT to return
		int rad = levelSize/2;
		for(int vOffset = 0; vOffset<oSize; vOffset+=levelSize){
			for(int uOffset = 0; uOffset<oSize; uOffset+=levelSize){
				for(int v = 0; v<rad; v++){
					for(int u = 0; u<rad; u++){
						int cv = vOffset+ v;
						int cu = uOffset+ u;
						fft[cv][cu] 			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv][cu+rad] 		= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu]			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu+rad] 	= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));	
					}
				}
			}
		}
		return fft;
	}
	
	public static double[][][] getReverseColorFastFourier(Complex[][][] fft){
		int N = fft.length;
		int nBits = (int) (Math.log(N)/Math.log(2));
		//Generate W array
		Complex[] w = new Complex[N];
		for(int v = 0; v<N; v++){
			w[v] = new Complex(0, 2.0 * Math.PI / (double) N * (double)v).exp();
		}
		Complex[][][] rFFT = rColorFFT(1,N,nBits,fft,w);
		double[][][] out = new double[N][N][3];
		for(int y=0; y<N; y++)
			for(int x=0; x<N; x++){
				for(int c=0; c<3; c++){
					out[y][x][c] = rFFT[y][x][c].real(); 
					if(CENTERED_FFT){
						out[y][x][c] *= Math.pow(-1, y+x);
					}
				}
			}
		return out;
	}
	
	private static Complex[][][] rColorFFT(int level, int oSize, int nBits, Complex[][][] img, Complex[] w){
		int levelSize = oSize/level;
		//At lowest possible level
		if(levelSize<=2){
			int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
			double c = 1;
			double exp = 2.0*Math.PI/(double)levelSize;
			Complex[][][] fft = new Complex[oSize][oSize][3];
			Complex F00, F01, F10, F11;
			for(int v = 0; v<fft.length; v+=2){
				for(int u = 0; u<fft.length; u+=2){
					x0 = bitReversal(u, 	nBits);
					x1 = bitReversal(u+1, 	nBits);
					y0 = bitReversal(v, 	nBits);
					y1 = bitReversal(v+1, 	nBits);
					for(int co=0; co<3; co++){
						F00 = img[y0][x0][co];
						F01 = img[y0][x1][co];
						F10 = img[y1][x0][co];
						F11 = img[y1][x1][co];
						x0=0;x1=1;y0=0;y1=1;
						// u=0;v=0
						fft[v][u][co] = (new Complex(c, 0)).times(F00.plus(F01).plus(F10).plus(F11));
						// u=1;v=0
						fft[v][u+1][co] = (new Complex(c, 0)).times(
								(F00.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
								(F01.times(new Complex(0,exp*(double)(x1)).exp())).plus //x1
								(F10.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
								(F11.times(new Complex(0,exp*(double)(x1)).exp())) //x1
								);
						// u=0;v=1
						fft[v+1][u][co] = (new Complex(c, 0)).times(
								(F00.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
								(F01.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
								(F10.times(new Complex(0,exp*(double)(y1)).exp())).plus //y1
								(F11.times(new Complex(0,exp*(double)(y1)).exp()))
								);
						// u=1;v=1
						fft[v+1][u+1][co] = (new Complex(c, 0)).times(
								(F00.times(new Complex(0,exp*(double)(x0+y0)).exp())).plus
								(F01.times(new Complex(0,exp*(double)(x1+y0)).exp())).plus
								(F10.times(new Complex(0,exp*(double)(x0+y1)).exp())).plus
								(F11.times(new Complex(0,exp*(double)(x1+y1)).exp()))
								);
					}
				}
			}
			return fft;
		}
		Complex[][][] fftLower = rColorFFT(level*2, oSize, nBits, img,w);//Lower Level FFT
		Complex[][][] fft = new Complex[oSize][oSize][3];//FFT to return
		int rad = levelSize/2;
		for(int vOffset = 0; vOffset<oSize; vOffset+=levelSize){
			for(int uOffset = 0; uOffset<oSize; uOffset+=levelSize){
				for(int v = 0; v<rad; v++){
					for(int u = 0; u<rad; u++){
						int cv = vOffset+ v;
						int cu = uOffset+ u;
						for(int c=0; c<3; c++){
							fft[cv][cu][c] 			= fftLower[cv][cu][c].plus(fftLower[cv][cu+rad][c].times(w[u*level])).plus(fftLower[cv+rad][cu][c].times(w[v*level])).plus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
							fft[cv][cu+rad][c] 		= fftLower[cv][cu][c].minus(fftLower[cv][cu+rad][c].times(w[u*level])).plus(fftLower[cv+rad][cu][c].times(w[v*level])).minus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
							fft[cv+rad][cu][c]		= fftLower[cv][cu][c].plus(fftLower[cv][cu+rad][c].times(w[u*level])).minus(fftLower[cv+rad][cu][c].times(w[v*level])).minus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
							fft[cv+rad][cu+rad][c] 	= fftLower[cv][cu][c].minus(fftLower[cv][cu+rad][c].times(w[u*level])).minus(fftLower[cv+rad][cu][c].times(w[v*level])).plus(fftLower[cv+rad][cu+rad][c].times(w[(v+u)*level]));
						}
					}
				}
			}
		}
		return fft;
	}
	
	/**
	 * Fast Cousine Transform (FCT)
	 */
	public static double[][] getFastCousine(byte[][] img){
		int N = img.length;
		byte[][] nImg = new byte[2*N][2*N];
		int nV = 0, nU = 0;
		for(int v=0; v<2*N; v++){
			for(int u=0; u<2*N; u++){
				if(u>=N) nU = 2*N-u-1;
				else nU = u;
				if(v>=N) nV = 2*N-v-1;
				else nV = v;
				nImg[v][u] = img[nV][nU];
			}
		}
		CENTERED_FFT = false;
		Complex[][] fctReturn = getFastFourier(nImg);
		CENTERED_FFT = true;
		double[][] fct = new double[N][N];
		double t1 = 4.0*Math.pow(N, 2)*Math.pow(2d/(double)N, 0.5);
		double t2 = 1d/Math.pow(2.0, 0.5);
		double eFactor = -Math.PI/(double)(2*N);
		Complex e;
		for (int v=0; v<N; v++){
			for (int u=0; u<N; u++){
				e = (new Complex(0,eFactor*(double)(u+v))).exp();
				fct[v][u] = t1*(fctReturn[v][u].times(e)).real();
				if(u==0) fct[v][u] *= t2;
				if(v==0) fct[v][u] *= t2;
			}
		}
		return fct;
	}

	//Not ready yet	
	public static double[][] getReverseFastCousine(double[][] fft){
		int N = fft.length;
		double[][] nImg = new double[2*N][2*N];
		int nV = 0, nU = 0;
		for(int v=0; v<2*N; v++){
			for(int u=0; u<2*N; u++){
				if(u>=N) nU = 2*N-u-1;
				else nU = u;
				if(v>=N) nV = 2*N-v-1;
				else nV = v;
				nImg[v][u] = fft[nV][nU];
			}
		}
		CENTERED_FFT = false;
		Complex[][] fctReturn = getIFCT(nImg);
		CENTERED_FFT = true;
		double[][] fct = new double[N][N];
		double t1 = 4.0*Math.pow(N, 2)*Math.pow(2d/(double)N, 0.5);
		double t2 = 1d/Math.pow(2.0, 0.5);
		double eFactor = -Math.PI/(double)(2*N);
		Complex e;
		for (int v=0; v<N; v++){
			for (int u=0; u<N; u++){
				e = new Complex(0,eFactor*(u+v)).exp();
				fct[v][u] = t1*(fctReturn[v][u].times(e)).real();
				if(u==0) fct[v][u] *= t2;
				if(v==0) fct[v][u] *= t2;
			}
		}
		return fct;
	}
	
	private static Complex[][] getIFCT(double[][] img){
		int N = img.length;
		int nBits = (int) (Math.log(N)/Math.log(2));
		//Generate W array
		Complex[] w = new Complex[N];
		for(int v = 0; v<N; v++){
			w[v] = new Complex(0, -2.0 * Math.PI / (double) N * (double)v).exp();
		}
		return iFCT(1,N,nBits,img,w);
	}

	private static Complex[][] iFCT(int level, int oSize, int nBits, double[][] img, Complex[] w){
		int levelSize = oSize/level;
		//At lowest possible level
		if(levelSize<=2){
			int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
			double c = Math.pow(oSize, -2);
			double exp = -2.0*Math.PI/(double)levelSize;
			Complex[][] fft = new Complex[oSize][oSize];
			Complex F00, F01, F10, F11;
			for(int v = 0; v<fft.length; v+=2){
				for(int u = 0; u<fft.length; u+=2){
					x0 = bitReversal(u, nBits);
					x1 = bitReversal(u+1, nBits);
					y0 = bitReversal(v, nBits);
					y1 = bitReversal(v+1, nBits);
					if(CENTERED_FFT){
						F00 = new Complex(img[y0][x0]*Math.pow(-1, y0+x0), 0);
						F01 = new Complex(img[y0][x1]*Math.pow(-1, y0+x1), 0);
						F10 = new Complex(img[y1][x0]*Math.pow(-1, y1+x0), 0);
						F11 = new Complex(img[y1][x1]*Math.pow(-1, y1+x1), 0);
					}else{
						F00 = new Complex(img[y0][x0], 0);
						F01 = new Complex(img[y0][x1], 0);
						F10 = new Complex(img[y1][x0], 0);
						F11 = new Complex(img[y1][x1], 0);
					}
					x0=0;x1=1;y0=0;y1=1;
					// u=0;v=0
					fft[v][u] = (new Complex(c, 0)).times(F00.plus(F01).plus(F10).plus(F11));
					// u=1;v=0
					fft[v][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F01.times(new Complex(0,exp*(double)(x1)).exp())).plus //x1
							(F10.times(new Complex(0,exp*(double)(x0)).exp())).plus //x0
							(F11.times(new Complex(0,exp*(double)(x1)).exp())) //x1
							);
					// u=0;v=1
					fft[v+1][u] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F01.times(new Complex(0,exp*(double)(y0)).exp())).plus //y0
							(F10.times(new Complex(0,exp*(double)(y1)).exp())).plus //y1
							(F11.times(new Complex(0,exp*(double)(y1)).exp()))
							);
					// u=1;v=1
					fft[v+1][u+1] = (new Complex(c, 0)).times(
							(F00.times(new Complex(0,exp*(double)(x0+y0)).exp())).plus
							(F01.times(new Complex(0,exp*(double)(x1+y0)).exp())).plus
							(F10.times(new Complex(0,exp*(double)(x0+y1)).exp())).plus
							(F11.times(new Complex(0,exp*(double)(x1+y1)).exp()))
							);
				}
			}
			return fft;
		}
		Complex[][] fftLower = FFT(level*2, oSize, nBits, img,w);//Lower Level FFT
		Complex[][] fft = new Complex[oSize][oSize];//FFT to return
		int rad = levelSize/2;
		for(int vOffset = 0; vOffset<oSize; vOffset+=levelSize){
			for(int uOffset = 0; uOffset<oSize; uOffset+=levelSize){
				for(int v = 0; v<rad; v++){
					for(int u = 0; u<rad; u++){
						int cv = vOffset+ v;
						int cu = uOffset+ u;
						fft[cv][cu] 			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv][cu+rad] 		= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).plus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu]			= fftLower[cv][cu].plus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).minus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));
						fft[cv+rad][cu+rad] 	= fftLower[cv][cu].minus(fftLower[cv][cu+rad].times(w[u*level])).minus(fftLower[cv+rad][cu].times(w[v*level])).plus(fftLower[cv+rad][cu+rad].times(w[(v+u)*level]));	
					}
				}
			}
		}
		return fft;
	}
	/**
	 * END Fast Cousine Transform	
	 */
	
	public static int bitReversal(int number, int nBits){
		int out = 0;
		for(int n = 0; n<nBits; n++){
			out += (number & (1 << n))>>n<<(nBits-1-n);
		}
		return out;
	}
	
	public static Complex[][] getLoG(Complex[][] fft, double rho){
		int N = fft.length;
		int c = N/2;
		Complex[][] out = new Complex[N][N];
		double sos;
		double factor;
		for(int v=0; v<N; v++){
			for(int u=0; u<N; u++){
				sos = Math.pow(u-c,2)+Math.pow(c-v,2);
				factor = -(sos)*Math.exp(-0.5*sos*Math.pow(rho,2));
				out[v][u] = fft[v][u].times(new Complex(factor,0));
			}
		}
		return out;
	}
	
	public static double[][] zeroPadding(double[][] img){
		int h = img.length;
		int w = img[0].length;
		int N0 = w;
		if(h>N0) N0=h;
		int nBits = (int) (Math.log(N0)/Math.log(2));
		int N = (int) Math.pow(2, nBits);
		if(N<N0)
			N = (int) Math.pow(2, nBits+1);
		double[][] out = new double[N][N];
		int y0 = (N-h)/2;
		int x0 = (N-w)/2;
		for(int y=0; y<h; y++)
			for(int x=0; x<w; x++)
				out[y0+y][x0+x] = img[y][x];		
		return out;
	}
	
	public static byte[][] zeroPadding(byte[][] img){
		int h = img.length;
		int w = img[0].length;
		int N0 = w;
		if(h>N0) N0=h;
		int nBits = (int) (Math.log(N0)/Math.log(2));
		int N = (int) Math.pow(2, nBits);
		if(N<N0)
			N = (int) Math.pow(2, nBits+1);
		byte[][] out = new byte[N][N];
		int y0 = (N-h)/2;
		int x0 = (N-w)/2;
		for(int y=0; y<h; y++)
			for(int x=0; x<w; x++)
				out[y0+y][x0+x] = img[y][x];		
		return out;
	}
	
	public static byte[][][] zeroPadding(byte[][][] img){
		int h = img.length;
		int w = img[0].length;
		int N0 = w;
		if(h>N0) N0=h;
		int nBits = (int) (Math.log(N0)/Math.log(2));
		int N = (int) Math.pow(2, nBits);
		if(N<N0)
			N = (int) Math.pow(2, nBits+1);
		byte[][][] out = new byte[N][N][img[0][0].length];
		int y0 = (N-h)/2;
		int x0 = (N-w)/2;
		for(int y=0; y<h; y++)
			for(int x=0; x<w; x++)
				for(int c = 0; c<img[0][0].length; c++)
					out[y0+y][x0+x][c] = img[y][x][c];		
		return out;
	}
	
	public static float[] getHistogram(byte[][] img){
		float[] hist = new float[256];
		int h = img.length, w = img[0].length;
		float weight = (float) (1.0/((float)h * (float)w));
		for(int y=0; y<h; y++)
			for (int x=0; x<w; x++)
				hist[img[y][x]&0x0ff] += weight;
		return hist;
	}
	
	public static float[][] getColorHistogram(byte[][][] img){
		float[][] hist = new float[img[0][0].length][256];
		int h = img.length, w = img[0].length, c = img[0][0].length;
		float weight = (float) (1.0/((float)h * (float)w));
		for(int y=0; y<h; y++)
			for (int x=0; x<w; x++)
				for(int chan = 0; chan< c; chan++)
					hist[chan][img[y][x][chan]&0x0ff] += weight;
		return hist;
	}
	
	public static byte[][] equalizeHistogram(byte[][] img){
		int h = img.length, w = img[0].length;
		byte[][] nImg = new byte[h][w];
		float[] hist = getHistogram(img);
		float[] prob = new float[256];
		prob[0] = hist[0];
		for(int i = 1; i<256; i++)
			prob[i] = prob[i-1]+hist[i];		
		for(int y=0; y<h; y++)
			for (int x=0; x<w; x++)
				nImg[y][x] = (byte) (prob[img[y][x] & 0x0FF] * 255);
		return nImg;
	}
	
	public static byte[][][] equalizeColorHistogram(byte[][][] img){
		int h = img.length, w = img[0].length, cols = img[0][0].length;
		byte[][][] nImg = new byte[h][w][cols];
		float[][] hist = getColorHistogram(img);
		float[][] prob = new float[cols][256];
		prob[0] = hist[0];
		for(int i = 1; i<256; i++)
			for(int c = 0; c<cols; c++)
				prob[c][i] = prob[c][i-1]+hist[c][i];		
		for(int y=0; y<h; y++)
			for (int x=0; x<w; x++)
				for(int c = 0; c<cols; c++)
					nImg[y][x][c] = (byte) (prob[c][img[y][x][c] & 0x0FF] * 255);
		return nImg;
	}
	
	public static byte[][] threshold(byte[][] img, int th){
		int h0 = img.length;
		int w0 = img[0].length;
		byte[][] out = new byte[h0][w0];
		for(int y=0; y<h0; y++)
			for(int x=0; x<w0; x++)
				if( (img[y][x] & 0x0FF) > th)
					out[y][x] = (byte) 255;
				else
					out[y][x] = 0;
		return out;
	}
	
	public static byte[][] threshold(byte[][] img, int thP, int thN){
		int N = img.length;
		byte[][] out = new byte[N][N];
		for(int y=0; y<N; y++)
			for(int x=0; x<N; x++)
				if( (img[y][x] & 0x0FF)-127 > thP || (img[y][x] & 0x0FF)-127 < thN)
					out[y][x] = (byte) 255;
				else
					out[y][x] = 0;
		return out;
	}
	
	public static double[][] getKirsh(double[][] img){
		int N = img.length;
		if((N==0) || (img[0].length != N)) return null;
		double[][] out = new double[N][N];
		int[] aX = new int[]{-1, 0, 1, 1, 1, 0,-1,-1};
		int[] aY = new int[]{-1,-1,-1, 0, 1, 1, 1, 0};
		int[] ind = new int[8];
		double s = 0, t = 0, g = 0, max = 0;
		for(int y = 0; y<N-3; y++){
			for(int x = 0; x<N-3; x++){
				g=0; max=0;
				for(int i = 0; i<8; i++){
					for(int in = 0; in<8; in++){
						ind[in] = i+in;
						if(ind[in]>7) ind[in] = ind[in] - 8; 
					}
					s = img[1+y+aY[ind[0]]][1+x+aX[ind[0]]] + img[1+y+aY[ind[1]]][1+x+aX[ind[1]]] + img[1+y+aY[ind[2]]][1+x+aX[ind[2]]];
					t = img[1+y+aY[ind[3]]][1+x+aX[ind[3]]] + img[1+y+aY[ind[4]]][1+x+aX[ind[4]]] + img[1+y+aY[ind[5]]][1+x+aX[ind[5]]] 
							+ img[1+y+aY[ind[6]]][1+x+aX[ind[6]]] + img[1+y+aY[ind[7]]][1+x+aX[ind[7]]];
					g = Math.abs(5*s - 3*t);
					if(g > max) max = g;
				}
				out[y+1][x+1] = max;
			}
		}		
		return out;
	}

	public static double[][] getSobelMag(double[][] img){
		double[][][] out = getSobel(img);
		if(out==null)return null;
		return out[0];
	}
	
	public static double[][][] getSobel(double[][] img){
		int N = img.length;
		if((N==0) || (img[0].length != N)) return null;
		double[][][] out = new double[2][N][N];
		float[][] s1 	= new float[][]{
				new float[]{-1, 0, 1},
				new float[]{-2, 0, 2},
				new float[]{-1, 0, 1}
				};
		float[][] s2 	= new float[][]{
				new float[]{ 1, 2, 1},
				new float[]{ 0, 0, 0},
				new float[]{-1,-2,-1}
				};
		double d1 = 0, d2 = 0;
		for(int y = 0; y<N-3; y++){
			for(int x = 0; x<N-3; x++){
				d1=0;d2=0;
				for(int v = 0; v<3; v++){
					for(int u = 0; u<3; u++){
						d1 += img[y+v][x+v]*s1[v][u];
						d2 += img[y+v][x+v]*s2[v][u];
					}
				}
				out[0][y+1][x+1] = Math.pow(Math.pow(d1, 2)+Math.pow(d2, 2), 0.5);
				out[1][y+1][x+1] = Math.atan(Math.abs(d2)/Math.abs(d1));
			}
		}		
		return out;
	}
	
	public static double[][] nonMaximumSupression(double[][][] sobelImage){
		if(sobelImage == null) return null;
		int[] s = new int[]{sobelImage[0].length, sobelImage[0][0].length};
		double[][] out = new double[s[0]][s[1]];
		double c = 180.0/(45.0*Math.PI);
		for(int y=1; y<s[0]-1; y++){
			for(int x=1; x<s[1]-1; x++){
				int a = (int)Math.round(sobelImage[1][y][x]*c);
				if(a>180){
					int t = (int)(a/180);
					a -= t*180;
				}
				double t = sobelImage[0][y][x];
				switch(a){
				case 0:	//0 deg
					if(t>sobelImage[0][y][x-1] && t>sobelImage[0][y][x+1])
						out[y][x] = t;
					else
						out[y][x] = 0;
					break;
				case 1:	//45 deg
					if(t>sobelImage[0][y-1][x+1] && t>sobelImage[0][y+1][x-1])
						out[y][x] = t;
					else
						out[y][x] = 0;
					break;
				case 2:	//90 deg
					if(t>sobelImage[0][y-1][x] && t>sobelImage[0][y+1][x])
						out[y][x] = t;
					else
						out[y][x] = 0;
					break;
				case 3:	//135 deg
					if(t>sobelImage[0][y-1][x-1] && t>sobelImage[0][y+1][x+1])
						out[y][x] = t;
					else
						out[y][x] = 0;
					break;
				case 4:	//180 deg
					if(t>sobelImage[0][y][x-1] && t>sobelImage[0][y][x+1])
						out[y][x] = t;
					else
						out[y][x] = 0;
					break;
				default:
					System.out.println("Unknown angle: "+a);
				}
			}
		}
		return out;
	}
	
	public static byte[][] rotate(byte[][] img, double theta){
		int h = img.length, w = img[0].length;
		double hm = h/2, wm = w/2;
		byte[][] out = new byte[h][w];
		double sin = Math.sin(-theta*Math.PI/180.0);
		double cos = Math.cos(-theta*Math.PI/180.0);
		int xN=0, yN=0; 
		for(int y=0; y<h; y++){
			for(int x=0; x<w; x++){
				xN = (int) (wm + ((x-wm)*cos-(hm-y)*sin));
				yN = (int) (hm - ((x-wm)*sin+(hm-y)*cos));
				if(xN>0 && xN<w && yN>0 && yN<h)
					out[y][x] = img[yN][xN];
			}
		}
		return out;
	}
	
	public static double[][] rotate(double[][] img, double theta){
		int h = img.length, w = img[0].length;
		double hm = h/2, wm = w/2;
		double[][] out = new double[h][w];
		double sin = Math.sin(-theta*Math.PI/180.0);
		double cos = Math.cos(-theta*Math.PI/180.0);
		int xN=0, yN=0; 
		for(int y=0; y<h; y++){
			for(int x=0; x<w; x++){
				yN = (int) (hm - ((x-wm)*sin+(hm-y)*cos));
				xN = (int) (wm + ((x-wm)*cos-(hm-y)*sin));
				if(xN>0 && xN<w && yN>0 && yN<h)
					out[y][x] = img[yN][xN];
			}
		}
		return out;
	}
	
	public static byte[][] crop(int x0, int y0, int x1, int y1, byte[][] img){
		byte[][] out = new byte[y1-y0][x1-x0];
		for(int y=y0; y<y1; y++)
			for(int x=x0; x<x1; x++)
				out[y-y0][x-x0] = img[y][x];
		return out;
	}
	

	public static double[][] crop(int x0, int y0, int x1, int y1, double[][] img){
		double[][] out = new double[y1-y0][x1-x0];
		for(int y=y0; y<y1; y++)
			for(int x=x0; x<x1; x++)
				out[y-y0][x-x0] = img[y][x];
		return out;
	}
	
	public static double[][] crop(int h, int w, double[][] img){
		int y0 = img.length/2-h/2, y1 = img.length/2+h/2;
		int x0 = img[0].length/2-w/2, x1 = img[0].length/2+w/2; 
		return crop(x0, y0, x1, y1, img);
	}
	
	public static byte[][] getNegative(byte[][] img){
		byte[][] out = new byte[img.length][img[0].length];
		for(int y=0; y<img.length; y++){
			for(int x=0; x<img[0].length; x++){
				out[y][x] = (byte) (255-img[y][x]);
			}
		}
		return out;
	}
	
	public static byte[][] hMotion(byte[][] img, double a){
		int[] size = new int[]{img.length, img[0].length};
		Complex[][] fft = getFastFourier(zeroPadding(img));
		byte[][] out = normilize(getReverseFastFourier(hMotion(fft, a)));
		int yOff=(out.length-size[0])/2, xOff=(out[0].length-size[1])/2;
		return crop(xOff, yOff, out[0].length-xOff, out.length-yOff, out);
	}
	
	public static Complex[][] hMotion(Complex[][] fft, double a){
		int[] size = new int[]{fft.length, fft[0].length};
		Complex[][] out = new Complex[size[0]][size[1]];	
		for(int v=0; v<size[0]; v++){
			for(int u=0; u<size[1]; u++){
				double hVal = sinc(a*(u-size[1]/2)*Math.PI);
				out[v][u] = fft[v][u].times(new Complex(hVal,0)).times(new Complex(0.0,-hVal).exp());
			}
		}
		return out;
	}
	
	public static double[][] atmosphericTurbulence(double[][] img, double a){
		int[] size = new int[]{img.length, img[0].length};
		double[][] out = new double[size[0]][size[1]];	
		double c1 = -Math.PI*Math.pow(a,2);
		for(int v=0; v<size[0]; v++){
			double vC = (size[0]/2-v);
			for(int u=0; u<size[1]; u++){
				double uC = (u-size[1]/2);
				out[v][u] = img[v][u]*Math.exp(c1*(Math.pow(vC,2)+Math.pow(uC,2)));
			}
		}
		return out;
	}
	
	public static byte[][] atmosphericTurbulence(byte[][] img, double a){
		int[] size = new int[]{img.length, img[0].length};
		Complex[][] fft = getFastFourier(zeroPadding(img));
		byte[][] out = normilize(getReverseFastFourier(atmosphericTurbulence(fft, a)));
		int yOff=(out.length-size[0])/2, xOff=(out[0].length-size[1])/2;
		return crop(xOff, yOff, out[0].length-xOff, out.length-yOff, out);
	}
	
	public static Complex[][] atmosphericTurbulence(Complex[][] fft, double a){
		int[] size = new int[]{fft.length, fft[0].length};
		Complex[][] out = new Complex[size[0]][size[1]];	
		double c1 = 1.0/Math.pow(a,2);
		double c2 = -Math.PI*c1;
		for(int v=0; v<size[0]; v++){
			double vC = (size[0]/2-v);
			for(int u=0; u<size[1]; u++){
				double uC = (u-size[1]/2);
				out[v][u] = fft[v][u].times(new Complex(c1*Math.exp(c2*(Math.pow(vC,2)+Math.pow(uC,2))),0));
			}
		}
		return out;
	}
		
	public static double sinc(double x){
		if(x==0)
			return 1;
		return Math.sin(x)/x;
	}
	
	public static double[][][] getLaplacian(double[][] dImg){
		if(dImg==null) return null;
		int[] s = new int[]{dImg.length, dImg[0].length};
		double[][][] out = new double[2][s[0]][s[1]];
//		double dy = 1.0/(double)s[0], dx = 1.0/(double)s[1];
		double ly, lx;
		for(int y=0; y<s[0]-1; y++){
			for(int x=0; x<s[1]-1; x++){
				ly = (dImg[y][x]-dImg[y+1][x]);
				lx = (dImg[y][x]-dImg[y][x+1]);
				out[0][y][x] = Math.sqrt(Math.pow(lx, 2)+Math.pow(ly, 2));
				out[1][y][x] = Math.atan(Math.abs(ly)/Math.abs(lx));
			}
		}
		return out;
	}
	
	public static float[][][] getLaplacian(float[][] fImg){
		if(fImg==null) return null;
		int[] s = new int[]{fImg.length, fImg[0].length};
		float[][][] out = new float[2][s[0]][s[1]];
//		double dy = 1.0/(double)s[0], dx = 1.0/(double)s[1];
		double ly, lx;
		for(int y=0; y<s[0]-1; y++){
			for(int x=0; x<s[1]-1; x++){
				ly = (fImg[y][x]-fImg[y+1][x]);
				lx = (fImg[y][x]-fImg[y][x+1]);
				out[0][y][x] = (float) Math.sqrt(Math.pow(lx, 2)+Math.pow(ly, 2));
				out[1][y][x] = (float) Math.atan(Math.abs(ly)/Math.abs(lx));
			}
		}
		return out;
	}

	public static double[][] getLaplacianMag(double[][] img){
		return getLaplacian(img)[0];
	}
	
	public static float[][] getLaplacianMag(float[][] img){
		return getLaplacian(img)[0];
	}
	
	public static double[][] getLaplacianMag2(double[][] img){
		return getLaplacian2(img)[0];
	}
	
	public static double[][][] getLaplacian2(double[][] img){
		if(img==null) return null;
		int[] s = new int[]{img.length, img[0].length};
		double[][][] out = new double[2][s[0]][s[1]];
//		double dy = 1.0/(double)s[0], dx = 1.0/(double)s[1];
		double ly, lx;
		for(int y=1; y<s[0]-1; y++){
			for(int x=1; x<s[1]-1; x++){
				ly = (img[y+1][x] - img[y-1][x]);
				lx = (img[y][x+1] - img[y][x-1]);
				out[0][y][x] = Math.sqrt(Math.pow(lx, 2)+Math.pow(ly, 2));
				out[1][y][x] = Math.atan(Math.abs(ly)/Math.abs(lx));
			}
		}
		return out;
	}
	
	public static double[][] byte2double(byte[][] a){
		if(a==null) return null;
		int[] s = new int[]{a.length, a[0].length};
		double[][] out = new double[s[0]][s[1]];
		for(int y=0; y<s[0]; y++){
			for(int x=0; x<s[1]; x++){
				out[y][x] = (a[y][x]&0x0FF);
			}
		}
		return out;
	}
	
	public static double[][][] byte2double(byte[][][] a){
		if(a==null) return null;
		int[] s = new int[]{a.length, a[0].length, a[0][0].length};
		double[][][] out = new double[s[0]][s[1]][s[2]];
		for(int c=0; c<s[0]; c++){
			for(int y=0; y<s[1]; y++){
				for(int x=0; x<s[2]; x++){
					out[c][y][x] = (a[c][y][x]&0x0FF);
				}
			}
		}
		return out;
	}
	
	public static float[][] byte2Float(byte[][] a){
		if(a==null) return null;
		int[] s = new int[]{a.length, a[0].length};
		float[][] out = new float[s[0]][s[1]];
		for(int y=0; y<s[0]; y++){
			for(int x=0; x<s[1]; x++){
				out[y][x] = (a[y][x]&0x0FF);
			}
		}
		return out;
	}
	
	public static float[][][] byte2Float(byte[][][] a){
		if(a==null) return null;
		int[] s = new int[]{a.length, a[0].length, a[0][0].length};
		float[][][] out = new float[s[0]][s[1]][s[2]];
		for(int c=0; c<s[0]; c++){
			for(int y=0; y<s[1]; y++){
				for(int x=0; x<s[2]; x++){
					out[c][y][x] = (a[c][y][x]&0x0FF);
				}
			}
		}
		return out;
	}

	public static double[][] convolve(double[][] img, double[][] k){
		int[] s1 = new int[]{img.length, img[0].length};
		int[] s2 = new int[]{k.length, k[0].length};
		double[][] out = new double[s1[0]][s1[1]];
		double t;
		int xOff = s2[0]/2, yOff = s2[1]/2;
		for(int y=0; y<s1[0]-s2[0]; y++){
			for(int x=0; x<s1[1]-s2[1]; x++){
				t = 0;
				for(int yk=0; yk<s2[0]; yk++){
					for(int xk=0; xk<s2[1]; xk++){
						t += img[y+yk][x+xk]*k[yk][xk];
					}
				}
				out[y+yOff][x+xOff] = t;
			}
		}
		return out;
	}
	
	public static double[][] convolveCell(double[][] img, double[][] k){
		int[] s1 = new int[]{img.length, img[0].length};
		int[] s2 = new int[]{k.length, k[0].length};
		double[][] out = new double[s1[0]][s1[1]];
		double t;
		int xOff = s2[0]/2, yOff = s2[1]/2;
		for(int y=0; y<s1[0]-yOff; y++){
			for(int x=0; x<s1[1]-xOff; x++){
				t = 0;
				for(int yk=0; yk<s2[0]; yk++){
					for(int xk=0; xk<s2[1]; xk++){
						t += img[y+yk][x+xk]*k[yk][xk];
					}
				}
				if(t<0) t=0;
				out[y+yOff][x+xOff] = t;
			}
		}
		return out;
	}
	
	public static double[][] abs(double[][] img){
		int[] s1 = new int[]{img.length, img[0].length};
		double[][] out = new double[s1[0]][s1[1]];
		for(int y=0; y<s1[0]; y++){
			for(int x=0; x<s1[1]; x++){
				out[y][x] = Math.abs(img[y][x]);
			}
		}
		return out;
	}

	public static double[][] nonMaximumSupression2(double[][][] img){
		if(img == null) return null;
		int[] s = new int[]{img[0].length, img[0][0].length};
		double[][] out = new double[s[0]][s[1]];
		double[][][] k = new double[19][][];
		for(int i=0; i<19; i++){
			k[i]=getKernel(10*i);
		}
		double c = 180.0/(Math.PI);
		for(int y=0; y<s[0]-2; y++){
			for(int x=0; x<s[1]-2; x++){
				double an = img[1][y+1][x+1]*c;
				if(an>180){
					int t = (int)(an/180);
					an -= t*180;
				}
				int a = (int)Math.round(an/10);
				double[][] uK = k[a];
				double m = img[0][y][x]*uK[0][0];
				for(int ky=0; ky<3; ky++){
					for(int kx=0; kx<3; kx++){
						double t = img[0][y+ky][x+kx]*uK[ky][kx];
						if(t>m) m=t;
					}
				}
				if(img[0][y+1][x+1]*uK[1][1] >= m){
					out[y+1][x+1] = img[0][y+1][x+1];
				}else{
					out[y+1][x+1] = 0;
				}
			}
		}
		return out;
	}
	
	public static double[][] getKernel(double a){
		if(a>=180){// Only use 180 degrees
			int t = (int) (a/180);
			a-=t*180;
		}
		double[][] out = new double[3][3];
		double w = Math.sqrt(0.5);
		double a1 = Math.atan(1.0/3.0)*180.0/Math.PI;
//		System.out.println("a1: "+a1);
		double x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5;
		if (a>=0 && a<a1){
			x4 = w/2;
			y4 = x4*Math.tan(a*Math.PI/180.0);
			x5 = 1.5*w;
			y5 = x5*Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[1][2] = Math.sqrt((x5-x4)*(x5-x4)+(y5-y4)*(y5-y4));
			out[1][0] = out[1][2];
		}
		if (a>=a1 && a<45){
			x4 = w/2;
			y4 = x4*Math.tan(a*Math.PI/180.0);
			y5 = w/2;
			x5 = y5/Math.tan(a*Math.PI/180.0);
			x2 = 1.5*w;
			y2 = x2*Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[1][2] = Math.sqrt((x5-x4)*(x5-x4)+(y5-y4)*(y5-y4));
			out[0][2] = Math.sqrt((x2-x5)*(x2-x5)+(y2-y5)*(y2-y5));
			out[1][0] = out[1][2];
			out[2][0] = out[0][2];
		}
		if (a>=45 && a<(90-a1)){
			y4 = w/2;
			x4 = y4/Math.tan(a*Math.PI/180.0);
			x1 = w/2;
			y1 = x1*Math.tan(a*Math.PI/180.0);
			y2 = 1.5*w;
			x2 = y2/Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[0][1] = Math.sqrt((x1-x4)*(x1-x4)+(x1-x4)*(x1-x4));
			out[0][2] = Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
			out[2][1] = out[0][1];
			out[2][0] = out[0][2];
		}
		if (a>=(90-a1) && a<90){
			y4 = w/2;
			x4 = y4/Math.tan(a*Math.PI/180.0);
			y1 = 1.5*w;
			x1 = y1/Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[0][1] = Math.sqrt((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4));
			out[2][1] = out[0][1];
		}
		if (a>=90 && a<(135-a1)){
			y4 = w/2;
			x4 = y4/Math.tan(a*Math.PI/180.0);
			y1 = 1.5*w;
			x1 = y1/Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[0][1] = Math.sqrt((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4));
			out[2][1] = out[0][1];
		}
		if (a>=(135-a1) && a<135){
			y4 = w/2;
			x4 = y4/Math.tan(a*Math.PI/180.0);
			x1 = w/2;
			y1 = x1*Math.tan(a*Math.PI/180.0);
			y0 = 1.5*w;
			x0 = y0/Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[0][1] = Math.sqrt((x1-x4)*(x1-x4)+(x1-x4)*(x1-x4));
			out[0][0] = Math.sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
			out[2][1] = out[0][1];
			out[2][2] = out[0][0];
		}
		if (a>=135 && a<180-a1){
			x4 = w/2;
			y4 = x4*Math.tan(a*Math.PI/180.0);
			y3 = w/2;
			x3 = y3/Math.tan(a*Math.PI/180.0);
			x0 = 1.5*w;
			y0 = x0*Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[1][2] = Math.sqrt((x3-x4)*(x3-x4)+(y3-y4)*(y3-y4));
			out[0][2] = Math.sqrt((x0-x3)*(x0-x3)+(y0-y3)*(y0-y3));
			out[1][0] = out[1][2];
			out[2][0] = out[0][2];
		}
		if (a>=(180-a1) && a<180){
			x4 = w/2;
			y4 = x4*Math.tan(a*Math.PI/180.0);
			x3 = 1.5*w;
			y3 = x3*Math.tan(a*Math.PI/180.0);
			out[1][1] = 2*Math.sqrt(x4*x4+y4*y4);
			out[1][0] = Math.sqrt((x3-x4)*(x3-x4)+(y3-y4)*(y3-y4));
			out[1][2] = out[1][0];
		}
		return out;
	}
	
	public static byte[][] getHisteresisEdges(double[][] img, double minTh, double maxTh){
		int h0 = img.length;
		int w0 = img[0].length;
		byte[][] out = new byte[h0][w0];
		byte[][] temp = new byte[h0][w0];
		double t = 0;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				t = img[y][x];
				if(t>=maxTh){
					temp[y][x] = 3; 
				}
				if(t>=minTh && t<maxTh){
					temp[y][x] = 2; 
				}
				if(t<minTh){
					temp[y][x] = 1; 
				}
			}
		}
		boolean hit = false;
		for (int y = 1; y < h0-1; y++) {
			for (int x = 1; x < w0-1; x++) {
				if(temp[y][x] == 3){
					out[y][x] = (byte)255; 
				}
				if(temp[y][x] == 2){
					hit = false;
					for(int y1 = -1; y1<2; y1++){
						for(int x1 = -1; x1<2; x1++){
							if(temp[y-y1][x-x1] == 3){
								hit = true;
							}
						}
					}
					if(hit){
						out[y][x] = (byte)255;
					}else{
						out[y][x] = (byte)0;
					}
				}
				if(temp[y][x] == 1){
					out[y][x] = (byte)0; 
				}
			}
		}
		return out;
	}
	
	public static byte[][] getHisteresisEdges(double[][] img, double minTh, double maxTh, int nPasses){
		int h0 = img.length;
		int w0 = img[0].length;
		byte[][] out = new byte[h0][w0];
		byte[][] temp = new byte[h0][w0];
		double t = 0;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				t = img[y][x];
				if(t>=maxTh){
					temp[y][x] = 3; 
				}
				if(t>=minTh && t<maxTh){
					temp[y][x] = 2; 
				}
				if(t<minTh){
					temp[y][x] = 1; 
				}
			}
		}
		boolean hit = false;
		for (int y = 1; y < h0-1; y++) {
			for (int x = 1; x < w0-1; x++) {
				if(temp[y][x] == 3){
					out[y][x] = (byte)255; 
				}
				if(temp[y][x] == 2){
					hit = false;
					for(int y1 = -1; y1<2; y1++){
						for(int x1 = -1; x1<2; x1++){
							if(temp[y-y1][x-x1] == 3){
								hit = true;
							}
						}
					}
					if(hit){
						out[y][x] = (byte)255;
						temp[y][x] = 3;
					}else{
						out[y][x] = (byte)0;
					}
				}
				if(temp[y][x] == 1){
					out[y][x] = (byte)0; 
				}
			}
		}
		if(nPasses>=0){
			for(int n = 1; n<nPasses; n++){
				for (int y = 1; y < h0-1; y++) {
					for (int x = 1; x < w0-1; x++) {
						if(temp[y][x] == 2){
							hit = false;
							for(int y1 = -1; y1<2; y1++){
								for(int x1 = -1; x1<2; x1++){
									if(temp[y-y1][x-x1] == 3){
										hit = true;
									}
								}
							}
							if(hit){
								out[y][x] = (byte)255;
								temp[y][x] = 3;
							}else{
								out[y][x] = (byte)0;
							}
						}
					}
				}
			}
		}else{
			int edgesChanges = 0;
			do{
				edgesChanges = 0;
				for (int y = 1; y < h0-1; y++) {
					for (int x = 1; x < w0-1; x++) {
						if(temp[y][x] == 2){
							hit = false;
							for(int y1 = -1; y1<2; y1++){
								for(int x1 = -1; x1<2; x1++){
									if(temp[y-y1][x-x1] == 3){
										hit = true;
									}
								}
							}
							if(hit){
								out[y][x] = (byte)255;
								temp[y][x] = 3;
								edgesChanges++;
							}else{
								out[y][x] = (byte)0;
							}
						}
					}
				}
			}while(edgesChanges>0);
		}
		return out;
	}

	public static int getOtsuTh(byte[][] img){
		int maxStd = 0;
		int out = 0;
		double std=0, w0=0, w1=0, u0=0, u1=0;
		float[] hist = getHistogram(img);
		for(int i=0; i<256; i++){
			w1 	+= hist[i];
			u1 	+= i*hist[i];
		}
		for(int t=1; t<256; t++){
			w0 += hist[t-1];
			w1 -= hist[t-1];
			u0 += (t-1)*hist[t-1];
			u1 -= (t-1)*hist[t-1];
			std = w0*w1*Math.pow(u0/w0-u1/w1, 2);
			if(std>maxStd){
				maxStd = (int)Math.round(std);
				out = t;
			}
		}
		return out;
	}
	
	public static double[][] normalizeImage(double[][] img){
		int[] l = new int[]{img.length, img[0].length};
		double[][] out = new double[l[0]][l[1]];
		double c = 1.0/255.0;
		for(int y=0; y<l[0]; y++){
			for(int x=0; x<l[1]; x++){
				out[y][x] = img[y][x]*c;
			}
		}
		return out;
	}
	
	/**
	 * Normalize the float image in any range to the 0 to 255 range.  
	 * @param img
	 * @return
	 */
	public static byte[][] normalize(float[][] img){
		int h0 = img.length;
		int w0 = img[0].length;
		
		byte[][] out = new byte[h0][w0];
		double max = 0, min = 0, range = 0;
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				if(img[y][x] > max) max = img[y][x];
				if(img[y][x] < min) min = img[y][x];  
			}
		}
		range = max-min;
//		System.out.println("Max: "+max+", Min: "+min);
		for (int y = 0; y < h0; y++) {
			for (int x = 0; x < w0; x++) {
				out[y][x] = (byte) ((img[y][x]-min)/range*255); 
			}
		}
		return out;
	}

	public static byte[][] toGrayscale(byte[][][] img){
		int[] l = new int[]{img[0].length, img[0][0].length};
		byte[][] out = new byte[l[0]][l[1]];
		for(int y=0; y<l[0]; y++){
			for(int x=0; x<l[1]; x++){
				out[y][x] = (byte) ((img[0][y][x]/3+img[1][y][x]/3+img[2][y][x]/3));
			}
		}
		return out;
	}
	
	public static double[][] getEdgeKernel(double a){
		double[][] out = new double[2][2];
		if(a>=360){
			int t = (int)(a/360);
			a -= t*360;
		}
		double w = Math.sqrt(0.5);
		double a0 = w*w;
		double dx, dy;
		if(a>=0 && a<45){
			dx = w;
			dy = dx*Math.tan(a*Math.PI/180.0);
			out[0][0] = +a0;
			out[0][1] = (a0-Math.abs(dx*dy));
			out[1][0] = (Math.abs(dx*dy)-a0);
			out[1][1] = -a0;
			return out;
		}
		if(a>=45 && a<135){
			dy = w;
			dx = dy/Math.tan(a*Math.PI/180.0);
			if(a<=90){
				out[0][0] = +a0;
				out[0][1] = (Math.abs(dx*dy)-a0);
				out[1][0] = (a0-Math.abs(dx*dy));
				out[1][1] = -a0;
			}else{
				out[0][0] = (a0-Math.abs(dx*dy));
				out[0][1] = -a0;
				out[1][0] = +a0;
				out[1][1] = (Math.abs(dx*dy)-a0);
			}
			return out;
		}
		if(a>=135 && a<225){
			dx = w;
			dy = dx*Math.tan(a*Math.PI/180.0);
			if(a<=180){
				out[0][0] = Math.abs(dx*dy)-a0;
				out[0][1] = -a0;
				out[1][0] = +a0;
				out[1][1] = a0-Math.abs(dx*dy);
			}else{
				out[0][0] = -a0;
				out[0][1] = Math.abs(dx*dy)-a0;
				out[1][0] = a0-Math.abs(dx*dy);
				out[1][1] = +a0;
			}
			return out;
		}
		if(a>=225 && a<315){
			dy = w;
			dx = dy/Math.tan(a*Math.PI/180.0);
			if(a<=270){
				out[0][0] = -a0;
				out[0][1] = a0-Math.abs(dx*dy);
				out[1][0] = Math.abs(dx*dy)-a0;
				out[1][1] = +a0;
			}else{
				out[0][0] = Math.abs(dx*dy)-a0;
				out[0][1] = +a0;
				out[1][0] = -a0;
				out[1][1] = a0-Math.abs(dx*dy);
			}
			return out;
		}
		if(a>=315 && a<360){
			dx = w;
			dy = dx*Math.tan(a*Math.PI/180.0);
			out[0][0] = a0-Math.abs(dx*dy);
			out[0][1] = +a0;
			out[1][0] = -a0;
			out[1][1] = Math.abs(dx*dy)-a0;
			return out;
		}
		return out;
	}
	
	public static double[][][] getNormColors(double[][][] img){
		if(img==null) return null;
		int[] l = new int[]{img[0].length, img[0][0].length};
		double t;
		double [][][] out = new double[4][l[0]][l[1]];
		for(int y = 0; y<l[0]; y++){
			for(int x = 0; x<l[1]; x++){
				t = img[0][y][x] + img[1][y][x] + img[2][y][x];
				out[0][y][x] = img[0][y][x]/t; 
				out[1][y][x] = img[1][y][x]/t;
				out[2][y][x] = img[2][y][x]/t;
				out[3][y][x] = t;
			}
		}
		return out;
	}
	
	public static java.awt.Color getColor(double d){
		double b, g, r, t;
		if(d>=0 && d<0.5){
			b = (0.5-d);
			g = d;
			r = 0;
			t = r+g+b;
			return new java.awt.Color( (int)(r/t*255), (int)(g/t*255), (int)(b/t*255));
		}
		if(d>=0.5 && d<=1.0){
			r = (d-0.5);
			g = (1.0-d);
			b = 0;
			t = r+g+b;
			return new java.awt.Color( (int)(r/t*255), (int)(g/t*255), (int)(b/t*255));
		}
		return new java.awt.Color(255, 255, 255);
	}
	
	public static java.awt.Color getNormalColorGrad(double d, double std){
		double b, g, r, t;
		double c1 = 1.0/Math.sqrt(2*std*Math.PI);
		double c2 = 1.0/(2.0*std);
		g = c1*Math.exp(-c2*Math.pow(d-.5,2));
		r = c1*Math.exp(-c2*Math.pow(d-1,2));
		b = c1*Math.exp(-c2*Math.pow(d,2));
		t = r+g+b;
		return new java.awt.Color( (int)(r/t*255), (int)(g/t*255), (int)(b/t*255));
	}
	
	public static void drawLine(byte[][][] img, int x0, int y0, int x1, int y1){
		double slope = (double)(x1-x0)/(double)(y1-y0);
		double intercept = x0-slope*y0;
		int minX, maxX, minY, maxY;
		if(x0>x1){
			minX = x1;
			maxX = x0;
		}else{
			minX = x0;
			maxX = x1;
		}
		if(y0>y1){
			minY = y1;
			maxY = y0;
		}else{
			minY = y0;
			maxY = y1;
		}
		if(minX<0 && maxX<0) return;
		if(minX>=img[0][0].length && maxX>=img[0][0].length) return;
		if(minY<0 && maxY<0) return;
		if(minY>=img[0].length && maxY>=img[0].length) return;
		if(minX<0) minX = 0;
		if(maxX>=img[0][0].length) maxX = img[0][0].length-1;
		if(minY<0) minX = 0;
		if(maxY>=img[0].length) maxY = img[0].length-1;
//		if(y0==y1){ //Horizontal line
//			for(int x=minX; x<=maxX; x++){
//				img[0][y0][x] = (byte) 255;
//				img[1][y0][x] = (byte) 0;
//				img[2][y0][x] = (byte) 0;
//			}
//		}else{ // Other lines
//			for(int y=minY; y<=maxY; y++){
//				int x = (int) (y*slope + intercept);
//				if(x<0 || x>=img[0][0].length) continue;
//				img[0][y][x] = (byte) 255;
//				img[1][y][x] = (byte) 0;
//				img[2][y][x] = (byte) 0;
//			}
//		}
		if(slope>-1 && slope<1){	// go by y index
			for(int y=minY; y<=maxY; y++){
				int x = (int) (Math.round(y*slope + intercept));
				if(x<0 || x>=img[0][0].length) continue;
				if(y<0 || y>=img[0].length) continue;
				img[0][y][x] = (byte) 255;
				img[1][y][x] = (byte) 0;
				img[2][y][x] = (byte) 0;
			}
		}else{ // go by x index
			slope = (double)(y1-y0)/(double)(x1-x0);
			intercept = y0-slope*x0;
			for(int x=minX; x<=maxX; x++){
				int y = (int) (Math.round(x*slope + intercept));
				if(x<0 || x>=img[0][0].length) continue;
				if(y<0 || y>=img[0].length) continue;
				img[0][y][x] = (byte) 255;
				img[1][y][x] = (byte) 0;
				img[2][y][x] = (byte) 0;
			}
		}
		return;
	}
}
