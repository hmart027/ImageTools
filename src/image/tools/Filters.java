package image.tools;

import java.util.Arrays;

import math2.Complex;

public class Filters {
	
	public static double[][] lowPass(double[][] freqDom, double r, boolean centered){
		int N = freqDom.length;
		int c = N/2;
		if(!centered) c = 0;
		double[][] out = new double[N][N];
		r = r*(N-c);
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++)
				if(Math.pow(Math.pow(x-c,2)+Math.pow(c-y, 2), 0.5)<=r)
					out[y][x] = freqDom[y][x];
				else
					out[y][x] = 0;
		return out;
	}
	
	public static Complex[][] lowPass(Complex[][] fft, double r){
		int N = fft.length;
		int c = N/2;
		Complex[][] out = new Complex[N][N];
		Complex zero = new Complex(0,0);
		r = r*N/2.0;
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++)
				if(Math.pow(Math.pow(x-c,2)+Math.pow(c-y, 2), 0.5)<=r)
					out[y][x] = fft[y][x];
				else
					out[y][x] = zero;
		return out;
	}
	
	public static double[][] highPass(double[][] fft, double r){
		int N = fft.length;
		int c = N/2;
		double[][] out = new double[N][N];
		r = r*2*N;
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++)
				if(Math.pow(Math.pow(x-c,2)+Math.pow(c-y, 2), 0.5)>=r)
					out[y][x] = fft[y][x];
				else
					out[y][x] = 0;
		return out;
	}
	
	public static Complex[][] highPass(Complex[][] fft, double r){
		int N = fft.length;
		int c = N/2;
		Complex[][] out = new Complex[N][N];
		Complex zero = new Complex(0,0);
		r = r*2*N;
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++)
				if(Math.pow(Math.pow(x-c,2)+Math.pow(c-y, 2), 0.5)>=r)
					out[y][x] = fft[y][x];
				else
					out[y][x] = zero;
		return out;
	}
	
	public static double[][] lowPassButterworth(double[][] fft, double coF, int order){
		int N = fft.length;
		int c = N/2;
		coF = coF*N/2.0;
		double[][] out = new double[N][N];
		double duv = 0;
		int exp = 2*order;
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++){
				duv = Math.pow(Math.pow(x-c,2)+Math.pow(c-y,2), 0.5);
				out[y][x] = fft[y][x]/(1+Math.pow(duv/coF,exp));
			}
		return out;
	}
	
	public static Complex[][] lowPassButterworth(Complex[][] fft, double coF, int order){
		int N = fft.length;
		int c = N/2;
		coF = coF*N/2.0;
		Complex[][] out = new Complex[N][N];
		double duv = 0;
		int exp = 2*order;
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++){
				duv = Math.pow(Math.pow(x-c,2)+Math.pow(c-y,2), 0.5);
				out[y][x] = new Complex(1.0/(1+Math.pow(duv/coF,exp)),0).times(fft[y][x]);
			}
		return out;
	}
	
	public static double[][] highPassButterworth(double[][] fft, double coF, int order){
		int N = fft.length;
		int c = N/2;
		coF = coF*N/2.0;
		double[][] out = new double[N][N];
		double duv = 0;
		int exp = 2*order;
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++){
				duv = Math.pow(Math.pow(x-c,2)+Math.pow(c-y,2), 0.5);
				out[y][x] = fft[y][x]/(1+Math.pow(coF/duv,exp));
			}
		return out;
	}
	
	public static Complex[][] highPassButterworth(Complex[][] fft, double coF, int order){
		int N = fft.length;
		int c = N/2;
		coF = coF*N/2.0;
		Complex[][] out = new Complex[N][N];
		double duv = 0;
		int exp = 2*order;
		for(int y = 0; y<N; y++)
			for(int x = 0; x<N; x++){
				duv = Math.pow(Math.pow(x-c,2)+Math.pow(c-y,2), 0.5);
				out[y][x] = new Complex(1.0/(1+Math.pow(coF/duv,exp)),0).times(fft[y][x]);
			}
		return out;
	}

	public static double[][] median(double[][] img, int kH, int kW){
		if(img == null) return null;
		if(kH<1 || kW<1) return null;
		int[] l = new int[]{img.length, img[0].length};
		double[][] out = new double[l[0]][l[1]];
		double[] k = new double[kH*kW];
		for(int y=0; y<l[0]-kH; y++){
			for(int x=0; x<l[1]-kW; x++){
				int i = 0;
				for(int yK=0; yK<kH; yK++){
					for(int xK=0; xK<kW; xK++){
						k[i++] = img[y+yK][x+xK];
					}
				}
				Arrays.sort(k);
				out[y+kH/2][x+kW/2] = k[k.length/2];
			}
		}
		return out;
	}

	public static double[][] gaussian(double[][] img, int size, double std){
		if(img == null) return null;
		if(size<1 || std<=0) return null;
		double[][] k = new double[size][size];
		double c1 = 1.0/(2.0*Math.PI*Math.pow(std, 2));
		double c2 = -c1*Math.PI;
		double off = size/2;
		for(int y=0; y<size; y++){
			for(int x=0; x<size; x++){
				k[y][x] = c1*Math.exp(c2*(Math.pow(off-y, 2)+Math.pow(x-off, 2)));
			}
		}
		return ITools.convolve(img, k);
	}
}
