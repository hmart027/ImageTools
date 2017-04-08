package image.tools;

public class Kernels {

	/**
	 * Kernels
	 */
	// NOT Ready
	public static double[][] getDCTKernel(int n){
		
		double[][] out = new double[n][n];
		
		for(int y = 0; y<n; y++){
			for(int x = 0; x<n; x++){
				out[y][x] = Math.pow(-1, 0);
				if(x==0 && y==0)
					out[y][x] *=  1d/(double)n;
				else 
					out[y][x] *= 2d/(double)n;
			}
		}
		return out;
	}
	
	public static double[][] getHadamardKernel(int n){
		
		double[][] out = new double[n][n];
		
		for(int y = 0; y<n; y++){
			for(int x = 0; x<n; x++){
				float e = 0;
				for(int i = 0; i<n; i++){
					e += (float)(((x)>>i) & 0x01) * (float)(((y)>>(i)) & 0x01);
				}
				out[y][x] = Math.pow(-1, e);
			}
		}
		return out;
	}
	
	public static double[][] getWalshKernel(int n){
		
		double[][] out = new double[n][n];
		int n2 = (int) (Math.log(n)/Math.log(2));
		
		for(int y = 0; y<n; y++){
			for(int x = 0; x<n; x++){
				double p = 1;
				for(int i = 0; i<n2; i++){
					float e = (float)((x>>i) & 0x01) * (float)((y>>(n2-1-i)) & 0x01);
					p *= Math.pow(-1.0, e);
				}
				out[y][x] = p;
			}
		}
		return out;
	}
	
	public static double[][] getHaarKernel(int n){
		
		double[][] out = new double[n][n];
		int maxR = (int) (Math.log(n)/Math.log(2));
		int maxM = 0;
		int y = 0;
		
		double cx = 0;
		double x0 = 0;
		double x1 = 0;
		double x2 = 0;
		
		for (int x = 0; x < n; x++) {
			out[y][x] = 1/Math.pow(n, 0.5);
		}
		y++;
		for(int r = 0; r<maxR; r++){
			maxM = (int) Math.pow(2, r);
			for(int m = 1; m<=maxM; m++){
				x0 = (m - 1)/(double)maxM;
				x1 = (m - 0.5)/(double)maxM;
				x2 = m/(double)maxM;
				for (int x = 0; x < n; x++) {
					cx = (double)x / (double)n;
					if(x0<=cx && cx<x1){
						out[y][x] = Math.pow(2,r/2d)/Math.pow(n, 0.5);
					}else if(x1<=cx && cx<x2){
						out[y][x] = -Math.pow(2,r/2d)/Math.pow(n, 0.5);
					}else{
						out[y][x] = 0;
					}
				}
				y++;
			}
		}
		
		return out;
	}

	public static double[][] convolve(double[][] img, double kernel){
		int N = img.length;
		double[][] out = new double[N][N];
		
		
		return out;
	}
	/**
	 * END Kernels
	 */
	
	
	public static double[][] getVisualDCT(int n){
		
		double[][] out = new double[n*n][n*n];
		double c = Math.PI/(2*(double)n);
		
		for(int v = 0; v<n; v++){
			for(int u = 0; u<n; u++){
				for(int y = 0; y<n; y++){
					for(int x = 0; x<n; x++){
						out[v*n+y][u*n+x] = Math.cos((2*x+1)*u*c)*Math.cos((2*y+1)*v*c);
					}
				}
			}
		}
		
		return out;
	}
	
	public static double[][] getVisualDFTCosine(int n){
		
		double[][] out = new double[n*n][n*n];
		double c = 2*Math.PI/((double)n);
		
		for(int v = 0; v<n; v++){
			for(int u = 0; u<n; u++){
				for(int y = 0; y<n; y++){
					for(int x = 0; x<n; x++){
						out[v*n+y][u*n+x] = Math.cos((u*x+v*y)*c);
					}
				}
			}
		}
		
		return out;
	}
	
	public static double[][] getVisualDFTSine(int n){
		
		double[][] out = new double[n*n][n*n];
		double c = 2*Math.PI/((double)n);
		
		for(int v = 0; v<n; v++){
			for(int u = 0; u<n; u++){
				for(int y = 0; y<n; y++){
					for(int x = 0; x<n; x++){
						out[v*n+y][u*n+x] = Math.sin((u*x+v*y)*c);
					}
				}
			}
		}
		
		return out;
	}
}
