package img;

import java.util.ArrayList;
import java.util.List;

import camera.CameraModel;

public class AffineTransforms {
	
	static int dxW = 1;
	static int dyW = 1;
	
	public static List<Point3D> getRotationX(List<Point3D> points, double a){
		List<Point3D> out = new ArrayList<>();
		double aR = Math.toRadians(a);
		for(int i=0; i<points.size(); i++){
			Point3D p0 = points.get(i);
			Point3D p1 = new Point3D(
					p0.x, 
					+p0.y*Math.cos(aR)+p0.z*Math.sin(aR), 
					-p0.y*Math.sin(aR)+p0.z*Math.cos(aR));
			out.add(p1);
		}
		return out;
	}
	
	public static List<Point3D> getRotationY(List<Point3D> points, double a){
		List<Point3D> out = new ArrayList<>();
		double aR = Math.toRadians(a);
		for(int i=0; i<points.size(); i++){
			Point3D p0 = points.get(i);
			Point3D p1 = new Point3D( 
					p0.x*Math.cos(aR)-p0.z*Math.sin(aR), 
					p0.y, 
					p0.x*Math.sin(aR)+p0.z*Math.cos(aR));
			out.add(p1);
		}
		return out;
	}
	
	public static List<Point3D> getRotationZ(List<Point3D> points, double a){
		List<Point3D> out = new ArrayList<>();
		double aR = Math.toRadians(a);
		for(int i=0; i<points.size(); i++){
			Point3D p0 = points.get(i);
			Point3D p1 = new Point3D( 
					+p0.x*Math.cos(aR)+p0.y*Math.sin(aR), 
					-p0.x*Math.sin(aR)+p0.y*Math.cos(aR), 
					p0.z);
			out.add(p1);
		}
		return out;
	}

	public static List<Point3D> translate(List<Point3D> points, double dx, double dy, double dz){
		List<Point3D> out = new ArrayList<>();
		for(int i=0; i<points.size(); i++){
			Point3D p0 = points.get(i);
			Point3D p1 = new Point3D(p0.x+dx, p0.y+dy, p0.z+dz);
			out.add(p1);
		}
		return out;
	}
	
	public static byte[][] getPerspectiveTransform(List<Point3D> points, CameraModel c){
		int h = c.getHeight();
		int w = c.getWidth();
		int hC = h/2, wC = w/2;
		double f = c.getFocalLength();
		byte[][] out = new byte[h][w];	
		int x, y;
		for(Point3D p: points){
			x = round(wC + c.getHPpm()*p.x*f/(p.z+f));
			y = round(hC - c.getVPpm()*p.y*f/(p.z+f));
			for(int dy=-dyW; dy<=dyW; dy++){
				for(int dx=-dxW; dx<=dxW; dx++){
					if((x+dx)>=0 && (y+dy)>=0 && (y+dy)<h && (x+dx)<w)
						out[y+dy][x+dx] = (byte) 255;
				}
			}
		}
		return out;
	}
	
	private static int round(double d){
		int d1 = (int)d;
		int decimal = (int) ((d-d1)*10.0);
		if(decimal>=5)
			return d1+1;
		return d1;
	}

}
