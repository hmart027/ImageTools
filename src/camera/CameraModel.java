package camera;

public class CameraModel {
	int width, height;
	double sensorHeight, sensorWidth;
	double vppm, hppm;
	double f;
	
	public int getHeight(){
		return this.height;
	}
	
	public int getWidth(){
		return this.width;
	}
	
	public double getFocalLength(){
		return this.f;
	}
	
	public double getVPpm(){
		return vppm;
	}
	
	public double getHPpm(){
		return hppm;
	}
}
