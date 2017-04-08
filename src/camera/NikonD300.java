package camera;

public class NikonD300 extends CameraModel {

	public NikonD300(){
		this.width  = 4608; //4288
		this.height = 3072; //2848
		this.sensorWidth  = 0.0231;
		this.sensorHeight = 0.0154;
		this.f = 0.055;
		
		this.vppm = (double)height/sensorHeight;
		this.hppm = (double)width/sensorWidth;
	}
	 
	public double setFocalLenght(double f){
		double tf = this.f;
		this.f=f;
		return tf;
	}
}
