package image.tools;

import java.awt.Dimension;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;

import javax.swing.JPanel;

public class ImagePanel extends JPanel implements MouseWheelListener, MouseMotionListener{

	private static final long serialVersionUID = 7613874066609166472L;
	private Dimension dim;
	private BufferedImage imgOriginal = null;
	
	private int xOff, yOff, imgW, imgH, zxOff, zyOff;
	private float zoomLvl = 1, minZoom = 1, maxZoom = 1;
	private boolean enableZoom = true;
	
	private int pointX, pointY;

	public ImagePanel(Dimension dim){
		this.dim = dim;
		if(dim!=null)
			this.setPreferredSize(dim);
		else
			this.setPreferredSize(new Dimension(500,500));
		this.addMouseWheelListener(this);
		this.addMouseMotionListener(this);
	}
	
	public void setImage(BufferedImage img){
		this.imgOriginal = img;
		setImgDims();
		computeMaxMinZomm();
		repaint();
	}
	
	private void setImgDims(){
		if(imgOriginal!=null){
			float dAR = (float)dim.width/(float)dim.height;
			float iAR = (float)imgOriginal.getWidth()/(float)imgOriginal.getHeight();
			if(iAR>dAR){
				imgW = dim.width;
				imgH = (int) (dim.width/iAR);
				xOff = 0;
				yOff = (dim.height - imgH)/2;
				if(yOff<0) yOff *= -1;
			}else{
				imgH = dim.height;
				imgW = (int) (dim.height*iAR);
				yOff = 0;
				xOff = (dim.width - imgW)/2;
				if(xOff<0) xOff *= -1;
			}
		}
	}
	
	private void computeMaxMinZomm(){
		this.zoomLvl = 1;
		this.minZoom = 1;
		this.maxZoom = 1;
		if(imgOriginal!=null){
			this.maxZoom = imgW;
			if(maxZoom<imgH)
				maxZoom = imgH;
		}
		this.zxOff = 0;
		this.zyOff = 0;
	}
	
	public BufferedImage getImage() {
		return this.imgOriginal;
	}
	
	@Override
	public void setPreferredSize(Dimension preferredSize) {
		super.setPreferredSize(preferredSize);
		this.dim = preferredSize;
		setImgDims();
		computeMaxMinZomm();
	}
	
	public void paint(java.awt.Graphics g){
		super.paintComponents(g);
		java.awt.Graphics2D g2 = (java.awt.Graphics2D)g;
		if(dim==null || !dim.equals(getSize())){
			this.dim = getSize();
			setImgDims();
		}
		g2.setColor(java.awt.Color.BLACK);
		g2.fillRect(0, 0, dim.width, dim.height);
		if(imgOriginal!=null){
			g.drawImage(imgOriginal, xOff-zxOff, yOff-zyOff, (int)(zoomLvl*imgW), (int)(zoomLvl*imgH), null);
		}
	}
	
	public void zoom(float zoomLvl, int xc, int yc){
		if(zoomLvl>maxZoom){
			zoomLvl = maxZoom;
		}
		if(zoomLvl<minZoom){
			zoomLvl = minZoom;
		}
//		zxOff = (int) ((zoomLvl-1)*imgW)/2;
//		zyOff = (int) ((zoomLvl-1)*imgH)/2;
		zxOff = (int) ((zoomLvl-1)*imgW/2f-pointX/zoomLvl-pointX*this.zoomLvl);
		zyOff = (int) ((zoomLvl-1)*imgH/2f-pointY/zoomLvl-pointY*this.zoomLvl);
		this.zoomLvl = zoomLvl;
		repaint();
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		zoom(this.zoomLvl-e.getWheelRotation()*0.01f, 0, 0);
	}

	@Override
	public void mouseDragged(MouseEvent e) {		
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		pointX = e.getX();
		pointY = e.getY();
	}
}