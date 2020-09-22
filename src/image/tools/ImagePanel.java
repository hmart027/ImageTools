package image.tools;

import java.awt.Dimension;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;

import javax.swing.JPanel;

public class ImagePanel extends JPanel implements MouseListener, MouseWheelListener, MouseMotionListener, KeyListener{

	private static final long serialVersionUID = 7613874066609166472L;
	private Dimension dim;
	private BufferedImage imgOriginal = null;
	
	private int xOff, yOff, imgW, imgH, zxOff, zyOff;
	private float zoomLvl = 1, minZoom = 1, maxZoom = 1;
	private boolean enableZoom = true;
	
	private int pointX, pointY;
	private int lx, ly;
	private boolean dragging;

	public ImagePanel(Dimension dim){
		this.dim = dim;
		if(dim!=null)
			this.setPreferredSize(dim);
		else
			this.setPreferredSize(new Dimension(500,500));
		this.addMouseListener(this);
		this.addMouseWheelListener(this);
		this.addMouseMotionListener(this);
		this.addKeyListener(this);
	}
	
	public void setImage(BufferedImage img){
		setImage(img, true);
	}
	
	public void setImage(BufferedImage img, boolean recenter){
		setImage(img, recenter, true);
	}
	
	public void setImage(BufferedImage img, boolean recenter, boolean repaint){
		this.imgOriginal = img;
		if(recenter){
			setImgDims();
			computeMaxMinZomm();
		}
		if(repaint)
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
		this.zxOff = xOff;
		this.zyOff = yOff;
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
	
	@Override
	public void paint(java.awt.Graphics g){
		super.paintComponents(g);
		java.awt.Graphics2D g2 = (java.awt.Graphics2D)g;
		if(dim==null || !dim.equals(getSize())){
			this.dim = getSize();
			setImgDims();
			computeMaxMinZomm();
		}
		g2.setColor(java.awt.Color.BLACK);
		g2.fillRect(0, 0, dim.width, dim.height);
		if(imgOriginal!=null){
			g.drawImage(imgOriginal, zxOff, zyOff, (int)(zoomLvl*imgW), (int)(zoomLvl*imgH), null);
		}
	}
	
	public void zoom(float zoomLvl, int xc, int yc){
		if(zoomLvl>maxZoom){
			zoomLvl = maxZoom;
		}
		if(zoomLvl<minZoom){
			zoomLvl = minZoom;
		}
		float z = zoomLvl/this.zoomLvl;
		zxOff = (int) -(-zxOff*z+(z-1f)*(pointX));
		zyOff = (int) -(-zyOff*z+(z-1f)*(pointY));
		
		this.zoomLvl = zoomLvl;
		if(zoomLvl==minZoom){
			zxOff = xOff;
			zyOff = yOff;
		}
		repaint();
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		zoom(this.zoomLvl-e.getWheelRotation()*0.01f, 0, 0);
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		zxOff += e.getX() - lx;
		zyOff += e.getY() - ly;
		lx = e.getX();
		ly = e.getY();
		repaint();
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		pointX = e.getX();
		pointY = e.getY();
	}

	
	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	

	@Override
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	

	@Override
	public void keyReleased(KeyEvent e) {
		switch(e.getKeyCode()){
			default:
				System.out.println(e.getKeyCode());
		}
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		
	}

	@Override
	public void mousePressed(MouseEvent e) {
		dragging = true;
		lx = e.getX();
		ly = e.getY();
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		dragging = false;
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		
	}

	@Override
	public void mouseExited(MouseEvent e) {
		
	}
}