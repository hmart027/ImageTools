package image.tools;

import java.awt.Dimension;
import java.awt.image.BufferedImage;

import javax.swing.JPanel;

public class ImagePanel extends JPanel {

	private static final long serialVersionUID = 7613874066609166472L;
	private Dimension dim;
	private BufferedImage img = null;
	
	private int xOff, yOff, imgW, imgH;

	public ImagePanel(Dimension dim){
		this.dim = dim;
		if(dim!=null)
			this.setPreferredSize(dim);
		else
			this.setPreferredSize(new Dimension(500,500));
	}
	
	public void setImage(BufferedImage img){
		this.img = img;
		setImgDims();
		repaint();
	}
	
	private void setImgDims(){
		if(img!=null){
			float dAR = (float)dim.width/(float)dim.height;
			float iAR = (float)img.getWidth()/(float)img.getHeight();
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
	
	@Override
	public void setPreferredSize(Dimension preferredSize) {
		super.setPreferredSize(preferredSize);
		this.dim = preferredSize;
		setImgDims();
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
		if(img!=null){
			g.drawImage(img, xOff, yOff, imgW, imgH, null);
		}
	}
}