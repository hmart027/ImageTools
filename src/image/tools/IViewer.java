package image.tools;

import java.awt.Dimension;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;

public class IViewer extends JFrame {
	
	private static final long serialVersionUID = -1513114930536087394L;
	private Dimension screen = null;
	private ImagePanel viewer;
	
	int heightGlobal = 800, widhtGlobal = 800;
	
	public IViewer(BufferedImage img){
		this("", img);
	}
	
	public IViewer(String title, BufferedImage img){		
		this.setTitle(title);
		this.setBounds(0, 0, heightGlobal, widhtGlobal);
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		this.screen = this.getToolkit().getScreenSize();
		
		int w = img.getWidth();
		int h = img.getHeight();
		float aspect = (float)w/(float)h;
		
		if(w>screen.getWidth()-20 || h>screen.getHeight()-20){
			float hResize = (float)screen.height/(float)h;
			float wResize = (float)screen.width/(float)w;
			
			if(hResize<=wResize){
				h = (int) (h*hResize - screen.height*0.1);
				w = (int) (h*aspect);
			}else{
				h = (int) (h/aspect);
				w = (int) (w*wResize - screen.width*0.1);
			}
		}
		
//		if(img.getWidth() < 512 || img.getHeight()<512){
			h=heightGlobal;w=widhtGlobal;
//		}
//		else{
//			w=img.getWidth();h=img.getHeight();
//		}
		
		viewer = new ImagePanel(new Dimension(w, h));
		viewer.setImage(img);
		viewer.addMouseListener(new MouseListener());
		viewer.addMouseWheelListener(new MouseWheelListener());
		
		this.getContentPane().add(viewer);
		this.pack();
//		this.setResizable(false);
		this.setVisible(true);
		
	}
	
	public void setImage(BufferedImage img){
		this.viewer.setImage(img);
	}
	
	private class MouseListener implements java.awt.event.MouseListener{

		@Override
		public void mouseClicked(MouseEvent e) {
//			System.out.println("Clicked");
		}

		@Override
		public void mousePressed(MouseEvent e) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void mouseReleased(MouseEvent e) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void mouseEntered(MouseEvent e) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void mouseExited(MouseEvent e) {
			// TODO Auto-generated method stub
			
		}
		
	}

	private class MouseWheelListener implements java.awt.event.MouseWheelListener{

		@SuppressWarnings("unused")
		@Override
		public void mouseWheelMoved(MouseWheelEvent e) {
			int rot = e.getWheelRotation();
//			System.out.println(rot);			
		}
		
	}
	
}
