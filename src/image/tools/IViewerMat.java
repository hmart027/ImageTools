package image.tools;

import java.awt.Dimension;
//import java.awt.event.MouseEvent;
//import java.awt.event.MouseWheelEvent;
import java.awt.image.BufferedImage;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class IViewerMat extends JFrame {
	
	private static final long serialVersionUID = -1513114930536087394L;
	private HorizontalPanel[] verticalPanels;
	
	int heightGlobal = 800, widhtGlobal = 800;
		
	public IViewerMat(String title, int hCount, int vCount){
		this.setTitle(title);
		this.setBounds(0, 0, heightGlobal, widhtGlobal);
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		int w = widhtGlobal;
		int h = heightGlobal;
		
		int vd = h/vCount;
		verticalPanels = new HorizontalPanel[vCount];
		
		JPanel contentPane = new JPanel();
		contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
		
		for(int y=0; y<vCount; y++){
			verticalPanels[y] = new HorizontalPanel(w, vd-2, hCount);
			contentPane.add(verticalPanels[y]);
		}
				
//		viewer = new ImagePanel(new Dimension(w, h));
//		viewer.addMouseListener(new MouseListener());
//		viewer.addMouseWheelListener(new MouseWheelListener());
		
		this.setContentPane(contentPane);
		this.pack();
//		this.setResizable(false);
		this.setVisible(true);
		
	}
	
	public void setImage(BufferedImage img, int y, int x){
		this.verticalPanels[y].viewers[x].setImage(img);
	}
	
	protected class HorizontalPanel extends JPanel{
		private static final long serialVersionUID = 6322801222975317867L;
		public ImagePanel[] viewers;
		int width, height;
		
		public HorizontalPanel(int width, int height, int viewerC){
			this.width = width-20;
			this.height = height;
			viewers = new ImagePanel[viewerC];
			Dimension pDim = new Dimension(this.width/viewerC, this.height);
			for(int x=0; x<viewerC; x++){
				viewers[x] = new ImagePanel(pDim);
				this.add(viewers[x]);
			}
		}
		
		@Override
		public void setBounds(int x, int y, int width, int height){
			super.setBounds(x, y, width, height);
			this.width = width-20;
			this.height = height;
			Dimension pDim = new Dimension(this.width/viewers.length-2, this.height-2);
			for(int i=0; i<viewers.length; i++){
				viewers[i].setPreferredSize(pDim);
			}
		}
	}
	
	/*
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
		
	}*/
	
}
