package image.tools;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.filechooser.FileNameExtensionFilter;

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
		
		int w = widhtGlobal;
		int h = heightGlobal;
		if(img!=null){
			w = img.getWidth();
			h = img.getHeight();
		}
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

		this.setMenuBar();
		this.getContentPane().add(viewer);
		this.pack();
//		this.setResizable(false);
		this.setVisible(true);
		
	}
	
	private void setMenuBar() {
		JMenuBar menuBar = new JMenuBar();
		JMenu fileMenu = new JMenu("File");
		JMenuItem saveItem = new JMenuItem("Save");
		saveItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if(viewer.getImage()!=null) {
					JFileChooser chooser = new JFileChooser();
					chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
					chooser.setFileFilter(new FileNameExtensionFilter("Images", new String[] {"jpg", "jpeg", "png", "tiff"}));
					if(chooser.showSaveDialog(null) == JFileChooser.APPROVE_OPTION){
						String path = chooser.getSelectedFile().getPath();
						String type = path.substring(path.lastIndexOf(".")+1, path.length());
					    try {
							ImageIO.write(viewer.getImage(), type, new File(path));
						} catch (IOException te) {
							System.out.println("Unable to save image");
						}
					}	
				}
			}
		});
		fileMenu.add(saveItem);
		menuBar.add(fileMenu);
		this.setJMenuBar(menuBar);
		
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
