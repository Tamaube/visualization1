/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {
    
    private Volume volume = null;
    private GradientVolume gradients = null;
    private int currentMode;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private boolean responsive;
    
    public static final int MODE_SLICER = 0;
    public static final int MODE_MIP = 1;
    public static final int MODE_COMPOSITING = 2;
    public static final int MODE_2DTRANSFER = 3;
    
    private boolean shading = false;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
        this.currentMode = MODE_SLICER;
        responsive = false;
    }

    public boolean isShading() {
        return shading;
    }

    public void setShading(boolean shading) {
        this.shading = shading;
    }

    public boolean isResponsive() {
        return responsive;
    }

    public void setResponsive(boolean responsive) {
        this.responsive = responsive;
    }

    
    public int getCurrentMode() {
        return currentMode;
    }

    public void setCurrentMode(int currentMode) {
        if(currentMode == MODE_SLICER || currentMode == MODE_MIP ||
            currentMode == MODE_COMPOSITING || currentMode == MODE_2DTRANSFER) {
            this.currentMode = currentMode;
        }
    } 

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     
    short getVoxel(double[] coord) {
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }
    
    int getGreyPixel(TFColor voxelColor, int val, double max) {
        // Map the intensity to a grey value by linear scaling
        voxelColor.r = val/max;
        voxelColor.g = voxelColor.r;
        voxelColor.b = voxelColor.r;
        voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
       
        return getPixelInInt(voxelColor);
    }
    
    int getPixelInInt(TFColor voxelColor){
         // BufferedImage expects a pixel color packed as ARGB in an int
        int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
        int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
        int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
        int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
        
        return pixelColor;
    }
     

    void slicer(double[] viewMatrix) {
        clearImage();

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord = calculatePixelCoordinates(uVec, vVec, viewVec, volumeCenter, imageCenter, i, j, 0);

                int val = getVoxel(pixelCoord);
                image.setRGB(i, j, getGreyPixel(voxelColor, val, max));
            }
        }
    }
    
    //method for computing the coordinates of the pixel
    double[] calculatePixelCoordinates (double[] uVec, double[] vVec, double[] viewVec, double[] volumeCenter, int imageCenter, int i, int j, int k) {
        double[]  pixelCoord = new double[3];
        
        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + k*viewVec[0];
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + k*viewVec[1];
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + k*viewVec[2];
        
        return pixelCoord;
    }
    
    //linear interpolation
    double linearInterpolation(double v0, double v1, double x, double x0, double x1) {
        return (1 - (x-x0)/(x1-x0)) * v0 + (x-x0)/(x1-x0) * v1;
    }
    
    //trilinear interpolation
    double tripleLinearInterpolation(double [] pixel) {
        //get the two points
        int x0 = (int) Math.floor(pixel[0]);
        int y0 = (int) Math.floor(pixel[1]);
        int z0 = (int) Math.floor(pixel[2]);
        int x1 = (int) Math.ceil(pixel[0]);
        int y1 = (int) Math.ceil(pixel[1]);
        int z1 = (int) Math.ceil(pixel[2]);    
        
        //check that the two points are inside the volume dimensions as well
        if(!checkPixelInVolume(x0,y0,z0) || !checkPixelInVolume(x1,y1,z1)) {
            return 0;
        }
        
        //compute the corners
        double c000 = volume.getVoxel(x0, y0, z0);
        double c001 = volume.getVoxel(x0, y0, z1);
        double c010 = volume.getVoxel(x0, y1, z0);
        double c011 = volume.getVoxel(x0, y1, z1);
        double c100 = volume.getVoxel(x1, y0, z0);
        double c101 = volume.getVoxel(x1, y0, z1);
        double c110 = volume.getVoxel(x1, y1, z0);
        double c111 = volume.getVoxel(x1, y1, z1);
        
        //interpolate along x
        double c00 = linearInterpolation(c000, c100, pixel[0], x0, x1);
        double c01 = linearInterpolation(c001, c101, pixel[0], x0, x1);
        double c10 = linearInterpolation(c010, c110, pixel[0], x0, x1);
        double c11 = linearInterpolation(c011, c111, pixel[0], x0, x1);
        
        //interpolate along y
        double c0 = linearInterpolation(c00, c10, pixel[1], y0, y1);
        double c1 = linearInterpolation(c01, c11, pixel[1], y0, y1);
        
        //interpolate along z
        double c = linearInterpolation(c0, c1, pixel[2], z0, z1);
        
        return  c;
    }
    
    //method for checking if the coordinates of the pixel are within the volume's dimensions
     boolean checkPixelInVolume(double x, double y, double z) {
        return ( x>=0 && x<volume.getDimX() 
                && y>=0 && y<volume.getDimY()
                && z>=0 && z<volume.getDimZ());
    }
    
     void clearImage() {
         // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
     }
    
     //Maximum intensity projection
    void MIP(double[] viewMatrix) {
        this.clearImage();

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord;
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {                
                int val=0;
                
                //use index k to go along the ray
                int k = 0;
                int step = 1;
                //increase with 10 instead of 1 to render faster
                if(this.responsive) {step = 10;}
                
                pixelCoord = calculatePixelCoordinates(uVec, vVec, viewVec, volumeCenter, imageCenter, i, j, k);
                
                while(checkPixelInVolume(pixelCoord[0],pixelCoord[1],pixelCoord[2])) {
                    val = Math.max(val, (short)tripleLinearInterpolation(pixelCoord));                   
                    k += step;
                    pixelCoord = calculatePixelCoordinates(uVec, vVec, viewVec, volumeCenter, imageCenter, i, j, k);
                }
                
                image.setRGB(i, j, getGreyPixel(voxelColor, val, max));
            }
        }
    }
    
    //apply the new color as an overlay to old color 
    public TFColor applyColor(TFColor oldColor, TFColor newColor){
        TFColor result = new TFColor();
        result.r = oldColor.r * (1 - newColor.a) + newColor.r * newColor.a;
        result.g = oldColor.g * (1 - newColor.a) + newColor.g * newColor.a;
        result.b = oldColor.b * (1 - newColor.a) + newColor.b * newColor.a;
        
        return result;
    }

    void compositing(double[] viewMatrix) {
        this.clearImage();

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord;
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                //initialize color with black
                TFColor voxelColor = new TFColor(0,0,0,1); 
                
                int val=0;
                //use index k to go along the ray
                int k = 0;
                int step = 1;
                if(this.responsive) {step = 10;}
                
                pixelCoord = calculatePixelCoordinates(uVec, vVec, viewVec, volumeCenter, imageCenter, i, j, k);
                
                while(checkPixelInVolume(pixelCoord[0],pixelCoord[1],pixelCoord[2])) {
                    val = (short)tripleLinearInterpolation(pixelCoord);
                    TFColor newColor = tFunc.getColor(val);
                    voxelColor = applyColor(voxelColor, newColor);                    
                    k += step;
                    pixelCoord = calculatePixelCoordinates(uVec, vVec, viewVec, volumeCenter, imageCenter, i, j, k);
                }
                
                image.setRGB(i, j, getPixelInInt(voxelColor));
            }
        }
    }
    
    //apply shading according to the Simplified Phong model
    void applyShading(VoxelGradient voxelGradient, TFColor output, double[] viewVec){ 
        double k_amb = 0.1;
        double k_diff = 0.7;
        double k_spec = 0.2;
        double alpha = 10.0;
        
        double[] halfwayVec = new double[3];
        double[] voxelGrad = {(double) voxelGradient.x, (double) voxelGradient.y, (double) voxelGradient.z};     
         
        //V is viewVec normalized
        double[] Vvec ={viewVec[0], viewVec[1], viewVec[2]};
        VectorMath.scale(Vvec, 1/VectorMath.length(Vvec));
        
        //N is voxel gradients vector normalized
        VectorMath.scale(voxelGrad, 1/VectorMath.length(voxelGrad)); 
        
        //We assume L=V so, H=(2*V)/length(2*V)
        VectorMath.setVector(halfwayVec, Vvec[0],Vvec[1],Vvec[2]);
        VectorMath.scale(halfwayVec, 2);
        VectorMath.scale(halfwayVec, 1/VectorMath.length(halfwayVec));

        //L*N (=V*N, since L is equal to V)
        double dotProduct1 = VectorMath.dotproduct(Vvec, voxelGrad);
        //N*H
        double dotProduct2 = VectorMath.dotproduct(voxelGrad, halfwayVec);

        //formula only applies when the dot products are positive
        if (dotProduct1 >= 0.0 && dotProduct2 >= 0.0) {
            //tried using in the formula light white source, doesn't work, hence it is commented out
            //double[] w = new double[] {255.0,255.0,255.0};
            
            //get surface color
            double[] surfaceCol = new double[] {this.getTF2DPanel().triangleWidget.color.r,this.getTF2DPanel().triangleWidget.color.g,this.getTF2DPanel().triangleWidget.color.b};
            
            //apply formula
            //VectorMath.scale(surfaceCol, k_amb+k_diff*dotProduct1);
            VectorMath.scale(surfaceCol, k_diff*dotProduct1);
            //VectorMath.scale(w, k_amb);
            //VectorMath.addVector(surfaceCol, w);
            VectorMath.add(surfaceCol, k_spec*(Math.pow(dotProduct2,alpha)));

            //write the obtained new color to the output
            output.r = surfaceCol[0];
            output.g = surfaceCol[1];
            output.b = surfaceCol[2];
        }
    }
    
    //apply the new color and opacity
    public TFColor applyColorAndOpacity(TFColor oldColor, TFColor newColor){
        TFColor result = new TFColor();
        result.r = oldColor.r * (1 - newColor.a) + newColor.r * newColor.a;
        result.g = oldColor.g * (1 - newColor.a) + newColor.g * newColor.a;
        result.b = oldColor.b * (1 - newColor.a) + newColor.b * newColor.a;
        
        //apply levoy's relation (p.32)
        result.a = 1 - ((1 - oldColor.a) * (1-newColor.a) );
        
        return result;
    }
    
    //compute opactity
    //apply levoy's relation (p.32)
    //f(xi) = intensity
    //|▼f(xi)| = gradient
    //fv = baseIntensity
    //r = radius
    public double computeOpacity(short intensity, float gradient) {
        short baseIntensity = this.getTF2DPanel().triangleWidget.baseIntensity;
        double radius = this.getTF2DPanel().triangleWidget.radius;
        
        if(gradient == 0 && intensity == baseIntensity) {
            return 1;
        } else if(Math.abs(gradient) > 0 && 
                  intensity - radius * Math.abs(gradient) <= baseIntensity &&
                  intensity + radius * Math.abs(gradient) >= baseIntensity 
                ) {
            return 1 - 1/radius * (Math.abs(baseIntensity - intensity) / Math.abs(gradient));
        }
        return 0;            
    }

    void transfer2D(double[] viewMatrix) {
        this.clearImage();
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord;
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                //initial color taken from the surface color
                TFColor voxelColor =  new TFColor(this.getTF2DPanel().triangleWidget.color.r,this.getTF2DPanel().triangleWidget.color.g,this.getTF2DPanel().triangleWidget.color.b,this.getTF2DPanel().triangleWidget.color.a);
                //temporary copy
                TFColor tempColor =  new TFColor(voxelColor.r,voxelColor.g,voxelColor.b,voxelColor.a);
                
                short val=0;
                VoxelGradient gradient;
                
                //use index k to go along the ray
                int k = 0;
                int step = 1;
                if(this.responsive) {step = 10;}
                
                pixelCoord = calculatePixelCoordinates(uVec, vVec, viewVec, volumeCenter, imageCenter, i, j, k);
                
                //set opacity to 0
                voxelColor.a = 0;
                
                while(checkPixelInVolume(pixelCoord[0],pixelCoord[1],pixelCoord[2])) {                  
                    val = (short)tripleLinearInterpolation(pixelCoord);
                    gradient = this.gradients.getGradient((int)pixelCoord[0],(int) pixelCoord[1],(int) pixelCoord[2]);
                    
                    //get gradient based opacity and store in temporary color
                    tempColor.a = computeOpacity(val, gradient.mag);
                    
                    if(tempColor.a > 0) {                       
                        //check if shading is selected
                        if(shading) {
                            applyShading(gradient, tempColor, viewVec);
                        }
                        
                        //apply color
                        voxelColor = applyColorAndOpacity(voxelColor, tempColor);
                    }
                    
                    k += step;
                    pixelCoord = calculatePixelCoordinates(uVec, vVec, viewVec, volumeCenter, imageCenter, i, j, k);
                }
                
                image.setRGB(i, j, getPixelInInt(voxelColor));
            }
        }
    }
    
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {

        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();

        if(this.currentMode == MODE_SLICER) {
            slicer(viewMatrix);    
        } else if(this.currentMode == MODE_MIP) {
            MIP(viewMatrix);
        } else if(this.currentMode == MODE_COMPOSITING){
            compositing(viewMatrix);
        } else {
            transfer2D(viewMatrix);
        }
            
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}