/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package util;

/**
 *
 * @author michel
 */
public class VectorMath {

    // assign coefficients c0..c2 to vector v
    public static void setVector(double[] v, double c0, double c1, double c2) {
        v[0] = c0;
        v[1] = c1;
        v[2] = c2;
    }

    // compute dotproduct of vectors v and w
    public static double dotproduct(double[] v, double[] w) {
        double r = 0;
        for (int i=0; i<3; i++) {
            r += v[i] * w[i];
        }
        return r;
    }

    // compute distance between vectors v and w
    public static double distance(double[] v, double[] w) {
        double[] tmp = new double[3];
        VectorMath.setVector(tmp, v[0]-w[0], v[1]-w[1], v[2]-w[2]);
        return Math.sqrt(VectorMath.dotproduct(tmp, tmp));
    }

    // compute dotproduct of v and w
    public static double[] crossproduct(double[] v, double[] w, double[] r) {
        r[0] = v[1] * w[2] - v[2] * w[1];
        r[1] = v[2] * w[0] - v[0] * w[2];
        r[2] = v[0] * w[1] - v[1] * w[0];
        return r;
    }
    
    // compute length of vector v
    public static double length(double[] v) {
        return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    
    //added to scale
    public static void scale(double[] vector, double factor){
        for(int i = 0; i < vector.length; ++i)
            vector[i] *= factor;
    }
    
    //find the most extreme value in the vector
    public static double extreme(double[] vector) {

        double extreme = 0; // The extreme value so far, including sign
        double absExtreme = 0; // The extreme value so far as a positive number

        double absCur; // The absolute value of the current value;

        for(double value : vector){
            absCur = Math.abs(value);
            if(absCur > absExtreme) {
                absExtreme = absCur;
                extreme = value;
            }
        }

        return extreme;
    }
    
    //added to normalize
    public static void normalize(double[] vector){
        scale(vector, 1/Math.abs(extreme(vector)));
    }
    
    public static void add(double[] u, double v){
        for(int i = 0 ; i < u.length; ++i)
            u[i] += v;
    }
}
