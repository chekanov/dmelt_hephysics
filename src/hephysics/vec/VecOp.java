package hephysics.vec;

import hephysics.matrix.MatrixOp;
import java.util.Formatter;
import java.util.Iterator;
import java.util.List;

/**
 * Utility methods for dealing with 3 and 4 vectors.
 * @version $Id: VecOp.java 10570 2007-03-08 17:42:05Z tonyj $
 */
public class VecOp
{
   private VecOp()
   {
   }
   public static Hep3Vector add(Hep3Vector v, Hep3Vector w)
   {
      return (new Hep3Vector(v.x() + w.x(), v.y() + w.y(), v.z() + w.z()));
   }
   public static Hep3Vector sub(Hep3Vector v, Hep3Vector  w)
   {
      return (new Hep3Vector(v.x() - w.x(), v.y() - w.y(), v.z() - w.z()));
   }
   public static Hep3Vector mult(double scalar, Hep3Vector  v)
   {
      return (new Hep3Vector(scalar*v.x(), scalar*v.y(), scalar*v.z()));
   }
   public static Hep3Vector mult(Hep3Matrix m, Hep3Vector v)
   {
      double w1 = v.x()*m.e(0,0) + v.y()*m.e(0,1) + v.z()*m.e(0,2);
      double w2 = v.x()*m.e(1,0) + v.y()*m.e(1,1) + v.z()*m.e(1,2);
      double w3 = v.x()*m.e(2,0) + v.y()*m.e(2,1) + v.z()*m.e(2,2);
      return (new Hep3Vector(w1, w2, w3));
   }
   public static Hep3Matrix mult(Hep3Matrix m1, Hep3Matrix m2)
   {
      double e0 = m1.e(0,0) * m2.e(0,0) + m1.e(0,1) * m2.e(1,0) + m1.e(0,2) * m2.e(2,0);
      double e1 = m1.e(0,0) * m2.e(0,1) + m1.e(0,1) * m2.e(1,1) + m1.e(0,2) * m2.e(2,1);
      double e2 = m1.e(0,0) * m2.e(0,2) + m1.e(0,1) * m2.e(1,2) + m1.e(0,2) * m2.e(2,2);

      double e3 = m1.e(1,0) * m2.e(0,0) + m1.e(1,1) * m2.e(1,0) + m1.e(1,2) * m2.e(2,0);
      double e4 = m1.e(1,0) * m2.e(0,1) + m1.e(1,1) * m2.e(1,1) + m1.e(1,2) * m2.e(2,1);
      double e5 = m1.e(1,0) * m2.e(0,2) + m1.e(1,1) * m2.e(1,2) + m1.e(1,2) * m2.e(2,2);

      double e6 = m1.e(2,0) * m2.e(0,0) + m1.e(2,1) * m2.e(1,0) + m1.e(2,2) * m2.e(2,0);
      double e7 = m1.e(2,0) * m2.e(0,1) + m1.e(2,1) * m2.e(1,1) + m1.e(2,2) * m2.e(2,1);
      double e8 = m1.e(2,0) * m2.e(0,2) + m1.e(2,1) * m2.e(1,2) + m1.e(2,2) * m2.e(2,2);

      return new Hep3Matrix(e0,e1,e2,e3,e4,e5,e6,e7,e8);
   }
   public static Hep3Matrix mult(double scalar, Hep3Matrix m)
   {
      double e0 = m.e(0,0) * scalar;
      double e1 = m.e(0,0) * scalar;
      double e2 = m.e(0,0) * scalar;
      
      double e3 = m.e(0,1) * scalar;
      double e4 = m.e(0,1) * scalar;
      double e5 = m.e(0,1) * scalar;
      
      double e6 = m.e(0,2) * scalar;
      double e7 = m.e(0,2) * scalar;
      double e8 = m.e(0,2) * scalar;
      
      return new Hep3Matrix(e0,e1,e2,e3,e4,e5,e6,e7,e8);
   }
   public static Hep3Matrix inverse(Hep3Matrix m) throws MatrixOp.IndeterminateMatrixException
   {
      Hep3Matrix result = new Hep3Matrix();
      MatrixOp.inverse(m,result);
      return result;
   }
   static Hep3Matrix transposed(Hep3Matrix m)
   {
      Hep3Matrix result =  new Hep3Matrix();
      MatrixOp.transposed(m,result);
      return result;
   }
   public static  Hep3Vector neg(Hep3Vector v)
   {
      return new Hep3Vector(-v.x(), -v.y(), -v.z());
   }
   public static double dot(Hep3Vector v, Hep3Vector w)
   {
      return v.x()*w.x() + v.y()*w.y() + v.z()*w.z();
   }
   // ww, 07/31/00: cross product added
   public static Hep3Vector cross(Hep3Vector v, Hep3Vector w)
   {
      double u1 = v.y()*w.z() - v.z()*w.y();
      double u2 = v.z()*w.x() - v.x()*w.z();
      double u3 = v.x()*w.y() - v.y()*w.x();
      return new Hep3Vector(u1, u2, u3);
   }
   // ww, 08/01/00: unit vector added
   /**
    * returns (0,0,0) vector if input vector has length 0
    */
   public static Hep3Vector unit(Hep3Vector v)
   {
      double mag = v.mag();
      if ( mag != 0 )
         return mult(1./mag,v);
      else
         return (new Hep3Vector(0., 0., 0.));
   }
   
  
  
   /**
    * Boost fourVector with boostVector.
    *
    * Note, that beta=abs(boostVector) needs to be 0 < beta < 1.
    */
   public static HepLorentzVector boost(HepLorentzVector  fourVector, Hep3Vector boostVector)
   {
      double beta = boostVector.mag();
      
      if ( beta >= 1.0 )
         throw new RuntimeException("Boost beta >= 1.0 !");
      
      double gamma = 1./Math.sqrt(1.-beta*beta);
      
      double     t = fourVector.t();
      Hep3Vector v = fourVector.v3();
      
      double     tp = gamma*(t-dot(boostVector,v));
      Hep3Vector  vp = add(v,add(mult((gamma-1.)/(beta*beta)*dot(boostVector,v),boostVector),mult(-gamma*t,boostVector)));
      
      return new HepLorentzVector(vp,tp);
   }
   
   /**
    * Boost fourVector into system of refFourVector.
    */
   public static HepLorentzVector boost(HepLorentzVector  fourVector, HepLorentzVector  refFourVector)
   {
      Hep3Vector refVector = refFourVector.v3();
      Hep3Vector boostVector = new Hep3Vector(refVector.x(),refVector.y(),refVector.z());
      boostVector = mult(1./refFourVector.t(),boostVector);
      return boost(fourVector, boostVector);
   }
   
   public static double cosTheta(Hep3Vector vector)
   {
      return vector.z()/vector.mag();
   }
   public static double phi(Hep3Vector vector)
   {
      return Math.atan2(vector.y(),vector.x());
   }
   public static String toString(Hep3Vector  v)
   {
      Formatter formatter = new Formatter();
      formatter.format("[%9.4g,%9.4g,%9.4g]",v.x(),v.y(),v.z());
      return formatter.out().toString();
   }
   public static String toString(HepLorentzVector v)
   {
      Formatter formatter = new Formatter();
      formatter.format("[%12.5g,%12.5g,%12.5g,%12.5g]",v.v3().x(),v.v3().y(),v.v3().z(),v.t());
      return formatter.out().toString();
   }
   public static String toString(Hep3Matrix m)
   {
      return MatrixOp.toString(m);
   }
}