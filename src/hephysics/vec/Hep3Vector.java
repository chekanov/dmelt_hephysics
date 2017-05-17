package hephysics.vec;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import net.jafama.FastMath;
/**
 * Basic implementation of a Hep3Vector
 * @author Gary Bower (grb@slac.stanford.edu)
 * @version $Id: BasicHep3Vector.java 9146 2006-10-16 19:22:42Z tonyj $
 */

public class Hep3Vector implements  Serializable
{
   static final long serialVersionUID = -52454965658870098L;
   protected double x;
   protected double y;
   protected double z;
   
   
   /**
    * New 3-vector
    */
   public Hep3Vector()
   {
      x = 0.;
      y = 0.;
      z = 0.;
   }
   
   /**
    * Define a 3-vector
    * @param dx
    * @param dy
    * @param dz
    */
   public Hep3Vector(double dx, double dy, double dz)
   {
      x = dx;
      y = dy;
      z = dz;
   }
   
   /**
    * Scalar product.
    * @param momentum
    * @return
    */
   
   double skp( Hep3Vector momentum)
   {
       return x * momentum.x() + y * momentum.y() + z * momentum.z();
   }

   /**
    * Magnitude
    * @return
    */
   public double abs()
   {
       return FastMath.sqrt(skp(this));
   }

   
   /**
    * New copy
    * @param p
    * @return
    */
   public Hep3Vector copy(Hep3Vector p) {
      return new Hep3Vector(p.x(), p.y(), p.z());
   }

   /**
    * cosTheta
    * @return
    */
   public double cosTheta() {
	    double ptot = mag();
	    return ptot == 0.0 ? 1.0 : z/ptot;
   }
   
   
   /**
    * Phi
    * @return
    */
   public double phi()
   {
	   return x == 0.0 && y == 0.0 ? 0.0 : FastMath.atan2(y,x);

   }
   
   
   /**
    * returns (0,0,0) vector if input vector has length 0
    */
   public Hep3Vector unit(Hep3Vector v)
   {
      double mag = v.mag();
      if ( mag != 0 )
         return mult(1./mag,v);
      else
         return (new Hep3Vector(0., 0., 0.));
   }

   /**
    * Invert
    * @param v
    * @return
    */
   public Hep3Vector neg(Hep3Vector v)
   {
      return new Hep3Vector(-v.x(), -v.y(), -v.z());
   }
 
  


   /**
    * Add 2 vectors and return a new object 
    * @param v
    * @param w
    * @return new object 
    */
   public Hep3Vector add(Hep3Vector v, Hep3Vector w)
   {
      return (new Hep3Vector(v.x() + w.x(), v.y() + w.y(), v.z() + w.z()));
   }
 

/**
    * Add a vector 
    * @param v
    * @return  
    */
   public void  add(Hep3Vector v)
   {
      setX(v.x() + x);
      setY(v.y() + y);
      setZ(v.z() + z);
   }

    /**
    * Subtract
    * @param v
    * @return
    */
   public void sub(Hep3Vector v)
   {
      setX(x-v.x());
      setY(y-v.y());
      setZ(z-v.z());

   }

 
   /**
    * Subtract
    * @param v
    * @param w
    * @return
    */
   public Hep3Vector sub(Hep3Vector v, Hep3Vector w)
   {
      return (new Hep3Vector(v.x() - w.x(), v.y() - w.y(), v.z() - w.z()));
   }
   
   
   /**
    * Multiply by a scaler 
    * @param scalar   
    * @param v
    * @return
    */
   public Hep3Vector mult(double scalar, Hep3Vector v)
   {
      return (new Hep3Vector(scalar*v.x(), scalar*v.y(), scalar*v.z()));
   }

  
 /**
    * Multiply by a scaler 
    * @param scalar   
    * @return
    */
   public void mult(double scalar)
   {

      setX(x*scalar);
      setY(y*scalar);
      setZ(z*scalar);
   }
 
   
   /**
    * Return a center of mass of a list of vectors
    * @param vecSet
    * @return
    */
   public Hep3Vector centerOfMass(List<Hep3Vector>  vecSet)
   {
      boolean empty = true;
      boolean threeVecSet = false;
      Hep3Vector cmVec = new Hep3Vector();
      for (Iterator<Hep3Vector>  i = vecSet.iterator();  i.hasNext();)
      {
         if (empty == true)
         {
            empty = false;
            Object e = i.next();
            if ( e instanceof Hep3Vector )
            {
               threeVecSet = true;
               cmVec = (Hep3Vector) e;
            }
            
            else
            {
               throw new RuntimeException("Element is not a 3- or 4-vector");
            }
            continue;

         }
       
         if ( threeVecSet )
         {
            try
            {
               cmVec = add(cmVec, (Hep3Vector) i.next());
            }
            catch ( ClassCastException ex )
            {
               throw new RuntimeException(
                  "Element of 3Vec enumeration is not a 3Vec object.");
            }
         }
      }
      if ( empty == false )
      {
         return cmVec;
      }
      else
      {
         throw new RuntimeException("CM:vector set is empty.");
      }
   }

  
   /**
    * Cross product of 2 vectors
    * @param v vector 1
    * @param w vector 2
    * @return cross product
    */
   public Hep3Vector cross(Hep3Vector v, Hep3Vector w)
   {
      double u1 = v.y()*w.z() - v.z()*w.y();
      double u2 = v.z()*w.x() - v.x()*w.z();
      double u3 = v.x()*w.y() - v.y()*w.x();
      return new Hep3Vector(u1, u2, u3);
   }

   
   /**
    * Dot product of 2 vectors
    * @param v vector 1
    * @param w vector 2
    * @return product
    */
   public double dot(Hep3Vector v, Hep3Vector w)
   {
      return v.x()*w.x() + v.y()*w.y() + v.z()*w.z();
   }


   
   /**
    * Transverse 
    * @return
    */
  public double perp() { return FastMath.sqrt(perp2()); }

  
  /**
   * Dot operation
   * @param p
   * @return
   */
  public double dot(Hep3Vector p)  {
	   return x*p.x() + y*p.y() + z*p.z();
  }
 
  
  /**
   * Transverse
   * @return
   */
   public double perp2()
   {
	   return x*x + y*y;

   }
   
   /**
    * Angle between 2 vectors
    * @param momentum
    * @return
    */
   public double angle(Hep3Vector  momentum)
   {
       if(abs() <= 0.0D || momentum.abs() <= 0.0D)
         return 0.0D;
       else
        return (double)(y * momentum.z() - z * momentum.y() <= 0.0D ? -1 : 1) * FastMath.acos(skp(momentum) / abs() / momentum.abs());
   }

   
  
   
   
   /**
    * Create a BasicHep3Vector from a double array
    * @param d An array {x,y,z}
    */
   public Hep3Vector(double[] d)
   {
      if (d.length != 3) throw new IllegalArgumentException("Illegal array length");
      x = d[0];
      y = d[1];
      z = d[2];
   }
   
   /**
    * Define a new vector
    * @param f
    */
   public Hep3Vector(float[] f)
   {
      if (f.length != 3) throw new IllegalArgumentException("Illegal array length");
      x = f[0];
      y = f[1];
      z = f[2];
   }
   
   /**
    * Set vector
    * @param dx
    * @param dy
    * @param dz
    */
   public void setV(double dx, double dy, double dz)
   {
      x = dx;
      y = dy;
      z = dz;
   }

   public double x()
   {
      return x;
   }
   public double y()
   {
      return y;
   }
   public double z()
   {
      return z;
   }
   
   /**
    * Set X component
    * @param x
    */
   public void setX(double x){
	   this.x=x;
   }
   
   
   /**
    * Boost fourVector with boostVector.
    *
    * Note, that beta=abs(boostVector) needs to be 0 < beta < 1.
    */
   public  HepLorentzVector boost(HepLorentzVector fourVector, Hep3Vector boostVector)
   {
      double beta = boostVector.mag();

      if ( beta >= 1.0 )
         throw new RuntimeException("Boost beta >= 1.0 !");

      double gamma = 1./FastMath.sqrt(1.-beta*beta);

      double     t = fourVector.t();
      Hep3Vector v = fourVector.v3();

      double     tp = gamma*(t-dot(boostVector,v));
      Hep3Vector vp = add(v,add(mult((gamma-1.)/(beta*beta)*dot(boostVector,v),boostVector),mult(-gamma*t,boostVector)));

      return new HepLorentzVector(vp,tp);
   }

   /**
    * Boost fourVector into system of refFourVector.
    */
   public HepLorentzVector boost(HepLorentzVector fourVector, HepLorentzVector refFourVector)
   {
      Hep3Vector refVector = refFourVector.v3();
      Hep3Vector boostVector = new Hep3Vector(refVector.x(),refVector.y(),refVector.z());
      boostVector = mult(1./refFourVector.t(),boostVector);
      return boost(fourVector, boostVector);
   }

   

   public Hep3Vector mult(Hep3Matrix m, Hep3Vector v)
   {
      double w1 = v.x()*m.e(0,0) + v.y()*m.e(0,1) + v.z()*m.e(0,2);
      double w2 = v.x()*m.e(1,0) + v.y()*m.e(1,1) + v.z()*m.e(1,2);
      double w3 = v.x()*m.e(2,0) + v.y()*m.e(2,1) + v.z()*m.e(2,2);
      return (new Hep3Vector(w1, w2, w3));
   }

   /**
    * Y component
    * @param y
    */
   public void setY(double y){
	   this.y=y;
   }
   
   /**
    * Set Z component
    * @param z
    */
   public void setZ(double z){
	   this.z=z;
   }
   
   public double mag2() { return x*x + y*y + z*z; }
   
   
   /**
    * angle theta (=atan2(perp(),z)
    * @return
    */
   public double theta() {
	   return x == 0.0 && y == 0.0 && z == 0.0 ? 0.0 : FastMath.atan2(perp(),z);
   }
   
   
   /**
    * Cosine Theta
    * @param vector
    * @return
    */
   public double cosTheta(Hep3Vector vector)
   {
      return vector.z()/vector.mag();
   }

  
   
   /**
    * Set phi
    * @param ph phi in radians
    */
   public void setPhi(double ph) {
	    double xy   = perp();
	    setX(xy*FastMath.cos(ph));
	    setY(xy*FastMath.sin(ph));
   }
   
   /**
    * Pseudorapidity. eta=-ln (tan (theta/2)) )
    * @return
    */
   public double pseudoRapidity()  {
	  
	   
	   double cosTheta = cosTheta();
	   if (cosTheta*cosTheta < 1) return -0.5* FastMath.log( (1.0-cosTheta)/(1.0+cosTheta) );
	   else        return -10e10;
	    
	   
   }

   
   /**
    * Rapidity. 0.5*log( (m+z)/(m-z) );
    * @return
    */
 public double rapidity()  {
	  
	 double m = mag();
	 if (mag()> z)  return 0.5*FastMath.log( (m+z)/(m-z) );
	 else return -10e10;
	    
	   
   }
   
   
   
   public double getEta()  { return pseudoRapidity();}

   
   
   public double mag()
   {
      return FastMath.sqrt(x*x + y*y + z*z);
   } 
   public double magnitudeSquared()
   {
      return x*x + y*y + z*z;
   }
   public double[] v()
   {
      return new double[] { x, y, z };
   } 

   public boolean equals(Object obj)
   {
      if (obj instanceof Hep3Vector)
      {
         Hep3Vector that = (Hep3Vector) obj;
         return x == that.x() && y == that.y() && z == that.z();
      }
      else return false;
   }

   public String toString()
   {
      return VecOp.toString(this);
   }

   public int hashCode()
   {
      return (int) (Double.doubleToLongBits(x) +
                    Double.doubleToLongBits(y) +
                    Double.doubleToLongBits(z));
   }








}
