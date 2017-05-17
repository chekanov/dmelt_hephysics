package hephysics.vec;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import net.jafama.FastMath;
/**
 * a HepLorentzVector (4-vector).
 * It can hold either space-coordinate and time (x,y,z,t) 
 * or 4-momentum (px,py,pz,energy).
 *  
 * @author Gary Bower (grb@slac.stanford.edu) and S.Chekanov (ANL)
 */

public class HepLorentzVector implements Serializable {
	static final long serialVersionUID = 1L;
	protected double t;
	protected double energy;
	protected Hep3Vector v;

	/**
	 * Default Lorentz vector constructor
	 */
	public HepLorentzVector() {
		this.t = 0.;
		this.energy=0;
		this.v = new Hep3Vector();
	}

	/**
	 * Define a Lorentz vector
	 * This can either be (x,y,z,time) or (px,py,pz,energy)
	 * @param x
	 *            X position (or Px)
	 * @param y
	 *            Y position (or Py)
	 * @param z
	 *            Z position (or Pz)
	 * @param t
	 *            time  or energy
	 */
	public HepLorentzVector(double x, double y, double z, double t) {
		this.t = t;
		this.energy=t;
		this.v = new Hep3Vector(x, y, z);
	}

	/**
	 * Define a Lorentz vector
	 * 
	 * @param x
	 *            3-vector position or momentum
	 * @param t
	 *            time  or energy
	 */
	public HepLorentzVector(double[] x, double t) {
		this.t = t;
		this.energy=t;
		this.v = new Hep3Vector(x);
	}

	/**
	 * Define a Lorentz vector
	 * 
	 * @param x
	 *            3-vector with position or 3-momentum
	 * @param t
	 *            time or energy
	 */
	public HepLorentzVector(float[] x, double t) {
		this.t = t;
		this.energy=t;
		this.v = new Hep3Vector(x);
	}

	/**
	 * Define Lorentz vector
	 * 
	 * @param x
	 *            3-vector with position or 3-momentum
	 * @param t
	 *            time or energy
	 */
	public HepLorentzVector(Hep3Vector v, double t) {
		this.t = t;
		this.v = v;
		this.energy=t;
	}

	/**
	 * Set 3-vector
	 * 
	 * @param v position 9x,y,z)
	 */
	public void setV3(Hep3Vector v) {
		this.v = v;
	}

	/**
	 * Set 3-vector position or Px,Py,Pz
	 * 
	 * @param x
	 *            X
	 * @param y
	 *            Y
	 * @param z
	 *            Z
	 */
	public void setV3(double x, double y, double z) {
		this.v.setV(x, y, z);
	}

	/**
	 * Get X position or (Px,Py,Pz)
	 * 
	 * @return X
	 */
	public double px() {
		return v.x();
	}

	/**
	 * Get Y position or Py
	 * 
	 * @return Y
	 */
	public double py() {
		return v.y();
	}

	/**
	 * Get Z position ot Pz
	 * 
	 * @return Z
	 */
	public double pz() {
		return v.z();
	}

	/**
	 * Get X position or Px
	 * 
	 * @return X
	 */
	public double x() {
		return v.x();
	}

	/**
	 * Get Y
	 * 
	 * @return Y
	 */
	public double y() {
		return v.y();
	}

	/**
	 * Get Z position or Pz
	 * 
	 * @return Z
	 */
	public double z() {
		return v.z();
	}
	
	/**
	 * Get energy
	 * 
	 * @return
	 */
	public double e() {
		return energy;
	}

	/**
	 * Get energy (as e())
	 * 
	 * @return
	 */
	public double getE() {
		return energy;
	}

	/**
	 * Set Px
	 * 
	 * @param px
	 */
	public void setPx(double px) {
		v.setX(px);
	}

	/**
	 * Set Py
	 * 
	 * @param py
	 */
	public void setPy(double py) {
		v.setY(py);
	}

	/**
	 * Set Pz
	 * 
	 * @param pz
	 */
	public void setPz(double pz) {
		v.setZ(pz);
	}

	
	
	
	/**
	 * Set X position.
	 * 
	 * @param x position.
	 */
	public void setX(double x) {
		v.setX(x);
	}

	/**
	 * Set Y position.
	 * 
	 * @param y Y position.
	 */
	public void setY(double y) {
		v.setY(y);
	}

	/**
	 * Set Z position.
	 * 
	 * @param z X position.
	 */
	public void setZ(double z) {
		v.setZ(z);
	}
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * Add 2 Lorentz vectors
	 * 
	 * @param parent
	 */
	public void add(HepLorentzVector parent) {
		v.setX(v.x() + parent.v3().x());
		v.setY(v.y() + parent.v3().y());
		v.setZ(v.z() + parent.v3().z());
		this.t = this.t + parent.getT();
		this.energy = this.energy + parent.getE();
	}

	/**
	 * Transverse mass: e()*e() - pz()*pz()
	 * 
	 * @return
	 */
	public double mt2() {
		return t() * t() - pz() * pz();
	}

	/**
	 * Magnitude squared: t()*t() - v.mag2() or e()*e() - v.mag2().
	 * This corresponds to the invariant mass squared
	 * 
	 * @return
	 */
	public double mag2() {

		return  (t() * t() - v.mag2());

	}

	/**
	 * Lorentz Boost using 3-vector
	 * 
	 * @param bx X
	 * @param by Y
	 * @param bz Z
	 */
	public void boost(double bx, double by, double bz) {

		double b2 = bx * bx + by * by + bz * bz;
		double gamma = 1.0 / FastMath.sqrt(1.0 - b2);
		double bp = bx * v.x() + by * v.y() + bz * v.z();
		double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

		v.setX(v.x() + gamma2 * bp * bx + gamma * bx * t());
		v.setY(v.y() + gamma2 * bp * by + gamma * by * t());
		v.setZ(v.z() + gamma2 * bp * bz + gamma * bz * t());
		setT(gamma * (t() + bp));

	}

	/**
	 * Set energy
	 * 
	 * @param e energy
	 */
	public void setE(double e) {
		this.energy = e;
	}

	/**
	 * Set time
	 * 
	 * @param t
	 *            time
	 */
	public void setT(double t) {
		this.t = t;
	}

	/**
	 * Get time
	 * 
	 * @return time
	 */
	public double t() {
		return t;
	}

	/**
	 * Get time (as t())
	 * 
	 * @return
	 */
	public double getT() {
		return t;
	}

	/**
	 * Get 3-vector
	 * 
	 * @return 3-vector
	 */
	public Hep3Vector v3() {
		return v;
	}

	
	/**
	 * Transverse energy squared.
	 */
	public double et2()  {
		   double pt2 = perp2();
		   return pt2 == 0 ? 0 : e()*e() * pt2/(pt2+v.z()*v.z());
		  }

	
	/**
	 * Transverse energy.
	 * @return
	 */
	 public double et()  {
	  double etet = et2();
      return e() < 0.0 ? -FastMath.sqrt(etet) : FastMath.sqrt(etet);
	 }

	 /**
	  * Rest mass squared -- same as m2()
	  * 
	  * @return
	  */
	 public  double restMass2() { return m2(); }
	 
	 
	 /**
	  * Transverse mass.
	  * @return
	  */
	 public double mt() { 
		 double mm = mt2();
		 return mm < 0.0 ? -FastMath.sqrt(-mm) : FastMath.sqrt(mm);	 
	 }

	 /**
	  * Same as m2(). 
	  * @return
	  */
	 public  double  invariantMass2() { return m2(); }

	 /**
	  * Invariant mass squared
	  * @return t**2 - x**2-y**2-z**2
	  */
	 public  double m2() { 
		 return energy*energy - mag2();
		}

	 /**
	  * Same as m().  If m2() is negative then -sqrt(-m2()) is returned.
	  * @return
	  */
	 public  double  invariantMass() { return m(); }
	 
	 
	 /**
	  * Invariant mass.  If m2() is negative then -sqrt(-m2()) is returned.
	  * @return invariant mass
	  */
	 public  double  m() { 
		 if (m2()>=0) return FastMath.sqrt(m2());
		 else  return -FastMath.sqrt(-m2());
		 
		 }
	 
	 
	 /**
	  * Is particle spacelike (i.e restMass2() smaller than 0)
	  * @return true if spacelike
	  */
	 public boolean isSpacelike() {
		  return restMass2() < 0;
		  }
	 
	 
		
	 /**
	  * is spacelike?
	  * @param epsilon precision
	  * @return
	  */
		public boolean isLightlike(double epsilon)  {
		 return FastMath.abs(restMass2()) < 2.0 * epsilon * energy * energy;
		 }
		 

	 
	 
	 /**
	  * Dot product
	  * @param v
	  * @param w
	  * @return
	  */
	 
		public double dot(HepLorentzVector v, HepLorentzVector w)
		   {
		      return v.t()*w.t() - dot(v.v3(),w.v3());
		   }

	
		
		/**
		 * Dot product for 3-vectors
		 * @param v
		 * @param w
		 * @return
		 */
		 private  double dot(Hep3Vector v, Hep3Vector w)
		   {
		      return v.x()*w.x() + v.y()*w.y() + v.z()*w.z();
		   }

	 
	 /**
	  * Add 2 vectors
	  * @param v
	  * @param w
	  * @return
	  */
		 public HepLorentzVector add(HepLorentzVector v, HepLorentzVector w)
		   {
		      return (new HepLorentzVector( add(v.v3(),w.v3()), v.t() + w.t() ));
		   }
		 	 
		 
		 /**
		  * Subtract 2 vectors
		  * @param v
		  * @param w
		  * @return
		  */
		   public  HepLorentzVector sub(HepLorentzVector v, HepLorentzVector w)
		   {
		      return (new HepLorentzVector(sub(v.v3(),w.v3()), v.t() - w.t()));
		   }
		   
		   /**
		    * Multiply  by a scaler
		    * @param scalar 
		    * @param v
		    * @return
		    */
		   public  HepLorentzVector mult( double scalar, HepLorentzVector v)
		   {
		      return (new HepLorentzVector(mult(scalar, v.v3()),scalar*v.t() ));
		   }
		   
		   /**
		    * Inverse
		    * @param v
		    * @return
		    */
		   public  HepLorentzVector neg(HepLorentzVector v)
		   {
		      return (new HepLorentzVector(neg(v.v3()),-v.t() ));
		   }
		   public  Hep3Vector neg(Hep3Vector v)
		   {
		      return new Hep3Vector(-v.x(), -v.y(), -v.z());
		   }
		 
		   
		   /**
		    * Add 3 vectors
		    * @param v
		    * @param w
		    * @return
		    */
		   private   Hep3Vector add(Hep3Vector v, Hep3Vector w)
		   {
		      return (new Hep3Vector(v.x() + w.x(), v.y() + w.y(), v.z() + w.z()));
		   }
		   private Hep3Vector sub(Hep3Vector v, Hep3Vector w)
		   {
		      return (new Hep3Vector(v.x() - w.x(), v.y() - w.y(), v.z() - w.z()));
		   }
		   private  Hep3Vector mult(double scalar, Hep3Vector v)
		   {
		      return (new Hep3Vector(scalar*v.x(), scalar*v.y(), scalar*v.z()));
		   }


		   
		   /**
		    * Return a center mass vector
		    * @param vecSet
		    * @return
		    */
		   public  Hep3Vector centerOfMass(List<HepLorentzVector> vecSet)
		   {
		      boolean empty = true;
		      boolean fourVecSet = false;
		      Hep3Vector cmVec = new Hep3Vector();
		      for (Iterator<HepLorentzVector> i = vecSet.iterator();  i.hasNext();)
		      {
		         if (empty == true)
		         {
		            empty = false;
		            Object e = i.next();
		            if ( e instanceof HepLorentzVector )
		            {
		               fourVecSet = true;
		               HepLorentzVector v = (HepLorentzVector) e;
		               cmVec = v.v3();
		            }
		            else
		            {
		               throw new RuntimeException("Element is not a 3- or 4-vector");
		            }
		            continue;

		         }
		         if ( fourVecSet )
		         {
		            try
		            {
		               HepLorentzVector v = (HepLorentzVector) i.next();
		               cmVec = add(cmVec, v.v3());
		            }
		            catch ( ClassCastException ex )
		            {
		               throw new RuntimeException(
		                  "Element of 4Vec enumeration is not a 4Vec.");
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

		   

		   
	
	/**
	 * cosTheta
	 * 
	 * @return
	 */
	public double cosTheta() {
		return v.cosTheta();
	}

	/**
	 * Phi
	 * 
	 * @return
	 */
	public double phi() {

		return v.phi();

	}

	/**
	 * Transverse momentum.
	 * 
	 * @return pt
	 */
	public double perp() {
		return v.perp();
	}

	/**
	 * Transverse momentum squared
	 * 
	 * @return pt*pt
	 */
	public double perp2() {
		return v.perp2();

	}

	
	

	
	
	/**
	 * Angle between 2 vectors
	 * 
	 * @param momentum
	 * @return
	 */
	public double angle(HepLorentzVector momentum) {
		return v.angle(momentum.v3());

	}

	/**
	 * Get theta angle
	 * 
	 * @return
	 */

	public double theta() {
		return v.theta();
	}

	
	/**
	    * Pseudorapidity. eta=-ln (tan (theta/2)) )
	    * @return
	    */
	   public double pseudoRapidity()  {
		  
		   return v.pseudoRapidity();
		  
		    
		   
	   }

	   
	   /**
	    * Rapidity. 0.5*log( (m+z)/(m-z) );
	    * @return
	    */
	 public double rapidity()  {
		  
		 double m = mag();
		 if (m> pz())  return 0.5*FastMath.log( (m+pz())/(m-pz()) );
		 else return -10e10;
		    
		   
	   }
	
	
	
	

	/**
	 * Get a pseudorapidity: eta=-ln (tan (theta/2)) )
	 * 
	 * @return pseudorapidity
	 */
	public double getEta() {
		return pseudoRapidity();
	}

	/**
	 * Get 3-vector
	 * 
	 * @return 3-vector
	 */
	public Hep3Vector getV3() {
		return v;
	}

	/**
	 * Magnitude sqrt(x**2+y**2+z**2)
	 * 
	 * @return
	 */
	public double mag() {
		return FastMath.sqrt(VecOp.dot(this.v, this.v));
	}

	/**
	 * Useful: for x1*x2+y1*y2+z1*z2
	 * 
	 * @param momentum
	 * @return
	 */
	public double skp(HepLorentzVector momentum) {

		return v.x() * momentum.v3().x() + v.y() * momentum.v3().y() + v.z()
				* momentum.v3().z();
	}

	/**
	 * Compre 2 vectors
	 */
	public boolean equals(Object obj) {
		if (obj instanceof HepLorentzVector) {
			HepLorentzVector that = (HepLorentzVector) obj;
			return v.equals(that.v3()) && that.t() == t;
		} else
			return false;
	}

	/**
	 * Make exact copy
	 * 
	 * @return new copy
	 */
	public HepLorentzVector copy() {
		HepLorentzVector tmp = new HepLorentzVector(px(), py(), pz(), t());
		tmp.setE(e());
		return tmp;

	}

	/**
	 * Hash code
	 */
	public int hashCode() {
		return v.hashCode() + (int) Double.doubleToRawLongBits(t);
	}

	/**
	 * Convert to string
	 */
	public String toString() {
		return VecOp.toString(this);
	}
}
