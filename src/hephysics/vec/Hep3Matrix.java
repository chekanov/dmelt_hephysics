package hephysics.vec;

import java.io.Serializable;
import hephysics.matrix.BasicMatrix;
import hephysics.matrix.Matrix;
import hephysics.matrix.MatrixOp;
import hephysics.matrix.MutableMatrix;
import hephysics.matrix.MatrixOp.InvalidMatrixException;
import net.jafama.FastMath;

/**
 * Hep 3x3 matrices
 * @see VecOp
 * @author Gary Bower and S.Chekanov
 */

public class Hep3Matrix extends  BasicMatrix {
   
   /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	 public Hep3Matrix(int nRows, int nCols)
	   {
	     super(nRows,nCols);
	   }
	   
	   /** Creates a new instance of BasicMatrix */
	   public Hep3Matrix(double[][] data)
	   {
		   super(data);
		  
	   }
	   
	   public Hep3Matrix()
	   {
		   super(3,3);
		   
	   }
	   
	   
	   public Hep3Matrix(double e11, double e12, double e13,
	                          double e21, double e22, double e23,
	                          double e31, double e32, double e33)
	   {
	    
		  super(e11, e12, e13,e21,e22,e23,e31,e32,e33);
	                          
		  
	   }
	   
	   
	   public Hep3Matrix(Matrix m) throws InvalidMatrixException
	   {
	    
		   super(m);
		   if (m.getNColumns() != 3 || m.getNRows() != 3) throw new InvalidMatrixException("Not 3x3 matrix");
	      m_dmat = new double[3][3];
	      m_dmat[0][0] = m.e(0,0);
	      m_dmat[0][1] = m.e(0,1);
	      m_dmat[0][2] = m.e(0,2);
	      m_dmat[1][0] = m.e(1,0);
	      m_dmat[1][1] = m.e(1,1);
	      m_dmat[1][2] = m.e(1,2);
	      m_dmat[2][0] = m.e(2,0);
	      m_dmat[2][1] = m.e(2,1);
	      m_dmat[2][2] = m.e(2,2);
	   }
	   /**
	    * Returns the (row, column) element
	    */
	   public double e(int row, int column)
	   {
	      return m_dmat[row][column];
	   }
	   /**
	    * Returns the determinent of the matrix.
	    */
	   public double det()
	   {
	      double cofact1 = m_dmat[1][1]*m_dmat[2][2] - m_dmat[1][2]*m_dmat[2][1];
	      double cofact2 = m_dmat[0][1]*m_dmat[2][2] - m_dmat[0][2]*m_dmat[2][1];
	      double cofact3 = m_dmat[0][1]*m_dmat[1][2] - m_dmat[0][2]*m_dmat[1][1];
	      return m_dmat[0][0]*cofact1 - m_dmat[1][0]*cofact2 + m_dmat[2][0]*cofact3;
	   }
	   /**
	    * Returns the trace of the matrix.
	    */
	   public double trace()
	   {
	      return m_dmat[0][0] + m_dmat[1][1] + m_dmat[2][2];
	   }
	   /**
	    * Sets the (row, column) element 
	    */
	   public void setElement( int row, int column, double value)
	   {
	      m_dmat[row][column] = value;
	   }
	   /**
	    * Defines a rotation matrix via Euler angles. A "passive" rotation
	    * matrix rotates the coordinate system, an "active" one rotates the
	    * vector(body). The angles are defined for a right handed coordinate
	    * system. They are defined by counterclockwise rotations about an
	    * axis by the right hand rule, ie, looking down the axis in the
	    * negative direction the transvers axes are seen to rotate
	    * counterclockwise. To define passive(active) angles first rotate the
	    * coordinates(body) about the z-axis by phi, then about the new x-axis
	    * by theta then about the new z-axis by psi.
	    * Angles in radians.
	    */
	   public void setPassiveEuler( double phi, double theta, double psi)
	   {
	      double cth = FastMath.cos(theta);
	      double sth = FastMath.sin(theta);
	      double cphi = FastMath.cos(phi);
	      double sphi = FastMath.sin(phi);
	      double cpsi = FastMath.cos(psi);
	      double spsi = FastMath.sin(psi);
	      m_dmat[0][0] =  cpsi*cphi - cth*sphi*spsi;
	      m_dmat[0][1] =  cpsi*sphi + cth*cphi*spsi;
	      m_dmat[0][2] =  spsi*sth;
	      m_dmat[1][0] = -spsi*cphi - cth*sphi*cpsi;
	      m_dmat[1][1] = -spsi*sphi + cth*cphi*cpsi;
	      m_dmat[1][2] =  cpsi*sth;
	      m_dmat[2][0] =  sth*sphi;
	      m_dmat[2][1] = -sth*cphi;
	      m_dmat[2][2] =  cth;
	   }
	   /**
	    * Defines a rotation matrix via Euler angles. A "passive" rotation
	    * matrix rotates the coordinate system, an "active" one rotates the
	    * vector(body). The angles are defined for a right handed coordinate
	    * system. They are defined by counterclockwise rotations about an
	    * axis by the right hand rule, ie, looking down the axis in the
	    * negative direction the transvers axes are seen to rotate
	    * counterclockwise. To define passive(active) angles first rotate the
	    * coordinates(body) about the z-axis by phi, then about the new x-axis
	    * by theta then about the new z-axis by psi.
	    * Angles in radians.
	    */
	   public void setActiveEuler( double phi, double theta, double psi)
	   {
	      double cth = FastMath.cos(theta);
	      double sth = FastMath.sin(theta);
	      double cphi = FastMath.cos(phi);
	      double sphi = FastMath.sin(phi);
	      double cpsi = FastMath.cos(psi);
	      double spsi = FastMath.sin(psi);
	      m_dmat[0][0] =  cpsi*cphi - cth*sphi*spsi;
	      m_dmat[1][0] =  cpsi*sphi + cth*cphi*spsi;
	      m_dmat[2][0] =  spsi*sth;
	      m_dmat[0][1] = -spsi*cphi - cth*sphi*cpsi;
	      m_dmat[1][1] = -spsi*sphi + cth*cphi*cpsi;
	      m_dmat[2][1] =  cpsi*sth;
	      m_dmat[0][2] =  sth*sphi;
	      m_dmat[1][2] = -sth*cphi;
	      m_dmat[2][2] =  cth;
	   }
	   
	   public static Hep3Matrix identity()
	   {
	      Hep3Matrix result = new Hep3Matrix();
	      result.m_dmat[0][0] = 1;
	      result.m_dmat[1][1] = 1;
	      result.m_dmat[2][2] = 1;
	      return result;
	   }
	   public String toString()
	   {
	      return VecOp.toString(this);
	   }

	   public int getNRows()
	   {
	      return 3;
	   }

	   public int getNColumns()
	   {
	      return 3;
	   }
	   
	   public void invert() throws MatrixOp.IndeterminateMatrixException
	   {
	      MatrixOp.inverse(this,this);
	   }

	   public void transpose()
	   {
	      double t = m_dmat[0][1];
	      m_dmat[0][1] = m_dmat[1][0];
	      m_dmat[1][0] = t;
	      
	      t = m_dmat[0][2];
	      m_dmat[0][2] = m_dmat[2][0];
	      m_dmat[2][0] = t;
	      
	      t = m_dmat[1][2];
	      m_dmat[1][2] = m_dmat[2][1];
	      m_dmat[2][1] = t;
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

}
