package hephysics.matrix;

import java.io.Serializable;

/**
 * A very simple matrix implementation
 * @author tonyj
 */
public class BasicMatrix implements MutableMatrix, Serializable
{
   static final long serialVersionUID = -3491275185124557222L;
   protected double[][] data;
   protected double[][] m_dmat;
   
   public BasicMatrix(int nRows, int nCols)
   {
      data = new double[nRows][nCols];
   }
   
   public BasicMatrix(double e11, double e12, double e13,
           double e21, double e22, double e23,
           double e31, double e32, double e33)
{


m_dmat = new double[3][3];
m_dmat[0][0] = e11;
m_dmat[0][1] = e12;
m_dmat[0][2] = e13;
m_dmat[1][0] = e21;
m_dmat[1][1] = e22;
m_dmat[1][2] = e23;
m_dmat[2][0] = e31;
m_dmat[2][1] = e32;
m_dmat[2][2] = e33;
}
   
   
   /** Creates a new instance of BasicMatrix */
   public BasicMatrix(double[][] data)
   {
      if (data.length <= 0) throw new IllegalArgumentException("Invalid data");
      int nCols = data[0].length;
      for (int i=0; i<data.length; i++)
      {
         if (data[i].length <= 0 || data[i].length != nCols) throw new IllegalArgumentException("Invalid data");
      }
      this.data = data;
   }
   public BasicMatrix(Matrix mIn)
   {
      int nRows = mIn.getNRows();
      int nCols = mIn.getNColumns();
      data = new double[nRows][nCols];
      for (int i=0; i<nRows; i++)
      {
         for (int j=0; j<nCols; j++)
         {
            data[i][j] = mIn.e(i,j);
         }
      }
   }

   public int getNRows()
   {
      return data.length;
   }

   public int getNColumns()
   {
      return data[0].length;
   }

   public double e(int row, int column)
   {
      return data[row][column];
   }
   
   public double det()
   {
      return MatrixOp.det(this);
   }
   
   public String toString()
   {
      return MatrixOp.toString(this);
   }

   public void setElement(int row, int column, double value)
   {
      data[row][column] = value;
   }

   public void invert() throws MatrixOp.IndeterminateMatrixException
   {
      MatrixOp.inverse(this,this);
   }
   
   public void transpose()
   {
      MatrixOp.transposed(this,this);
   }
}
