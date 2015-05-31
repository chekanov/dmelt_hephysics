package hephysics.matrix;

import java.util.Formatter;

/**
 * Simple operations on matrices
 * @author tonyj
 */
public class MatrixOp
{
   private MatrixOp()
   {
   }
   /**
    * Invert matrix mIn and write it to matrix mOut.
    * This method allows both arguments to be the same, e.g. <code>inverse(this,this);</code>
    * This method currently only supports square matrices.
    */
   public static void inverse(Matrix mIn, MutableMatrix mOut) throws InvalidMatrixException
   {
      int order = mIn.getNRows();
      if (order != mIn.getNColumns()) throw new InvalidMatrixException("Matrix.inverse only supports square matrices");
      if (order != mOut.getNColumns() && order != mOut.getNRows()) throw new InvalidMatrixException("mOut must be same size as mIn");
      
      int[] ik = new int[order];
      int[] jk = new int[order];
      double[][] array = new double[order][order];
      for (int i=0;i<order;i++)
      {
         for (int j=0; j<order; j++)
         {
            array[i][j] = mIn.e(i,j);
         }
      }
      
      for (int k=0; k<order; k++)
      {
         // Find largest element in rest of matrix
         double amax = 0;
         for (int i=k; i<order; i++)
         {
            for (int j=k; j<order; j++)
            {
               if (Math.abs(array[i][j]) > Math.abs(amax))
               {
                  amax = array[i][j];
                  ik[k] = i;
                  jk[k] = j;
               }
            }
         }
         
         // Interchange rows and columns to put max in array[k][k]
         
         if (amax == 0) throw new IndeterminateMatrixException();
         
         {
            int i = ik[k];
            assert(k <= i);
            if (i > k)
            {
               for (int j=0; j<order; j++)
               {
                  double save = array[k][j];
                  array[k][j] = array[i][j];
                  array[i][j] = -save;
               }
            }
         }
         {
            int j = jk[k];
            assert(k <= j);
            if (j > k)
            {
               for (int i=0; i<order; i++)
               {
                  double save = array[i][k];
                  array[i][k] = array[i][j];
                  array[i][j] = -save;
               }
            }
         }
         
         // Accumulate elements of inverse matrix
         
         for (int i=0; i<order; i++)
         {
            if (i == k) continue;
            array[i][k] = -array[i][k]/amax;
         }
         for (int i=0; i<order; i++)
         {
            if (i == k) continue;
            for (int j=0; j<order; j++)
            {
               if (j == k) continue;
               array[i][j] += array[i][k]*array[k][j];
            }
         }
         for (int j=0; j<order; j++)
         {
            if (j == k) continue;
            array[k][j] = array[k][j]/amax;
         }
         array[k][k] = 1/amax;
      }
      
      // restore ordering of matrix
      
      for (int l=0; l<order; l++)
      {
         int k = order - l - 1;
         {
            int j = ik[k];
            if (j>k)
            {
               for (int i=0; i<order; i++)
               {
                  double save = array[i][k];
                  array[i][k] = -array[i][j];
                  array[i][j] = save;
               }
            }
         }
         {
            int i = jk[k];
            if (i>k)
            {
               for (int j=0; j<order; j++)
               {
                  double save = array[k][j];
                  array[k][j] = -array[i][j];
                  array[i][j] = save;
               }
            }
         }
      }
      for (int i=0;i<order;i++)
      {
         for (int j=0; j<order; j++)
         {
            mOut.setElement(i,j,array[i][j]);
         }
      }
   }
   public static String toString(Matrix m)
   {
      Formatter formatter = new Formatter();
      formatter.format("[");
      for (int i=0;;)
      {
         formatter.format("[");
         for (int j=0;;)
         {
            formatter.format("%12.5g",m.e(i,j));
            if (++j>=m.getNColumns()) break;
            formatter.format(",");
         }
         if (++i>=m.getNRows()) break;
         formatter.format("]\n ");
      }
      formatter.format("]");
      return formatter.out().toString();
   }
   // ToDo: Clean up the code here cut and pasted from invert().
   public static double det(Matrix mIn)
   {
      int order = mIn.getNRows();
      if (order != mIn.getNColumns()) throw new InvalidMatrixException("Matrix.det only supports square matrices");
      
      int[] ik = new int[order];
      int[] jk = new int[order];
      double[][] array = new double[order][order];
      for (int i=0;i<order;i++)
      {
         for (int j=0; j<order; j++)
         {
            array[i][j] = mIn.e(i,j);
         }
      }
      
      double det = 1;
      
      for (int k=0; k<order; k++)
      {
         // Find largest element array[i][k] in rest of matrix
         double amax = 0;
         for (int i=k; i<order; i++)
         {
            for (int j=k; j<order; j++)
            {
               if (Math.abs(array[i][j]) > Math.abs(amax))
               {
                  amax = array[i][j];
                  ik[k] = i;
                  jk[k] = j;
               }
            }
         }
         
         // Interchange rows and columns to put max in array[k][k]
         
         if (amax == 0) return 0;
         
         {
            int i = ik[k];
            assert (k <= i);
            if (i > k)
            {
               for (int j=0; j<order; j++)
               {
                  double save = array[k][j];
                  array[k][j] = array[i][j];
                  array[i][j] = -save;
               }
            }
         }
         {
            int j = jk[k];
            assert (k <= j);
            if (j > k)
            {
               for (int i=0; i<order; i++)
               {
                  double save = array[i][k];
                  array[i][k] = array[i][j];
                  array[i][j] = -save;
               }
            }
         }
         
         // Accumulate elements of inverse matrix
         
         for (int i=0; i<order; i++)
         {
            if (i == k) continue;
            array[i][k] = -array[i][k]/amax;
         }
         for (int i=0; i<order; i++)
         {
            if (i == k) continue;
            for (int j=0; j<order; j++)
            {
               if (j == k) continue;
               array[i][j] += array[i][k]*array[k][j];
            }
         }
         for (int j=0; j<order; j++)
         {
            if (j == k) continue;
            array[k][j] = array[k][j]/amax;
         }
         array[k][k] = 1/amax;
         det *= amax;
      }
      return det;
   }
   /**
    * Traspose matrix mIn and write it to matrix mOut.
    * This method allows both arguments to be the same, e.g. <code>transposed(this,this);</code>
    * This method currently only supports square matrices.
    */
   public static void transposed(Matrix mIn, MutableMatrix mOut)
   {
      int order = mIn.getNRows();
      if (order != mIn.getNColumns()) throw new InvalidMatrixException("Matrix.transposed only supports square matrices");
      if (order != mOut.getNColumns() && order != mOut.getNRows()) throw new InvalidMatrixException("mOut must be same size as mIn");
 
      for (int i=0; i<order; i++)
      {
         for (int j=0; j<i; j++)
         {
            double t1 = mIn.e(i,j); // In case mIn == mOut
            mOut.setElement(i,j,mIn.e(j,i));
            mOut.setElement(j,i,t1);
         }
         mOut.setElement(i,i,mIn.e(i,i));
      }
   }
   public static Matrix mult(Matrix m1, Matrix m2)
   {
      int nAdd = m1.getNColumns();
      if (nAdd != m2.getNRows()) throw new InvalidMatrixException("Incompatible matrices for multiplication");
      int nRows = m1.getNRows();
      int nCols = m2.getNColumns();
      BasicMatrix result = new BasicMatrix(nRows,nCols);
      for (int i=0; i<nRows; i++)
      {
         for (int j=0; j<nCols; j++)
         {
            double sum = 0;
            for (int k=0; k<nAdd; k++)
            {
               sum += m1.e(i,k)*m2.e(k,j);
            }
            result.setElement(i,j,sum);
         }
      }
      return result;
   }

   public static class IndeterminateMatrixException extends InvalidMatrixException
   {
      public IndeterminateMatrixException()
      {
         super("Matrix is indeterminate");
      }
   };
   public static class InvalidMatrixException extends RuntimeException
   {
      public InvalidMatrixException(String message)
      {
         super(message);
      }
   }
}
