package hephysics.jet;

import java.util.Random;


import hephysics.particle.LParticle;

/*
 * An example of how to embed a JAIDA IPlotter into your own application.
 */
public class Main 
{
	
	  
   public static void main(String[] args)
   {


	   double P  = 0.0;  // GeV/c
	   
	   
	   /*
	  LParticle rho = new LParticle("rho",0.77540000);
	  LParticle pi0 = new LParticle("pi0",0.13497660);
	  LParticle pic = new LParticle("pi+",0.13957018);
	  LParticle gamma1 = new LParticle("gamma1",0);
	  LParticle gamma2 = new LParticle("gamma2",0);
	  */
	   
	      LParticle G = new LParticle("Graviton",2000);
		  LParticle top1 = new LParticle("top",170);
		  LParticle top2 = new LParticle("topBar",170);
		  LParticle W1 = new LParticle("W1",80);
		  LParticle W2 = new LParticle("W2",80);
		  LParticle b1 = new LParticle("b1",50);
		  LParticle b2 = new LParticle("b2",50);
		  LParticle q1 = new LParticle("q1",0.004);
		  LParticle q2 = new LParticle("q2",0.004);
		  LParticle q3 = new LParticle("q3",0.004);
		  LParticle q4 = new LParticle("q4",0.004);
		  
		  
		     System.out.println("Initial Graviton momenta:"+Double.toString(P));
		     System.out.println("Initial parameters:");
		     G.setV3( 0, 0, P );
		     
		     
		     G.print(); 
		     top1.print(); // System.out.println(pic.calcMass());
		     top2.print();  // System.out.println(pi0.calcMass());
		     W1.print(); // System.out.println(gamma1.calcMass());
		     W2.print(); // System.out.println(gamma2.calcMass());
		     b1.print(); // System.out.println(gamma1.calcMass());
		     b2.print(); // System.out.println(gamma2.calcMass());
		  
	  Random r= new Random();
	  for(int event = 0; event < 10; event++){

		     System.out.println("Event="+Integer.toString(event));
		     
		  
		     //***
		     //*** First decay:
		     //***
		     G.twoBodyDecay(top1, top2,true);
		     
		     double ran1= r.nextDouble(); 
		     double ran2= r.nextDouble();
		    // Isotropic angles to give a random direction 
		     double theta = Math.acos( 2.0*ran1 - 1.0 );
		     double phi   = 2.0 *Math.PI*ran2;

		     // rho 4-momentum components in lab frame
                     G.setThetaPhiP(theta,phi,P);   
		     // boost the pions
		     top1.boost( G );
		     top2.boost( G );

		   //***
		     //*** Second decay: pi0 --> gamma1 + gamma2 (pi0 rest frame)
		     //***
		     top1.twoBodyDecay(b1,W1, true);
		     // boost the gammas
		     b1.boost( top1 );
		     W1.boost( top1 );
		     
		     W1.twoBodyDecay(q1,q2, true);
		     // boost the gammas
		     q1.boost( W1 );
		     q2.boost( W1 );
		     
		     
		     top2.twoBodyDecay(b2,W2, true);
		     // boost the gammas
		     b2.boost( top2 );
		     W2.boost( top2 );
		     
		     W2.twoBodyDecay(q3,q4, true);
		     // boost the gammas
		     q3.boost( W2 );
		     q4.boost( W2 );
		     
		     G.print();
		     top1.print(); // System.out.println(pic.calcMass());
		     top2.print();  // System.out.println(pi0.calcMass());
		     W1.print(); // System.out.println(gamma1.calcMass());
		     W2.print(); // System.out.println(gamma2.calcMass());
		     b1.print(); // System.out.println(gamma1.calcMass());
		     b2.print(); // System.out.println(gamma2.calcMass());
                     q1.print(); 
                     q2.print(); 
                     q3.print(); 
                     q4.print();  
		    
		     double deg= (1.0 / 3.1415926535897931D) * 180D;
		     System.out.println( top1.angle(W1)*deg );
		     System.out.println( top1.angle(b1)*deg );
		     System.out.println( top1.angle(q1)*deg );
		     System.out.println( top1.angle(q2)*deg );
		     System.out.println( "Second decay");
		     System.out.println( top2.angle(W2)*deg );
		     System.out.println( top2.angle(b2)*deg );
		     System.out.println( top2.angle(q3)*deg );
		     System.out.println( top2.angle(q4)*deg );
 
                     System.out.println("Checking-----");   
                     LParticle XX = new LParticle("",0.0);
                     XX.add(q1); XX.add(q2); 
                     XX.add(q3); XX.add(q4); 
                     XX.add(b1); XX.add(b2); 
         	     System.out.println(XX.calcMass());
 
	  }
	  
 
  }

   
   
}


   
   
 
  
