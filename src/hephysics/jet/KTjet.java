package hephysics.jet;

import hephysics.jet.KTjet;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;
import java.text.*;

/**
 * Longitudinally-invariant kT, anti-KT and  Cambridge/Aachen clustering algorithms (light-weight implementation).  
 * This class uses float values for fast computation and {@link hephysics.jet.ParticleF} with float definitions. The algorithm uses phi-pseudorapidity to define distances, similar to the Tevatron approach. <p></p> To speed-up calculations and to have a low memory footprint, it uses {@link hephysics.jet.ParticleF} with float definitions. Use slower {@link hephysics.jet.SCJet} class for with  double types and rapidity, similar to the LHC approach. This class uses E-scheme to combine particles (p1+p2).
 * More details in http://arxiv.org/pdf/hep-ph/0210022v1.pdf.
 * 
 * @author S.Chekanov
 * 
 */
public class KTjet {

	private int recom = 1;
	private static float R;
	private static float R2;
	private int[] is_consider;
	private float[] ktdistance1;
	private float[][] ktdistance12;
	private ArrayList<ParticleF> jets;
	static private final float PI2 = (float)(Math.PI * 2);
	private boolean debug=false;
	private double minpt=0;
	private static int mode=1;
	private DecimalFormat formatter = new DecimalFormat("#.#####");

	/**
	 * Initialize calculations of the longitudinally invariant kT algorithm in inclusive mode. 
	 * Jet can be clustered using Cambridge/Aachen or anti-kT approaches, depending on the "mode" parameter. 
	 * The distance parameters are Eta (pseudorapidity) and Phi. 
	 * 
	 * @param R
	 *            distance measure
	 * @param recom
	 *            recombination scheme.<br>
	 *            1: The E-scheme Simple 4-vector addition. <br>
	 *            2: The pT-scheme. <br>
	 *            3: The pT^2 scheme. <br>
	 *            Currently only E-scheme is implemented.
	        * @param mode 
	        *          clustering mode dij=min(kT_i^{2* mode},kT_j^{2* mode})). <br>
	        *          mode=1 means inclusive KT jet algorithm <br> 
	        *          mode=0 means Cambridge/Aachen jet algorithm <br> 
	        *          mode=-1 means anti-KT jet algorithm <br> 
	        * @param minpt
	        *            min pT for final jets.
	 */
	public KTjet(double R, int recom, int mode, double minpt) {
		this.R = (float)R;
		this.R2 = (float)(R * R);
		this.recom = recom;
		this.debug=false;
		this.minpt=minpt;
		this.mode=mode;
		DecimalFormat formatter = new DecimalFormat("#0.00");
		String rs=formatter.format(this.R);
		System.out.println("Initialization of Java jet algorithm. S.Chekanov (ANL)");
		System.out.println("Inclusive mode using the E-scheme recombination and R="+rs);
                System.out.println("Distance is defined in pseudo-rapidity, phi, Et");
		if (mode==1) System.out.println("Longitudinally invariant kT algorithm");
		else if (mode==0) System.out.println("Cambridge/Aachen algorithm" );
		else if (mode==-1) System.out.println("Longitudinally invariant anti-kT algorithm");
		else  System.out.println("Not correct mode:  Fallback to the inclusive kT algorithm using E-scheme and R="+rs);
	}


	/** Initialize calculations of the KT algorithm. Meaningful values are R=0.2- 1.
	* Jets are clustered in Eta (pseudorapidity) and Phi space. 
	* 
	* @param R
	*            distance measure
	* @param recom
	*            recombination scheme.<br>
	*            1: The E-scheme Simple 4-vector addition. <br>
	*            2: The pT-scheme. <br>
	*            3: The pT^2 scheme. <br>
	*            Currently only E-scheme is implemented.
	* @param minpt
	*            min pT for final jets.
	*/
	public KTjet(double R, int recom, double minpt) {
		this(R,recom,1, minpt);
	}


	/** Initialize calculations of the KT algorithm. Meaningful values are R=0.2- 1.
	* Jets are clustered in Eta (pseudorapidity) and Phi space. The The E-scheme with 4-vector addition is used. 
	* 
	* @param R
	*            distance measure
	* @param minpt
	*            min pT for final jets.
	*/
	public KTjet(double R, double minpt) {
		this(R,1, 1, minpt);
	}




	/**
	  * Run the jet algorithm using the list of particles 
	  * 
	  * @param list
	  *            list with particles
	  * @return final jets without sorting.
	  */
	public  ArrayList<ParticleF>  buildJets(ArrayList<ParticleF> list) {

		jets = new ArrayList<ParticleF>();
		int size = list.size();

		long startTime = 0;
		if (debug) startTime=System.currentTimeMillis();

		int j1 = -1;
		int j2 = -1;
		double min12 = Double.MAX_VALUE;

		ktdistance1 = new float[size];
		is_consider = new int[size];
		for (int m = 0; m < size; m++) {
			is_consider[m] = 1;
			ParticleF p1 = (ParticleF) list.get(m);
			ktdistance1[m] = getKtDistance1(p1);
		}


		ktdistance12 = new float[size][size];
		for (int i = 0; i < size - 1; i++) {
			ParticleF p1 = (ParticleF) list.get(i);
			for (int j = i + 1; j < size; j++) {
				ParticleF p2 = (ParticleF) list.get(j);
				ktdistance12[i][j] = getKtDistance12(p1, p2);
				//System.out.println(ktdistance12[i][j]);
			}
		}

		if (debug) {
			long stopTime = System.currentTimeMillis();
			long runTime = stopTime - startTime;
			System.out.println("--->  Run time after creating initial pair distances (ms): " + runTime);
		}

	
                boolean merged=false; 
		int Nstep = size;
		while (Nstep > 0) {

			min12 = Double.MAX_VALUE;
			// this is after reseting to a new jet
			if (!merged) {
				for (int i = 0; i < size-1; i++) {
					if (is_consider[i]<=0) continue;
					for (int j = i+1; j < size; j++) {
						if (is_consider[j]<=0) continue;
						if (ktdistance12[i][j] < min12) {
							min12 = ktdistance12[i][j];
							j1 = i;
							j2 = j;
						}
					}
				}
			} else {
				// find another minimum around this jet  when j1>0
				for (int j = 0; j < size; j++) {
					if (is_consider[j]<=0 || j == j1) continue;
					if (ktdistance12[j1][j] < min12) {
						min12 = ktdistance12[j1][j];
						j1 = j1;
						j2 = j;
					}
				}


			} // end of min finding


                        if (merged==false && Nstep==1) break;

			// find min distance to the beam
			double min1 = ktdistance1[j1];
			if (ktdistance1[j2]<min1) {min1 = ktdistance1[j2];};


                    // make the decision about this particle
                        merged=false;
                        if (min12<min1) merged=true;

                        if (merged) {
                                ParticleF p1 = (ParticleF) list.get(j1);
                                ParticleF p2 = (ParticleF) list.get(j2);
                                if (j1 != j2) p1.add(p2); 
                                Nstep--;
                                list.set(j1, p1); // replace with p1+p2
                                is_consider[j2] = 0;
                                is_consider[j1]=is_consider[j1]+1;
                                // recalculate distance for this particle
                                ktdistance1[j1]  = getKtDistance1(p1);
                                for (int i = 0; i < size; i++) {
                                        if (is_consider[i]<=0 || i==j1)  continue;
                                        ParticleF pp1 = (ParticleF)list.get(i);
                                        ktdistance12[j1][i] = getKtDistance12(p1, pp1);
                                }

                        }


                        if (!merged) {   // add this to the jet
                                is_consider[j1] = -1;
                                ParticleF pj = (ParticleF) list.get(j1);
                                Nstep--;
                                if (pj.getEt() > minpt) {
                                        jets.add(pj); // fill jets
                                }
                        }


                        // end loop
                }
	


        if (debug) {
        // attempt to deal with unmeargable particle
        int ins=-1;
        for (int i = 0; i < size; i++)
        if (is_consider[i]==1) {ins=i;};

        if (ins>-1) {
                ParticleF p2 = list.get(ins);
                if (debug) System.out.println("Unmerged particle id="+Integer.toString(ins));
                min12 = Double.MAX_VALUE;
                for (int j = 0; j < jets.size(); j++) {
                        ParticleF lp = jets.get(j);
                        float d=getDistance(p2, lp);
                        if (d<min12) { j1=j; min12 =d; };
                }
                if (debug)  System.out.println("Distance R to closest jet="+Double.toString(min12));
                if (min12<R) {
                        if (debug)  System.out.println(" --> Particle merged");
                        ParticleF lp = jets.get(j1);
                        lp.add(p2);
                        is_consider[ins] = 0;
                }
        }

                // sanity test. All particles were merged?
                int nn=0; ins=-1;
                for (int i = 0; i < size; i++)
                if (is_consider[i]==1) {nn++; ins=i;};
                if (nn != 0)   System.out.println( "--> WARNING: particle with ID="+ Integer.toString(ins)+" unmerged");


                        long stopTime2 = System.currentTimeMillis();
                        long runTime = stopTime2 - startTime;
                        System.out.println("  --> Final time for calculation (ms): "
                                        + runTime);
                        System.out.println("  --> Nr of jets : " + jets.size());

        } // end debug mode


		if (debug) {
			long stopTime2 = System.currentTimeMillis();
			long runTime = stopTime2 - startTime;
			System.out.println("  --> Final time for calculation (ms): " + runTime);
		}

		is_consider=null;
		ktdistance12=null;
		ktdistance1=null;
		return jets;

	}




	/**
	* Get jets after  sorting in jet pT. Run  buildJets before calling this method. 
	* 
	* @return list with sorted jets
	*/

	public ArrayList<ParticleF> getJetsSorted() {

		Collections.sort(jets);

		return jets;
	}



        /**
         * Print the kT jets for debugging.
         */
        public void printJets() {

                ArrayList<ParticleF> sjets=getJetsSorted();

                System.out.println("# Nr of jets=" + Integer.toString(sjets.size()));
                System.out.format("%5s %14s %14s %14s \n","jet #", "Eta", "phi", "Et");
                for (int i = 0; i < sjets.size(); i++) {
                        ParticleF lp = sjets.get(i);
                        double phi = lp.getPhi();
                        if (phi<0) phi=PI2+phi;
                        String s1=String.format("%15.8f", lp.getEta());
                        String s2=String.format("%15.8f", phi);
                        String s3=String.format("%15.8f",lp.getEt());
                        System.out.format("%5s%15s%15s%15s\n",Integer.toString(i),s1,s2,s3);
                }

        }


	/**
	* Print the kT jets for debugging to a string.
	* @return String representing a jet 
	*/
	public String toString() {
		ArrayList<ParticleF> sjets=getJetsSorted();
		String tmp="# Nr of jets=" + Integer.toString(sjets.size())+"\n";
		for (int i = 0; i < sjets.size(); i++) {
			ParticleF lp = sjets.get(i);
			String spx=formatter.format(lp.getEta());
			String spy=formatter.format(lp.getPhi());
			String spz=formatter.format(lp.getEt());
			tmp=tmp+"n="+Integer.toString(i)+" eta="+spx+" phi="+spy+" et="+spz+"\n";
		}

		return tmp;
	}

	/**
	 * Calculate delta R distance.
	 * 
	 * @param a
	 *            input particle
	 * @param b
	 *            input particle
	 * @param p
	 *            power parameter
	 * @return Kt distance
	 */
          public float getKtDistance12(ParticleF a, ParticleF b) {
                double rsq, esq, deltaEta, deltaPhi;
                deltaEta = a.getEta() - b.getEta();
                double phi1 = a.getPhi();
                double phi2 = b.getPhi();
                deltaPhi = phi2 - phi1;
                if (deltaPhi>Math.PI) deltaPhi=PI2-deltaPhi;
                if (deltaPhi<-Math.PI) deltaPhi=PI2+deltaPhi;
                rsq = (deltaEta*deltaEta + deltaPhi*deltaPhi);
                esq = 0;
                if      (mode==1) esq=Math.min(a.getEt2(), b.getEt2());         // kT
                else if (mode==0) esq=Math.min(a.getEt(), b.getEt());           // C-A
                else if (mode==-1) esq=Math.min(1.0/a.getEt2(), 1.0/b.getEt2()); // antiKT
                else esq=Math.min(a.getEt2(), b.getEt2());         // kT

                return (float)(esq * rsq/R2);
         }


          /**
         * Calculate  R distance in eta-phi.
         * 
         * @param a
         *            input particle
         * @param b
         *            input particle
         * @return eta-phi distance
         */
        public float getDistance(ParticleF a, ParticleF b) {
                double rsq, deltaEta, deltaPhi;
                deltaEta = a.getEta() - b.getEta();
                double phi1 = a.getPhi();
                double phi2 = b.getPhi();
                deltaPhi = phi2 - phi1;
                if (deltaPhi > Math.PI)
                        deltaPhi = PI2 - deltaPhi;
                if (deltaPhi < -Math.PI)
                        deltaPhi = PI2 + deltaPhi;
                rsq = (deltaEta * deltaEta + deltaPhi * deltaPhi);
                return (float)Math.sqrt(rsq);
        }


	/**
	* This is the KT distance to the beam (assuming Z=Y=0). 
	* The distance measure depends on the mode parameter. 
	* @param a
	*            particle
	* @return kT distance
	*/
	public float getKtDistance1(ParticleF a) {
		if (mode==1) return (float)(a.getEt2());
		else if (mode==0) return (float)(a.getEt());
		else if (mode==-1) return (float)(1.0/a.getEt2());
		return (float)(a.getEt2());
	}



	/**
	 * Print debugging information. It shows how much time spend to make jets in ms.
	 * @param debug true if printing benchmark information. 
	 */
	public void setDebug(boolean debug){

		this.debug=debug;
	}

	/**
	 * Main class for testing.
	 * 
	 * @param args
	 */

	public static void main(String[] args) {



                // for correct benchmark with C++ (after just-in-time compiler)
                for (int i=0; i<3; i++){

		ArrayList<ParticleF> list = new ArrayList<ParticleF>();
		try {
			File file = new File("jets/single-event.dat");
			System.out
			.println("Reading test file with jets. Number of particles:");
			FileReader fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			String line;
			while ((line = bufferedReader.readLine()) != null) {

				StringTokenizer st = new StringTokenizer(line);
				int j = 0;
				double[] mom = new double[4];
				while (st.hasMoreElements()) {
					Double d = Double.parseDouble(st.nextElement().toString());
					mom[j] = d;
					j++;
				}

				// px,py,pz,e
				ParticleF  pp = new  ParticleF();
				pp.setPxPyPzE((float)mom[0], (float)mom[1], (float)mom[2], (float)mom[3]);
				list.add(pp);

			}
			fileReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// System.out.println(list.size());
		KTjet ktjet = new KTjet(0.6, 1, -1, 5.0);
		ktjet.setDebug(true);
		ktjet.buildJets(list);
		ktjet.printJets();

                } // end multiple runs



	}

}
