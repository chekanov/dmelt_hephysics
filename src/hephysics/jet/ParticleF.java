package hephysics.jet;
import java.io.Serializable;
import java.text.*;
import jhplot.gui.HelpBrowser;
import hephysics.particle.LParticle;

/**
 * A class representing a jet or particle with pre-computed et, et2, phi, eta (float values). It uses floats to speedup calculations for jet fundings. The class is used by light-weight  {@link hephysics.jets.KTjet} algorithm.  
 * "F" means float calculations The class has a minimum dynamic computation to minimize CPU. Use
 * {@link hephysics.particle.LParticle}  for (slower) dynamic calculations. 
 * To use double precision calculations, use {@link hephysics.jets.ParticleD} class. 
 * 
 * @author sergei
 * 
 */
public class ParticleF implements Comparable<ParticleF>, Serializable  {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private  float et2;

	private float eta;

	private float phi;

	private float et;

	private DecimalFormat formatter = new DecimalFormat("0.###E0");


	private final double PI2 = Math.PI * 2;

	/**
	 * Initialize pseudoparticle.
	 * 
	 * @param p
	 *            initialisation particle. 
	 */
	public ParticleF(LParticle p) {

		et2 = (float)p.et2();
		eta = (float)p.pseudoRapidity();
		phi = (float)p.phi();
		et =  (float)p.et();
	}

	/**
	 * Initialize pseudoparticle.
	 * 
	 */
	public ParticleF() {
		et2 =0;
		eta =0;
		phi =0;
		et = 0;
	}

	/**
	 * Initialize fast particle from 4-momenta. The methods precomputes
	        * internal variables et, et, phi, eta. 
	 * @param px
	 * @param py
	 * @param pz
	 * @param energy
	 */
	public  ParticleF(float px, float py, float  pz, float energy) {
		//http://particle.physics.ucdavis.edu/seminars/data/media/2006/dec/soper.pdf
		LParticle p = new LParticle(px,py,pz,energy);
		et2 = (float)p.et2();
		eta = (float)p.pseudoRapidity();
		phi = (float)p.phi();
		et =  (float)p.et();
	}


	/**
	 * Initialize from 4-momenta. Recompute et2, et, eta, and phi.
	 * @param px
	 * @param py
	 * @param pz
	 * @param energy
	 */
	public void setPxPyPzE(float px, float py, float  pz, float energy) {
		LParticle p = new LParticle(px,py,pz,energy);
		et2 = (float)p.et2();
		eta = (float)p.pseudoRapidity();
		phi = (float)p.phi();
		et =  (float)p.et();

	}


	public float getEta() {
		return this.eta;
	}


	public float getPhi() {
		return this.phi;
	}

	public float getEt2() {
		return this.et2;
	}

	public float getEt() {
		return (float)et;
	}

       /**
        * Same as getEt(). Returns perp without calculations. 
        **/
        public float perp() {
                return (float)et;
        }

    
        /**
        * Same as getEta(). Returns pseudo-rapidity without calculations. 
        **/
        public float eta() {
                return (float)eta;
        }

	public void setEta(float eta) {
		this.eta = eta;
	}

	public void setEt2(float et2) {
		this.et2 = et2;
	}

	public void setPhi(float phi) {
		this.phi = phi;
	}


	public LParticle getLParticle() {

		double apt = Math.abs(this.et);
		double px = apt * Math.cos(this.phi);
		double py = apt * Math.sin(this.phi);
		double pz = apt * Math.sinh(this.eta);
		double theta=Math.atan2(Math.sqrt(et2),pz);
		double energy = this.et/Math.sin(theta);
		LParticle pp = new LParticle(px,py,pz,energy);
		return pp;

	}



         /**
         * Comparator. using perp2  for comparison (in increasing order)
         * 
         * @param o
         * @return
         */
        public int compareTo(ParticleF o) {
                if (et<o.getEt()) return 1;
                if (et>o.getEt()) return -1;
                return 0;
        }


	/**
	 * Show online documentation.
	 */
	public void doc() {

		String a = this.getClass().getName();
		a = a.replace(".", "/") + ".html";
		new HelpBrowser(HelpBrowser.JHPLOT_HTTP + a);

	}

	/**
	* Convert a particle to a string. 
	* @return a string with the particle 
	*/
	public String toString() {
		String se=formatter.format(et);
		String srap=formatter.format(eta);
		String sphi=formatter.format(phi);
		return "et="+se+" eta="+srap+" phi="+sphi;
	}


	/**
	 * Add to this particle and recalculate all characteristics.
	 * @param a
	 */
	public  void add(ParticleF a) {

		double et2_2=a.getEt2();
		double et_2= a.getEt();
		double eta_2= a.getEta();
		double phi_2= a.getPhi();

		double et2_1=this.et2;
		double et_1= this.et;
		double eta_1= this.eta;
		double phi_1= this.phi;

		this.et=(float)(et_1+et_2);
		this.et2=(float)(et2_1+et2_2);
		this.eta=(float)(((et_1*eta_1) + (et_2*eta_2))/this.et);
		this.phi=(float)(((et_1*phi_1) + (et_2*phi_2))/this.et);


	}

	/**
	    * Get a hash code
	    */
	public int hashCode() {
		return hashCode() + (int) Double.doubleToRawLongBits(et2);
	}

}
