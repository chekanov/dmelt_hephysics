package hephysics.particle;

import hephysics.vec.*;
import java.io.Serializable;
import java.util.Formatter;
import jhplot.gui.HelpBrowser;

/**
 * a HEP-type basic particle based on 4-Lorentz vector.
 * It is characterized by status, particle code, and coordinate position in (X,Y,Z,time).
 * 
 * @author Sergei Chekanov
 */

public class HEParticle extends LParticle implements Serializable {
	static final long serialVersionUID = 1L;
	protected int status;
	protected int m_pType;
	protected float isospin;
	protected float spin;
	protected float spaceparity;
	protected float chargeparity;
	protected int pdgcode;
	protected HEParticle parent;
	protected HepLorentzVector position;

	/**
	 * Define a particle with a name and mass
	 * 
	 * @param name
	 *            particle name
	 * @param mass
	 *            mass
	 */
	public HEParticle(String name, double mass) {
		super(name,mass);
		this.status = 0;
		this.m_pType = 0;
		this.chargeparity=0;
		this.isospin=0;
		this.spaceparity=0;
		this.pdgcode=0;
		this.spin=0;
	}

	/**
	 * Define a 3-momentum for a given particle. Energy is
	 * set to 0.
	 * 
	 * @param px
	 *            Px
	 * @param py
	 *            Py
	 * @param pz
	 *            Pz
	 */
	public HEParticle(double px, double py, double pz) {
		super(px, py, pz);
		this.status = 0;
		this.m_pType = 0;
		this.chargeparity=0;
		this.isospin=0;
		this.spaceparity=0;
		this.pdgcode=0;
		
	}

	/**
	 * Define a 4-momentum and energy particle
	 * 
	 * @param px
	 *            Px (or X position)
	 * @param py
	 *            Py (or Y position)
	 * @param pz
	 *            Pz (or Z position)
	 * @param energy
	 *            energy (or time)
	 * @param mass
	 *            mass
	 */
	public HEParticle(double px, double py, double pz, double energy, double mass) {
		super(px, py, pz, energy,mass);
		this.status = 0;
		this.m_pType = 0;
		this.chargeparity=0;
		this.isospin=0;
		this.spaceparity=0;
		this.pdgcode=0;
		this.spin=0;
		
	}


	/**
	 * Define a particle in momentum and position
	 * @param momentum 4-momentum
	 * @param position postion in X,Y,Z,time
	 * @param mass mass
	 */
	public HEParticle(HepLorentzVector momentum, HepLorentzVector position, double mass) {
		super(momentum.px(), momentum.py(), momentum.pz(), momentum.e(),mass);
		this.position = position;
		this.status = 0;
		this.m_pType = 0;
		this.chargeparity=0;
		this.isospin=0;
		this.spaceparity=0;
		this.pdgcode=0;
		this.spin=0;
	}

	
	/**
	 * Set position of particle in X,Y,Z,time.
	 * @param position
	 */
	
	public void setPosition(HepLorentzVector position) {
		this.position = position;
	}
	
	
	/**
	 * Get position in X,Y,Z,time
	 * @return
	 */
	public HepLorentzVector getPosition() {
		return position;
	}
	

	
	/**
	 * Define a particle in momentum space.  Mass is set to 0
	 * 
	 * @param px
	 *            Px or X position
	 * @param py
	 *            Py or Y position
	 * @param pz
	 *            Pz or Z position
	 * 
	 * @param energy
	 *            energy or time
	 */
	public HEParticle(double px, double py, double pz, double energy) {
		super(px, py, pz, energy);
		this.status = 0;
		this.m_pType = 0;
		this.chargeparity=0;
		this.isospin=0;
		this.spaceparity=0;
		this.pdgcode=0;
		this.spin=0;
		
	}

          /**
         * Compute rapidity. 0.5*log( (m+z)/(m-z) );
         * @return
         */
        public double rapidity()  {
                double rapidity=-10e10;
                if (e()>pz())  rapidity=0.5*Math.log( (e()+pz())/(e()-pz()) );
                return rapidity; 
        }


	/**
	 * Define a Lorentz particle in momentum space.
	 * 
	 * @param name
	 *            Name of particle
	 * @param px
	 *            px (or X)
	 * @param py
	 *            py (or Y)
	 * @param pz
	 *            pz (or Z)
	 * @param energy
	 *            energy or time
	 * @param mass
	 *            mass
	 */
	public HEParticle(String name, double px, double py, double pz,
			double energy, double mass) {
		super(px, py, pz, energy,mass);
		this.status = 0;
		this.m_pType = 0;
		this.chargeparity=0;
		this.isospin=0;
		this.spaceparity=0;
		this.pdgcode=0;
		this.spin=0;
	}

	
	
	
	

	/**
	 * Add 2 particles
	 * 
	 * @param another
	 *            particle to be added
	 */
	public void add(HEParticle another) {

		add(another);
		this.chargeparity +=another.getChargeParity();
		this.spaceparity +=another.getSpaceparity();
		this.isospin +=another.getIsospin();
		this.spin +=another.getSpin();
		
	}

	
	

	
	/**
	 * Set isospin 
	 * @param isospin
	 */
	public void setIsospin(float isospin) {
		this.isospin = isospin;
	}
	
	
	/**
	 * Get isospin
	 * @return
	 */
	public float getIsospin() {
		return isospin;
	}
	

	/**
	 * Set C charge conjugation parity
	 * @param chargeparity charge parity
	 */
	public void setChargeParity(float chargeparity) {
		this.chargeparity = chargeparity;
	}
	
	/**
	 * Returns C charge conjugation parity
	 * @return
	 */
	public float getChargeParity() {
		return chargeparity;
	}
	
	
	/**
	 * Set J total spin
	 * @param spin
	 */
	public void setSpin(float spin) {
		this.spin = spin;
	}
	
	/**
	 * Get J total spin
	 * @return
	 */
	public float getSpin() {
		return spin;
	}
	
	
	/**
	 * Set PDG code
	 * @param pdgcode
	 */
	public void setPdgcode(int pdgcode) {
		this.pdgcode = pdgcode;
	}
	
	
	/**
	 * Get PDG code
	 * @return
	 */
	public int getPdgcode() {
		return pdgcode;
	}
	
	
	/**
	 * Set P space parity
	 * @param spaceparity
	 */
	public void setSpaceParity(float spaceparity) {
		this.spaceparity = spaceparity;
	}
	
	
	/**
	 * get P space parity
	 * @return
	 */
	public float getSpaceparity() {
		return spaceparity;
	}
	
	/**
	 * Status Monte Carlo code 
	 * 
	 * @return status
	 */
	public int getStatus() {
		return status;
	}

	/**
	 * Set Monte Carlo status code
	 * 
	 * @param status
	 */
	public void setStatus(int status) {
		this.status = status;
	}

	/**
	 * Get barcode 
	 * 
	 * @return barcode 
	 */
	public int getBarcode() {
		return m_pType;
	}

	/**
	 * Set barcode ID
	 * 
	 * @param barcode 
	 */
	public void setBarcode(int barcode) {
		this.m_pType = barcode;
	}

	
	/**
	 * Convert to string
	 */
	public String toString() {

		Formatter formatter = new Formatter();
		formatter.format("M=%9.4g,E=%9.4g,T=%9.4g,S=%9.4g,Type=%9.4,Isospin=%6.3,ChargeParity=%6.3,SpaceParity=%6.3,Spin=%6.3", mass, energy, t, status, m_pType,isospin,chargeparity,spaceparity,spin);
		String sname = String.format("%10s", name);

		return new String(sname + " " + v.toString() + ", "
				+ formatter.out().toString());

	}

	/**
	 * Make an exact copy of this particle
	 * 
	 * @return new copy
	 */
	public HEParticle copy() {
		HEParticle tmp = new HEParticle(getName(), pz(), pz(), pz(), e(), getMass());
		tmp.setParent(getParent());
		setT(t());
		setCharge(getCharge());
		setBarcode(getBarcode());
		setStatus(getStatus());
		setIsospin(getIsospin());
		setSpaceParity(getSpaceparity());
		setChargeParity(getChargeParity());
		setSpin(getSpin());
		setPdgcode(getPdgcode());
		return tmp;
	}

	
	/**
	 * Returns a light-weight particle class
	 * @return
	 */
	
	public LParticle getLParticle(){
		return this;
		
	}
	
	/**
	 * Evaluates 4-vector of decay product in the rest frame of parent.
	 * 
	 * @param prod1
	 *            first decay product
	 * @param prod2
	 *            second decay product
	 * @param randomRotate
	 *            is Phi randomly rotated?
	 * @return
	 */
	public void twoBodyDecay(HEParticle prod1, HEParticle prod2,
			boolean randomRotate) {

		twoBodyDecay( prod1.getLParticle(),prod2.getLParticle(),randomRotate );
	}
	
	/**
	 * Lorentz Boost
	 * 
	 * @param parent
	 *            parent particle
	 */
	public void boost(HEParticle  parent) {
		 boost(parent.getLParticle());
	}
	
	/**
	 * Print particle
	 */

	public void print() {
		System.out.println(toString());
	}

           /**
         * Show online documentation.
         */
        public void doc() {

                String a = this.getClass().getName();
                a = a.replace(".", "/") + ".html";
                new HelpBrowser(HelpBrowser.JHPLOT_HTTP + a);

        }

}
