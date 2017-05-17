package hephysics.particle;

import hephysics.vec.Hep3Vector;
import hephysics.vec.HepLorentzVector;
import java.io.Serializable;
import java.util.Formatter;
import java.util.Random;
import net.jafama.FastMath;
import jhplot.gui.HelpBrowser;

/**
 * A class representing a Lorentz particle. It is based on four-Lorentz vector
 * and extends HepLorentzVector. A particle can be represented by either
 * (px,py,pz,E) or (x,y,z,time). In addition, you can define name, mass,charge
 * and a parent particle as an object.
 * 
 * 
 * @author Sergei Chekanov
 */
public class LParticle extends HepLorentzVector implements
		Comparable<LParticle>, Serializable {

	static final long serialVersionUID = -6544699016896436061L;
	protected double mass;
	protected String name;
	protected LParticle parent;
	protected double charge;

	/**
	 * Define a dummy particle. Set Px,Py,Pz,mass to zero.
	 * 
	 */
	public LParticle() {
		super(0.0, 0.0, 0.0, 0.0);
		this.name = "";
		this.mass = 0;
		this.charge = 0;
		this.energy = 0;
	}

	/**
	 * Define a particle with a name. Set Px,Py,Pz,mass to zero.
	 * 
	 * @param name
	 *            particle name
	 */
	public LParticle(String name) {
		super(0.0, 0.0, 0.0, 0.0);
		this.name = name;
		this.mass = 0;
		this.charge = 0;
		this.energy = 0;
	}

	/**
	 * Define a particle with a name and mass. Set Px,Py,Pz to zero.
	 * 
	 * @param name
	 *            particle name
	 * @param mass
	 *            mass
	 */
	public LParticle(String name, double mass) {
		super(0.0, 0.0, 0.0, 0.0);
		this.name = name;
		this.mass = mass;
		this.charge = 0;
		this.energy = 0;
	}

	/**
	 * Define a 3-momentum or X,Y,Z position. Energy and mass are set to 0.
	 * 
	 * @param px
	 *            Px (or X_)
	 * @param py
	 *            Py (or Y)
	 * @param pz
	 *            Pz (or Z)
	 */
	public LParticle(double px, double py, double pz) {
		super(0.0, px, py, pz);
		this.mass = 0;
		this.charge = 0;
		this.energy = 0;
	}

	/**
	 * Set Px,Py,Pz (or x,y,z)
	 * 
	 * @param px
	 * @param py
	 * @param pz
	 */
	public void setPxPyPz(double px, double py, double pz) {
		setPx(px);
		setPy(py);
		setPz(pz);
	}

	/**
	 * Set Px,Py,Pz and energy of a particle
	 * 
	 * @param px
	 *            Px momentum
	 * @param py
	 *            Py momentum
	 * @param pz
	 *            Pz momentum
	 * @param e
	 *            energy
	 */
	public void setPxPyPzE(double px, double py, double pz, double e) {
		setPx(px);
		setPy(py);
		setPz(pz);
		setE(e);
	}

	/**
	 * Set x,y,z and time for a particle (vertex position)
	 * 
	 * @param x
	 *            X position
	 * @param y
	 *            Y position
	 * @param z
	 *            Z position
	 * @param t
	 *            time
	 */
	public void setXYZT(double x, double y, double z, double t) {
		setPx(x);
		setPy(y);
		setPz(z);
		setT(t);
	}

	/**
	 * Define a particle in momentum or coordinate space.
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
	public LParticle(double px, double py, double pz, double energy, double mass) {
		super(px, py, pz, energy);
		this.mass = mass;
		this.charge = 0;
	}

	/**
	 * Define a particle in momentum or 4-space time. Mass is set to 0
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
	public LParticle(double px, double py, double pz, double energy) {
		super(px, py, pz, energy);
		this.mass = 0;
		this.charge = 0;
	}

	/**
	 * Define a Lorentz particle in momentum or coordinate space.
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
	public LParticle(String name, double px, double py, double pz,
			double energy, double mass) {
		super(px, py, pz, energy);
		this.name = name;
		this.mass = mass;
		this.charge = 0;
	}

	/**
	 * Set a parent particle
	 * 
	 * @param parent
	 *            parent particle
	 */
	public void setParent(LParticle parent) {
		this.parent = parent;
	}

	/**
	 * Scale all components (Px,Py,Px,E) by a factor C. Keep mass unchanged.
	 * This means all components are multiplied by this factor.
	 * 
	 * @param scale
	 *            factor used to multiply Px,Py,Pz,E components.
	 * 
	 */
	public void scalePxPyPzE(double scale) {
		setPx(px() * scale);
		setPy(py() * scale);
		setPz(pz() * scale);
		setE(e() * scale);
	}

	/**
	 * Scale all components (X,Y,Z,T) by a factor This means all components are
	 * multiplied by this factor.
	 * 
	 * @param scale
	 *            factor used to multiply X,Y,Z,T components.
	 * 
	 */
	public void scaleXYZT(double scale) {
		setX(x() * scale);
		setY(y() * scale);
		setZ(z() * scale);
		setT(t() * scale);
	}

	/**
	 * Divide all components (Px,Py,Px,E) by a factor C. Useful for varints
	 * (integer) representation.
	 * 
	 * @param scale
	 *            factor used to divide Px,Py,Pz,E
	 * 
	 */
	public void dividePxPyPzE(double scale) {
		setPx(px() / scale);
		setPy(py() / scale);
		setPz(pz() / scale);
		setE(e() / scale);
	}

	/**
	 * Get a parent particle
	 * 
	 * @return parent
	 */
	public LParticle getParent() {
		return this.parent;
	}

	/**
	 * Add 2 particle
	 * 
	 * @param another
	 *            particle to be edded
	 */
	public void add(LParticle another) {

		setPx(v.x() + another.v3().x());
		setPy(v.y() + another.v3().y());
		setPz(v.z() + another.v3().z());
		this.energy = this.energy + another.getE();
		this.t = this.t + another.getT();
		this.name = this.name + "+" + another.getName();
		this.charge = this.charge + another.getCharge();
	}

	/**
	 * Get a mass
	 * 
	 * @return mass
	 */
	public double getMass() {
		return mass;
	}

	/**
	 * Get a mass
	 * 
	 * @return mass
	 */
	public double mass() {
		return mass;
	}

	/**
	 * Get a name
	 * 
	 * @return name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Get a 3-vector woth (Px,Py,Pz) or (X,Y,Z)
	 */
	public Hep3Vector getV3() {
		return this.getV3();
	}

	/**
	 * Set a mass
	 * 
	 * @param mass
	 *            Mass
	 */
	public void setMass(double mass) {
		this.mass = mass;
	}

	/**
	 * Get a hash code
	 */
	public int hashCode() {
		return hashCode() + (int) Double.doubleToRawLongBits(energy);
	}

	/**
	 * Set HepLorentzVector using theta angle, phi and total momentum P. e =
	 * sqrt(P*P + mass*mass); <br>
	 * pX = P*sin(theta)*cos(phi); <br>
	 * pY = P*sin(theta)*sin(phi); <br>
	 * pZ = P*cos(theta);
	 * <p>
	 * 
	 * @param theta
	 *            theta
	 * @param phi
	 *            phi
	 * @param P
	 *            momentum
	 * @return
	 */
	public void setThetaPhiP(double theta, double phi, double P) {
		double e = FastMath.sqrt(P * P + getMass() * getMass());
		double pX = P * FastMath.sin(theta) * FastMath.cos(phi);
		double pY = P * FastMath.sin(theta) * FastMath.sin(phi);
		double pZ = P * FastMath.cos(theta);
		setV3(pX, pY, pZ);
		setE(e);
	}

	/**
	 * Set 4-momentum of a particle using transverse momentum (pt),
	 * pseudorapidity (Eta) and azimuthal angle (phi) and energy.
	 * <p>
	 * 
	 * @param pt
	 *            transverse momentum of a particle.
	 * @param eta
	 *            eta (pseudorapidity) .
	 * @param phi
	 *            azimuthal angle
	 * @param e
	 *            energy
	 */
	public void setPtEtaPhiE(double pt, double eta, double phi, double e) {
		double apt = FastMath.abs(pt);
		double px = apt * FastMath.cos(phi);
		double py = apt * FastMath.sin(phi);
		double pz = apt * FastMath.sinh(eta);
		setPxPyPzE(px, py, pz, e);
	}

	/**
	 * Set 4-momentum of a particle using transverse momentum (pt),
	 * pseudorapidity (eta) and azimuthal angle (phi) and mass (m). In this
	 * case, the energy is FastMath.sqrt(px*px+py*py+pz*pz-m*m).
	 * <p>
	 * 
	 * @param pt
	 *            transverse momentum of a particle.
	 * @param eta
	 *            eta (pseudorapidity) .
	 * @param phi
	 *            azimuthal angle
	 * @param m
	 *            mass
	 */
	public void setPtEtaPhiM(double pt, double eta, double phi, double m) {
		double apt = FastMath.abs(pt);
		double px = apt * FastMath.cos(phi);
		double py = apt * FastMath.sin(phi);
		double pz = apt * FastMath.sinh(eta);
		this.mass = m;
		double ee = px * px + py * py + pz * pz - m * m;
		if (ee < 0)
			ee = 0;
		else
			ee = FastMath.sqrt(ee);
		setPxPyPzE(px, py, pz, ee);
	}

	/**
	 * Magnitude
	 * 
	 * @return
	 */
	public double abs() {
		return FastMath.sqrt(skp(this));
	}

	/**
	 * Set charge
	 * 
	 * @param charge
	 *            charge
	 */
	public void setCharge(double charge) {
		this.charge = charge;
	}

	/**
	 * Get charge
	 * 
	 * @return
	 */
	public double getCharge() {
		return charge;
	}

	/**
	 * Angle between 2 vectors in rad
	 * 
	 * @param momentum
	 *            parent particle
	 * @return angle in rad
	 */
	public double angle(LParticle momentum) {

		if (abs() <= 0.0D || momentum.abs() <= 0.0D)
			return 0.0D;
		else
			return FastMath.acos(skp(momentum) / abs() / momentum.abs());

	}

	/**
	 * Get calculated mass as sqrt(e*e-px**2-py**2-pz**2)
	 * 
	 * @return
	 */
	public double calcMass() {

		double s = energy * energy - skp(this);

		if (s > 0)
			return FastMath.sqrt(s);
		else
			return -1 * FastMath.sqrt(-1.0 * s);

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
	public void twoBodyDecay(LParticle prod1, LParticle prod2,
			boolean randomRotate) {

		double m1 = prod1.getMass();
		double m2 = prod2.getMass();
		double m = mass;

		// check if particle m can decay
		if (m < m1 + m2) {
			System.out
					.println("twoBodyDecay  parent mass is less than sum of products.");
			prod1.setV3(0., 0., 0.);
			prod2.setV3(0., 0., 0.);
			return;
		}

		// CM energies and momentum
		double e1 = (m * m + m1 * m1 - m2 * m2) / (2.0 * m);
		double e2 = (m * m - m1 * m1 + m2 * m2) / (2.0 * m);
		double P = FastMath.sqrt(e1 * e1 - m1 * m1);

		double ran1 = 0;
		double ran2 = 0;

		if (randomRotate) {
			Random r = new Random();
			ran1 = r.nextDouble();
			ran2 = r.nextDouble();
		}

		// Isotropic random angles
		double theta = FastMath.acos(2.0 * ran1 - 1.0);
		double phi = 2.0 * FastMath.PI * ran2;

		double pX = P * FastMath.sin(theta) * FastMath.cos(phi);
		double pY = P * FastMath.sin(theta) * FastMath.sin(phi);
		double pZ = P * FastMath.cos(theta);

		// set the 4-momenta
		prod1.setV3(pX, pY, pZ);
		prod1.setE(e1);
		prod2.setV3(-pX, -pY, -pZ);
		prod2.setE(e2);

	}

	/**
	 * Lorentz Boost
	 * 
	 * @param parent
	 *            parent particle
	 */
	public void boost(LParticle parent) {

		// beta and gamma values
		double betax = parent.v3().x() / parent.e();
		double betay = parent.v3().y() / parent.e();
		double betaz = parent.v3().z() / parent.e();
		double beta2 = betax * betax + betay * betay + betaz * betaz;
		double gamma = 1.0 / FastMath.sqrt(1.0 - beta2);
		double dot = betax * v.x() + betay * v.y() + betaz * v.z();
		double prod = gamma * (gamma * dot / (1.0 + gamma) + energy);

		double pX = v.x() + betax * prod;
		double pY = v.y() + betay * prod;
		double pZ = v.z() + betaz * prod;
		double e = gamma * (energy + dot);

		setV3(pX, pY, pZ);
		this.energy = e;
		this.parent = parent;
	}

	/**
	 * Convert to string
	 */
	public String toString() {

		Formatter formatter = new Formatter();
		formatter.format("M=%9.4g,E=%9.4g,T=%9.4g", mass, energy, t);
		String sname = String.format("%10s", name);

		return new String(sname + " " + v.toString() + ", "
				+ formatter.out().toString());

	}

	/**
	 * Make an exact copy of this particle
	 * 
	 * @return new copy
	 */
	public LParticle copy() {
		LParticle tmp = new LParticle(getName(), pz(), pz(), pz(), e(),
				getMass());
		tmp.setParent(getParent());
		setT(t());
		setCharge(getCharge());
		return tmp;
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

	/**
	 * Comparator. using perp2  for comparison (in increasing order)
	 * 
	 * @param o
	 * @return
	 */
	public int compareTo(LParticle o) {

		int tmp =0;
		if (perp2() < o.perp2())
			tmp = 1;
		if (perp2()> o.perp2())
			tmp = -1;
		return tmp;
	}

}
