package hephysics.hepsim;

import hephysics.jet.*;
import hephysics.particle.*;
import hephysics.vec.*;
import java.util.*;
import net.jafama.FastMath;
import proto.*;
import promc.io.*;

/**
 * Utils classes to fill particles from ProMC data format of the HepSim repository with Monte Carlo predictions.
 * 
 * @author S.Chekanov
 * 
 */
public class PromcUtil {


      /**
      *  For a given momentum unit, create a Lorentz particle object at the position j.
      *  It keeps 4-momentum, postion, PDG ID, status etc. 
      *  @param pa 
      *           particle record from a ProMC file 
      *  @param energy_unit momentum / energy unit used to create a 4-momentum from a varint (int64)  
      *   @param length_unit length unit used to create (x,y,z,t) from a varint (int64)
      *  @param j   position of the particle in the ProMC record
      *  @return a  HEP Lorentz particle of the record.
      * 
      **/
      static public HEParticle getParticle(ProMC.ProMCEvent.Particles pa, int energy_unit, int length_unit, int j){

        HEParticle lp = new HEParticle(       pa.getPx(j) / (double) energy_unit,
                                                pa.getPy(j) / (double) energy_unit,
                                                pa.getPz(j) / (double) energy_unit,
                                                pa.getEnergy(j) / (double) energy_unit);
       lp.setMass(pa.getMass(j) / (double) energy_unit); 
       lp.setPdgcode(pa.getPdgId(j));
       lp.setStatus(pa.getStatus(j));
       lp.setCharge(pa.getCharge(j));
       lp.setBarcode(pa.getBarcode(j));

       // set position 
       HepLorentzVector pos= new HepLorentzVector(
                   pa.getX(j) / (double) length_unit,
                   pa.getY(j) / (double) length_unit,
                   pa.getZ(j) / (double) length_unit,
                   pa.getT(j) / (double) length_unit );
       lp.setPosition(pos);

       return lp;

      }



     
      /**
      *  For a given space/time unit, create a particle with postions X, Y,Z,Time.
      *  
      *  @param pa 
      *           particle record from a ProMC file. 
      *  @param unit space/time lenght unit used to convert a varint to a value. 
      *  @param j   position of particle in the record
      *  @return  particle position 
      * 
      **/
      static public LParticle getPosition(ProMC.ProMCEvent.Particles pa, int unit, int j){

        LParticle lp = new LParticle();
        lp.setXYZT(pa.getX(j) / (double) unit,
                   pa.getY(j) / (double) unit,
                   pa.getZ(j) / (double) unit,
                   pa.getT(j) / (double) unit);
       return lp;

      }
 


        /**
         * Fill particles from ProMC and apply some cuts if available. 
         * Also you can apply some basic cuts. 
         * @param header 
         *           header from ProMC file 
         * @param pa 
         *           particle record from ProMC file 
         * @param status 
         *            status code for particles to be accepted. For final state, use status==1 
         * @param pTmin 
         *            minimum pT cut for filled particles. If pTmin less than, do not apply any pT cut 
         * @param etaMax 
         *            maximum Eta (pseudorapidity) cut on filled particles. If etaMax greater than 1000, no any Eta cut is applied. 
         *                   
         **/
	static public List<LParticle> getLParticleList(ProMCHeaderFile.ProMCHeader header, ProMC.ProMCEvent.Particles pa, int status, double pTmin, double etaMax ) {

                int unit = header.getMomentumUnit();
                int lunit = header.getLengthUnit();
		List<LParticle> par = new ArrayList<LParticle>();

                if (pTmin<=0 && etaMax>=1000) { 
                for (int j = 0; j < pa.getPxCount(); j++) {
                if ( pa.getStatus(j)==status ) {
                   LParticle lp = new LParticle(pa.getPx(j) / (double) unit,
                                                pa.getPy(j) / (double) unit,
                                                pa.getPz(j) / (double) unit,
                                                pa.getEnergy(j) / (double) unit); 
                   par.add(lp);
                }
                }
                return par;
                }

          
               // when pT cut is applied
               if (pTmin>0 || etaMax<1000) {              
                for (int j = 0; j < pa.getPxCount(); j++) {
                if ( pa.getStatus(j)==status ) {
                   LParticle lp = new LParticle(pa.getPx(j) / (double) unit,
                                                pa.getPy(j) / (double) unit,
                                                pa.getPz(j) / (double) unit,
                                                pa.getEnergy(j) / (double) unit);
                   if (lp.perp()>pTmin && FastMath.abs(lp.pseudoRapidity())<etaMax) par.add(lp);
                }
                }
                return par;
                }




                return par;
		}

         /**
         * Fill particles from ProMC and apply into array after some cuts if available. 
         * Also you can apply some basic cuts. 
         * @param header 
         *           header from ProMC file 
         * @param pa 
         *           particle record from ProMC file 
         * @param status 
         *            status code for particles to be accepted. For final state, use status==1 
         * @param pTmin 
         *            minimum pT cut for filled particles. If pTmin less than, do not apply any pT cut 
         * @param etaMax 
         *            maximum Eta (pseudorapidity) cut on filled particles. If etaMax greater than 1000, no any Eta cut is applied. 
         *                   
         **/
static public LParticle[] getLParticleArray(ProMCHeaderFile.ProMCHeader header, ProMC.ProMCEvent.Particles pa, int status, double pTmin, double etaMax ) {

                List<LParticle> list=getLParticleList(header, pa, status,  pTmin, etaMax);
                
                LParticle[] lp = new  LParticle[list.size()];

                for (int j = 0; j < list.size(); j++) {
                 
                 lp[j]=(LParticle)list.get(j);
                
                }

          


                return lp;
		}




/**
 * Fill particles from ProMC and apply some cuts if available.  Use float precisions for particle definitions.
 * Also you can apply some basic cuts. 
 * @param header 
 *           header from ProMC file 
 * @param pa 
 *           particle record from ProMC file 
 * @param status 
 *            status code for particles to be accepted. For final state, use status==1 
 * @param pTmin 
 *            minimum pT cut for filled particles. If pTmin less than, do not apply any pT cut 
 * @param etaMax 
 *            maximum Eta (pseudorapidity) cut on filled particles. If etaMax greater than 1000, no any Eta cut is applied. 
 *                   
 **/

static public List<ParticleF> getParticleFList(ProMCHeaderFile.ProMCHeader header, ProMC.ProMCEvent.Particles pa, int status, double pTmin, double etaMax ) {

        int unit = header.getMomentumUnit();
        int lunit = header.getLengthUnit();
        List<ParticleF> par = new ArrayList<ParticleF>();

        if (pTmin<=0 && etaMax>=1000) { 
        for (int j = 0; j < pa.getPxCount(); j++) {
        if ( pa.getStatus(j)==status ) {
        	ParticleF lp = new ParticleF(
        			                    (float)(pa.getPx(j) / (double) unit),
        			                    (float)(pa.getPy(j) / (double) unit),
        			                    (float)(pa.getPz(j) / (double) unit),
        			                    (float)(pa.getEnergy(j) / (double) unit)); 
           par.add(lp);
        }
        }
        return par;
        }

  
       // when pT cut is applied
       if (pTmin>0 || etaMax<1000) {              
        for (int j = 0; j < pa.getPxCount(); j++) {
        if ( pa.getStatus(j)==status ) {
        	ParticleF lp = new ParticleF(
                    (float)(pa.getPx(j) / (double) unit),
                    (float)(pa.getPy(j) / (double) unit),
                    (float)(pa.getPz(j) / (double) unit),
                    (float)(pa.getEnergy(j) / (double) unit)); 
        	
   
           if (lp.getEt()>pTmin && FastMath.abs(lp.getEta())<etaMax) par.add(lp);
        }
        }
        return par;
        }

        return par;
}



/**
 * Fill particles from ProMC from a typical NLO program, and apply some cuts if available.  Use float precisions for particle definitions.
 * Also you can apply some basic cuts. 
 * @param header 
 *           header from ProMC file 
 * @param pa 
 *           particle record from ProMC file 
 * @param pTmin 
 *            minimum pT cut for filled particles. If pTmin less than, do not apply any pT cut 
 * @param etaMax 
 *            maximum Eta (pseudorapidity) cut on filled particles. If etaMax greater than 1000, no any Eta cut is applied. 
 *                   
 **/
static public List<ParticleF> getParticleFList(pronlo.io.ProMCHeaderFile.ProMCHeader header, pronlo.io.ProMC.ProMCEvent.Particles pa, double pTmin, double etaMax ) {

        int unit = header.getMomentumUnit();
        int lunit = header.getLengthUnit();
        List<ParticleF> par = new ArrayList<ParticleF>();

        if (pTmin<=0 && etaMax>=1000) { 
        for (int j = 0; j < pa.getPxCount(); j++) {
        	ParticleF lp = new ParticleF(
        			                    (float)(pa.getPx(j) / (double) unit),
        			                    (float)(pa.getPy(j) / (double) unit),
        			                    (float)(pa.getPz(j) / (double) unit),
        			                    (float)(pa.getEnergy(j) / (double) unit)); 
           par.add(lp);
        }
        return par;
        }

  
       // when pT cut is applied
       if (pTmin>0 || etaMax<1000) {              
        for (int j = 0; j < pa.getPxCount(); j++) {
        	ParticleF lp = new ParticleF(
                    (float)(pa.getPx(j) / (double) unit),
                    (float)(pa.getPy(j) / (double) unit),
                    (float)(pa.getPz(j) / (double) unit),
                    (float)(pa.getEnergy(j) / (double) unit)); 
           if (lp.getEt()>pTmin && FastMath.abs(lp.getEta())<etaMax) par.add(lp);
        }
        return par;
        }

        return par;
}

/**
 * Fill particles from ProMC from a typical NLO program, and apply some cuts if available.  Use double precisions for particle definitions.
 * Also you can apply some basic cuts. 
 * @param header 
 *           header from ProMC file 
 * @param pa 
 *           particle record from ProMC file 
 * @param pTmin 
 *            minimum pT cut for filled particles. If pTmin less than, do not apply any pT cut 
 * @param etaMax 
 *            maximum Eta (pseudorapidity) cut on filled particles. If etaMax greater than 1000, no any Eta cut is applied. 
 *                   
 **/
static public List<ParticleD> getParticleDList(pronlo.io.ProMCHeaderFile.ProMCHeader header, pronlo.io.ProMC.ProMCEvent.Particles pa, double pTmin, double etaMax ) {

        int unit = header.getMomentumUnit();
        int lunit = header.getLengthUnit();
        List<ParticleD> par = new ArrayList<ParticleD>();

        if (pTmin<=0 && etaMax>=1000) { 
        for (int j = 0; j < pa.getPxCount(); j++) {
        	ParticleD lp = new ParticleD(
        			                    (pa.getPx(j) / (double) unit),
        			                    (pa.getPy(j) / (double) unit),
        			                    (pa.getPz(j) / (double) unit),
        			                    (pa.getEnergy(j) / (double) unit)); 
           par.add(lp);
        }
        return par;
        }

  
       // when pT cut is applied
       if (pTmin>0 || etaMax<1000) {              
        for (int j = 0; j < pa.getPxCount(); j++) {
        	ParticleD lp = new ParticleD(
                    (pa.getPx(j) / (double) unit),
                    (pa.getPy(j) / (double) unit),
                    (pa.getPz(j) / (double) unit),
                    (pa.getEnergy(j) / (double) unit)); 
           if (lp.getPt()>pTmin && FastMath.abs(lp.eta())<etaMax) par.add(lp);
        }
        return par;
        }

        return par;
}












/**
 * Fill particles from ProMC and apply some cuts if available. Use double precision particle definition.
 * Also you can apply some basic cuts. 
 * @param header 
 *           header from ProMC file 
 * @param pa 
 *           particle record from ProMC file 
 * @param status 
 *            status code for particles to be accepted. For final state, use status==1 
 * @param pTmin 
 *            minimum pT cut for filled particles. If pTmin less than, do not apply any pT cut 
 * @param etaMax 
 *            maximum Eta (pseudorapidity) cut on filled particles. If etaMax greater than 1000, no any Eta cut is applied. 
 *                   
 **/

static public List<ParticleD> getParticleDList(ProMCHeaderFile.ProMCHeader header, ProMC.ProMCEvent.Particles pa, int status, double pTmin, double etaMax ) {

        int unit = header.getMomentumUnit();
        int lunit = header.getLengthUnit();
        List<ParticleD> par = new ArrayList<ParticleD>();

        if (pTmin<=0 && etaMax>=1000) { 
        for (int j = 0; j < pa.getPxCount(); j++) {
        if ( pa.getStatus(j)==status ) {
        	ParticleD lp = new ParticleD(
        			                    (float)(pa.getPx(j) / (double) unit),
        			                    (float)(pa.getPy(j) / (double) unit),
        			                    (float)(pa.getPz(j) / (double) unit),
        			                    (float)(pa.getEnergy(j) / (double) unit)); 
           par.add(lp);
        }
        }
        return par;
        }

  
       // when pT cut is applied
       if (pTmin>0 || etaMax<1000) {              
        for (int j = 0; j < pa.getPxCount(); j++) {
        if ( pa.getStatus(j)==status ) {
        	ParticleD lp = new ParticleD(
                    (float)(pa.getPx(j) / (double) unit),
                    (float)(pa.getPy(j) / (double) unit),
                    (float)(pa.getPz(j) / (double) unit),
                    (float)(pa.getEnergy(j) / (double) unit)); 
        	
   
           if (lp.getPt()>pTmin && FastMath.abs(lp.eta())<etaMax) par.add(lp);
        }
        }
        return par;
        }

        return par;
}













}
