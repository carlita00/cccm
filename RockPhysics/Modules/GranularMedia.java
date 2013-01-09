package Modules;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.dsp.*;


/**
 * Effective elastic media:  bounds and mixing laws
 * Reference Rock physics handbook
 * @author Carla Carvajal, Colorado School of Mines.
 * @version 2012.12.21
 */


public class GranularMedia {
   /**
   * Elastic Media constructor.
   */
  public GranularMedia(float n , float phiC,float gs,float ks, float prs, float[] p) {
    _n    = n;
    _phiC = phiC;
    _gs   = gs;
    _ks   = ks;
    _prs  = prs;
    _p    = p;
  }
  /**
   * Hertz mindlin
   * Serial version
   * Page 177 
   * PARAMETERS
   * @ n     = coordination number
   * @ phiC  = critical porosity
   * @ gs    = shear modulus of the min
   * @ p     = effective pressure vector
   * @ prs   = poisson moduli
   */
  public float[] applyKhm() {
    int n1 = _p.length;
    float [] k = new float[n1];
    for (int i=0; i<n1; ++i){
      k[i] = (_n*_n*pow(1-_phiC,2)*_gs*_gs);
      k[i] /= (18*pow(PI,2)*pow(1-_prs,2));  
      k[i] *= _p[i]; 
      k[i] = pow(k[i],(1.0f/3.0f));
    }
    return k;
  }

  public float[] applyGhm() {
    int n1 = _p.length;
    float [] g = new float[n1];
    for (int i=0; i<n1; ++i){
      g[i] = (3.0f*_n*_n*pow(1f-_phiC,2)*_gs*_gs);
      g[i]/= (2.0f*pow(PI,2.0)*pow(1f-_prs,2));
      g[i]*= _p[i]; 
      g[i] = ((5.0f-4.0f*_prs)/(10f-5f*_prs))*pow(g[i],0.333f);
    }
    return g;
  }

  public float[] applyZhm(){
    int n1 = _p.length;
    float [] z = new float[n1];
    float[] ghm = applyGhm();
    float[] khm = applyKhm();
    for (int i=0; i<n1; ++i){
      z[i] = ghm[i]/6.0f*(9.0f*khm[i]+8*ghm[i]);
      z[i]/= khm[i]+2*ghm[i];  
    }
    return z;
  }

  public float applyZ(){
    float z = _gs/6.0f*(9f*_ks+8f*_gs);
    z/=_ks+2*_gs;  
    return z;

  }


  ////////////////////////////Private///////////////////////////////////

  private float _n; 
  private float _phiC;    
  private float _gs;
  private float _ks;  
  private float _prs;       
  private float[] _p; 



}

