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


public class ElasticMedia {
   /**
   * Elastic Media constructor.
   */
  public ElasticMedia() {
  }

  /**
   * Voigt-Reuss-Hill average moduli estimate
   * Serial version
   * Page 177 
   * PARAMETERS
   * @ k = vector of size n containing the values of the moduli
           of the minerals
   * @ f = vector of size n containing the values of the composition
   *       of the minerals, have to add 1.
   */
  public float hillAverage(float[] k, float[] f) {
    int n = k.length;
    float mv = 0;
    float mr = 0;
    for (int i=0; i<n; ++i){
      mv += k[i]*f[i];
      mr += f[i]/k[i];
    }
    float mvr = 0.5f*(mv+1.0f/mr);
    return mvr;
  }
}

