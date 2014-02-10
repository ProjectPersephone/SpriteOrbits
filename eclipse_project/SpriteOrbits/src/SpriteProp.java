/**
 * SpriteProp - for KickSat Sprites
 */

import java.util.List;
import java.util.ArrayList;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.util.MathUtils;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;
import org.orekit.errors.OrekitException;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Transform;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.apache.commons.math3.util.FastMath;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.time.AbsoluteDate;

/*
 * @author Luc Maisonobe
 * @author Fabien Maussion
 * @author Pascal Parraud
 * @author Romain Di Costanzo
 * @author Michael Turner
 *
 */
public class SpriteProp {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
         int n;
         AbsoluteDate releaseDate = AbsoluteDate.J2000_EPOCH;
         double dv;
         NumericalPropagator kickSatPropagator;
         KeplerianOrbit orbit;

         // Keplerian parameters taken from http://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html
         try {
			orbit = new org.orekit.orbits.KeplerianOrbit(
					 6794126.00, //  a, semi-major axis -  size.
					 .0014140, // meters: e, eccentricity -  shape.
					 Math.toRadians(51.83436), // i, inclination -  orientation w.r.t equator.
					 Math.toRadians(82.92739), // pa, perigee
					 Math.toRadians(31.76736), // degrees: raan, Right Ascension of the Ascending Node
					 Math.toRadians(46.64131), // anomaly, TA/MA True/Mean Anomaly
					 PositionAngle.valueOf("TRUE"), // type,
					 FramesFactory.getITRF (IERSConventions.IERS_2010, true), // frame,
					 AbsoluteDate.J2000_EPOCH, // date
					 Constants.EIGEN5C_EARTH_MU // central attraction coefficient MU
					 );

         n = 200; // egregious MAGIC - MT
         
// copied from https://www.orekit.org/static/architecture/propagation.html 
      // steps limits
         final double minStep  = 0.001;
         final double maxStep  = 1000;
         final double initStep = 60;

         // error control parameters (absolute and relative)
         final double positionError = 10.0;
         final double[][] tolerances = NumericalPropagator.tolerances(positionError, orbit, orbit.getType());

         // set up mathematical integrator
         AdaptiveStepsizeIntegrator integrator =
             new DormandPrince853Integrator(minStep, maxStep, tolerances[0], tolerances[1]);
         integrator.setInitialStepSize(initStep);

         // set up space dynamics propagator
         kickSatPropagator = new NumericalPropagator(integrator);
// End of copy 
         
// Luc Maisonobe, from e-mail
//		- start with the initial orbit of KickSat
//		 - at the prescribed time of sprites release, compute position and
//		   attitude of KickSat using any orbit propagator, and a Sun pointing
//		   attitude (i.e. using CelestialBodyPointing for the attitude provider),
//		 - supposing KickSat will release n sprites uniformly spread with some
//		   known radial velocity in its plane orthogonal to Sun pointing, you can
//		   compute the sprites initial positions and velocities as follows:

		    SpacecraftState kickSat = kickSatPropagator.propagate(releaseDate);
		    Transform satToInertial = kickSat.toTransform().getInverse();
		    List<PVCoordinates> spritesPV = new ArrayList<PVCoordinates>(n);
		    
		    for (int i = 0; i < n; i++) {

		      // each sprite has a relative position of 0, 0, 0 wrt KicSat
		      // and a relative velocity at some angle in the release plane
		      // orthogonal to KickSat Z axis
		    	// [Turner]: guessing that Luc's "dv" is that "relative velocity".
		    	dv = 1.0; 
		      double alpha = i * MathUtils.TWO_PI / n;
		      PVCoordinates relativePV =
		          new PVCoordinates(Vector3D.ZERO,
		                            new Vector3D(dv * FastMath.cos(alpha),
		                                         dv * FastMath.sin(alpha),
		                                         0.0));

		      // convert the relative PV (i.e. position velocity expressed
		      // in spacecraft frame) into an absolute PV
		      spritesPV.add(satToInertial.transformPVCoordinates(relativePV));
		    // End of Luc Maisonobe e-mail contribution
		    }
 		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (OrekitException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
         
	}

}
