/**
 * SpriteProp - for KickSat Sprites
 */

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
// import java.util.Date;
import java.util.List;
import java.util.Locale;

import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.attitudes.CelestialBodyPointed;
import org.orekit.attitudes.LofOffset;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.ZipJarCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.drag.IsotropicDrag; // import org.orekit.forces.SphericalSpacecraft;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.atmosphere.HarrisPriester;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.LOFType;
import org.orekit.frames.Transform;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateComponents;
import org.orekit.time.TimeComponents;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;

/*
 * @author Luc Maisonobe
 * @author Fabien Maussion
 * @author Pascal Parraud
 * @author Romain Di Costanzo
 * @author Michael Turner
 *
 */
public class SpriteProp {

    /** Orbit propagator for KickSat itself. */
    private final Propagator kickSatPropagator;

    /** Sprites propagators. */
    private final List<Propagator> spritesPropagators;

    /** Earth model. */
    private OneAxisEllipsoid earth;

    /** Sun model. */
    private CelestialBody sun;

    /**
	 * @param args
	 */
	public static void main(String[] args) {

        PrintStream out = null;

        try {

            // set up orekit data
            // for this to work, the orekit-data.zip file avaiable at the following URL
            // must be put in the current working directory
            // URL:     https://www.orekit.org/forge/projects/orekit/files
            File userDir = new File(System.getProperty("user.dir"));
            File orekitZip = new File(userDir, "orekit-data.zip");
            DataProvidersManager.getInstance().addProvider(new ZipJarCrawler(orekitZip));

            // reference models
            final Frame eme2000   = FramesFactory.getEME2000(); // this is the frame labeled J2K in the NASA page
            final Frame itrf      = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
            final TimeScale utc   = TimeScalesFactory.getUTC();
            final double mu       = Constants.EIGEN5C_EARTH_MU; // central attraction coefficient MU

            // set up some data, ideally, this should be provided as input to the program
            // ISS orbit is from http://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html
            //
            //   Coasting Arc #15 (Orbit 3386)
            //   ---------------------------------------
            //   
            //   Vector Time (GMT): 2014/055/16:27:06.921
            //   Vector Time (MET): N/A
            //   Weight (LBS)     : 911651.1
            //
            //   ...
            //
            //                             J2K Cartesian  
            //                       --------------------------------
            //                       X    =         2998767.75
            //                       Y    =        -6097451.56  meter
            //      ...              Z    =         -141448.92
            //                       XDOT =        4323.077242
            //                       YDOT =        1994.291706  meter/sec
            //                       ZDOT =        6000.774574
            final int numberOfSprites = 12; // we start small, we can increase this later
            final double relativeReleaseVelocity = 1.0; // (m/s)
            final AbsoluteDate releaseDate = new AbsoluteDate(2014,3, 1,  // year, month, day
                                                              12, 0, 0.0, // hours, minutes, seconds
                                                              utc);

            final CartesianOrbit kickSatOrbit =
                    new CartesianOrbit(new PVCoordinates(new Vector3D (2998767.75, -6097451.56, -141448.92),    // position (m)
                    									 new Vector3D (4323.077242, 1994.291706, 6000.774574)), // velocity (m/s)
                                       eme2000,  
                                       new AbsoluteDate(new DateComponents(2014, 55),      // year, day in year as NASA page above
                                                        new TimeComponents(16, 27, 6.921), // hour in day
                                                        utc),
                                       mu);
            final double kickSatMass         = 10.0;   // kg
            final double kickSatCrossSection = 0.03;   // m^2
            final double kickSatDragCoeff    = 2.2;    // no units
            final double spriteMass          = 0.01;   // kg
            final double spriteCrossSection  = 2.5e-3; // m^2
            final double spriteDragCoeff     = 2.2;    // no units

            SpriteProp spriteProp = new SpriteProp(numberOfSprites, kickSatOrbit,
                                                   kickSatMass, kickSatCrossSection, kickSatDragCoeff,
                                                   spriteMass, spriteCrossSection, spriteDragCoeff,
                                                   relativeReleaseVelocity, releaseDate,
                                                   itrf);

            final double propagationDuration = 0.2;  // days after release
            final double step                = 60.0; // seconds
            //out = new PrintStream(new File(userDir, "sprites-prop.txt"));
            out = new PrintStream(new File(userDir, "orbits.json"));
            out.format(Locale.US, "[");
            spriteProp.run(out, propagationDuration, step, utc);
            out.format(Locale.US, "]");

        } catch (IOException e) {
            e.printStackTrace();
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        } catch (OrekitException e) {
            e.printStackTrace();
        } finally {
            if (out != null) {
                out.close();
            }
        }
         
	}

	/** Simple constructor.
	 * @param n number of sprites
	 * @param kickSatOrbit Kicksat orbit
     * @param kickSatMass mass of KickSat
     * @param kickSatCrossSection cross-section of KickSat
     * @param kickSatDragCoeff drag coefficient for KickSat
     * @param spriteMass mass of sprites
     * @param spriteCrossSection cross-section of sprites
     * @param spriteDragCoeff drag coefficient for sprites
	 * @param relativeReleaseVelocity release velocity (sprite wrt. KickSat) in m/s
	 * @param releaseDate release date
	 * @param itrf Earth frame
	 */
	private SpriteProp(final int n, final Orbit kickSatOrbit,
	                   final double kickSatMass, final double kickSatCrossSection, final double kickSatDragCoeff,
	                   final double spriteMass, final double spriteCrossSection, final double spriteDragCoeff,
                       final double relativeReleaseVelocity, final AbsoluteDate releaseDate,
                       final Frame itrf)
	    throws OrekitException {

        sun   = CelestialBodyFactory.getSun();
        earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                     Constants.WGS84_EARTH_FLATTENING,
                                     itrf);

        // Sun-pointing attitude
        final AttitudeProvider sunPointing =
                new CelestialBodyPointed(kickSatOrbit.getFrame(), sun,
                                         kickSatOrbit.getPVCoordinates().getMomentum(),
                                         Vector3D.PLUS_K, Vector3D.PLUS_J);
        kickSatPropagator = createPropagator(new SpacecraftState(kickSatOrbit, kickSatMass),
                                             sunPointing, kickSatCrossSection, kickSatDragCoeff);

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

		    spritesPropagators = new ArrayList<Propagator>(n);
		    for (int i = 0; i < n; i++) {

		      // each sprite has a relative position of 0, 0, 0 wrt KicSat
		      // and a relative velocity at some angle in the release plane
		      // orthogonal to KickSat Z axis
		      double alpha = i * MathUtils.TWO_PI / n;
		      PVCoordinates relativePV =
		          new PVCoordinates(Vector3D.ZERO,
		                            new Vector3D(relativeReleaseVelocity * FastMath.cos(alpha),
		                                         relativeReleaseVelocity * FastMath.sin(alpha),
		                                         0.0));

		      // convert the relative PV (i.e. position velocity expressed
		      // in spacecraft frame) into an absolute PV
              final Orbit spriteOrbit =
                      new CartesianOrbit(satToInertial.transformPVCoordinates(relativePV),
                                         kickSatOrbit.getFrame(), releaseDate, kickSatOrbit.getMu());
              final double roll  = 0.0; // this should be replaced by some random draw
              final double pitch = 0.0; // this should be replaced by some random draw
              final double yaw   = 0.0; // this should be replaced by some random draw
              spritesPropagators.add(createPropagator(new SpacecraftState(spriteOrbit, spriteMass),
                                                      new LofOffset(spriteOrbit.getFrame(), LOFType.VNC,
                                                                    RotationOrder.XYZ, roll, pitch, yaw),
                                                      spriteCrossSection, spriteDragCoeff));

// End of Luc Maisonobe e-mail contribution

		    }

	}

	/** Create a numerical propagator for a state.
	 * @param state state to propagate
	 * @param attitudeProvider provider for the attitude
	 * @param crossSection cross section of the object
	 * @param dragCoeff drag coefficient
	 */
	private Propagator createPropagator(final SpacecraftState state, final AttitudeProvider attitudeProvider,
	                                    final double crossSection, final double dragCoeff)
	  throws OrekitException {

	    // see https://www.orekit.org/static/architecture/propagation.html 
        // steps limits
         final double minStep  = 0.001;
         final double maxStep  = 1000;
         final double initStep = 60;

         // error control parameters (absolute and relative)
         final double positionError = 10.0;
         final OrbitType orbitType = OrbitType.CARTESIAN; // we will propagate in Cartesian coordinates
         final double[][] tolerances = NumericalPropagator.tolerances(positionError, state.getOrbit(), orbitType);

         // set up mathematical integrator
         AdaptiveStepsizeIntegrator integrator =
             new DormandPrince853Integrator(minStep, maxStep, tolerances[0], tolerances[1]);
         integrator.setInitialStepSize(initStep);

         // set up space dynamics propagator
         NumericalPropagator propagator = new NumericalPropagator(integrator);
         propagator.setOrbitType(orbitType);

         // add gravity field force model
         final NormalizedSphericalHarmonicsProvider gravityProvider =
                 GravityFieldFactory.getNormalizedProvider(8, 8);
         propagator.addForceModel(new HolmesFeatherstoneAttractionModel(earth.getBodyFrame(), gravityProvider));

         // add atmospheric drag force model
         propagator.addForceModel(new DragForce(new HarrisPriester(sun, earth),
                                                new IsotropicDrag(crossSection, dragCoeff, 0.0, 0.0)));

         // set attitude mode
         propagator.setAttitudeProvider(attitudeProvider);

         propagator.setInitialState(state);

         return propagator;

	}

	/** run the application.
	 * @param out output file
	 * @param propagationDuration duration of the propagation
	 * @param step fixed step between output lines
	 * @param utc UTC time scale
	 */
	public void run(final PrintStream out,
	                final double propagationDuration, final double step,
	                final TimeScale utc)
	    throws OrekitException {

	    final AbsoluteDate start = spritesPropagators.get(0).getInitialState().getDate();
	    final AbsoluteDate end   = start.shiftedBy(propagationDuration * Constants.JULIAN_DAY);
        /*out.format(Locale.US, "# file generated on %s%n",
                   new AbsoluteDate(new Date(), utc).toString(utc));
        out.format(Locale.US, "# propagating %d sprites from %s to %s%n",
                   spritesPropagators.size(), start.toString(utc), end.toString(utc));
        out.format(Locale.US, "# column 1:    date (UTC)%n");
        out.format(Locale.US, "# column 2:    date offset since start (seconds)%n");
        out.format(Locale.US, "# column 3:    KickSat geodetic latitude (degrees)%n");
        out.format(Locale.US, "# column 4:    KickSat geodetic longitude (degrees)%n");
        out.format(Locale.US, "# column 5:    KickSat geodetic altitude (meters)%n");
        out.format(Locale.US, "# column 3i+3: sprite i geodetic latitude (degrees)%n");
        out.format(Locale.US, "# column 3i+4: sprite i geodetic longitude (degrees)%n");
        out.format(Locale.US, "# column 3i+5: sprite i geodetic altitude (meters)%n");
        */

        // in order to speed up computation, we let the numerical propagator choose its
        // steps, and create ephemerides, then we will use the ephemerides with fixed
        // steps for output
        List<Propagator> ephemerides = new ArrayList<Propagator>(spritesPropagators.size());
        for (final Propagator spritePropagator : spritesPropagators) {
            spritePropagator.setEphemerisMode();
            spritePropagator.propagate(end);
              ephemerides.add(spritePropagator.getGeneratedEphemeris());
        }

        boolean firstTime = true;
        for (AbsoluteDate date = start; date.compareTo(end) < 0; date = date.shiftedBy(step)) {
            if(!firstTime){
            	out.format(Locale.US, ",");
            }
            else{
            	firstTime = false;
            }

            /*out.format(Locale.US, "%s %9.1f", date.toString(utc), date.durationFrom(start));

            final GeodeticPoint kickSatGP = geodeticPosition(kickSatPropagator, date);
            out.format(Locale.US, " %8.3f %8.3f %8.1f",
                       FastMath.toDegrees(kickSatGP.getLatitude()),
                       FastMath.toDegrees(kickSatGP.getLongitude()),
                       kickSatGP.getAltitude());

            for (final Propagator ephemeride : ephemerides) {
                final GeodeticPoint spriteGP = geodeticPosition(ephemeride, date);
                out.format(Locale.US, " %8.3f %8.3f %8.1f",
                           FastMath.toDegrees(spriteGP.getLatitude()),
                           FastMath.toDegrees(spriteGP.getLongitude()),
                           spriteGP.getAltitude());
            }
            */
        	out.format(Locale.US, "{\"date\":\"%s\",\"offset\":%.1f,", date.toString(utc), date.durationFrom(start));

            final GeodeticPoint kickSatGP = geodeticPosition(kickSatPropagator, date);
            out.format(Locale.US, "\"kicksat\":{\"lat\":%.3f,\"lng\":%.3f,\"alt\":%.1f},",
                       FastMath.toDegrees(kickSatGP.getLatitude()),
                       FastMath.toDegrees(kickSatGP.getLongitude()),
                       kickSatGP.getAltitude());

            out.format(Locale.US, "\"sprites\":[");
            boolean firstSprite = true;
            for (final Propagator ephemeride : ephemerides) {
                final GeodeticPoint spriteGP = geodeticPosition(ephemeride, date);
                if(!firstSprite){
                	out.format(Locale.US, ",");
                }
                else{
                	firstSprite = false;
                }
                out.format(Locale.US, "{\"lat\":%.3f,\"lng\":%.3f,\"alt\":%.1f}",
                           FastMath.toDegrees(spriteGP.getLatitude()),
                           FastMath.toDegrees(spriteGP.getLongitude()),
                           spriteGP.getAltitude());
            }
            out.format(Locale.US, "]}");
            out.println();

        }
	}

	/** Get geodetic position of one object.
	 * @param propagator propagator managing the object
	 * @param date date at which the position is requested
	 * @return geodetic position
	 */
	private GeodeticPoint geodeticPosition(final Propagator propagator, final AbsoluteDate date)
	    throws OrekitException {
	    final SpacecraftState state = propagator.propagate(date);
	    return earth.transform(state.getPVCoordinates().getPosition(), state.getFrame(), date);
	}

}
