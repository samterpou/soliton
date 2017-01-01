import java.util.*;
/**
 * SOLITON
 * 
 * A miniscule program for attempts to observe emergent behavior from a system of chemical
 * species undergoing reaction and diffusion according to Gray-Scott mechanics.
 * 
 * A "soliton" is in fact "a quantum or quasiparticle propagated as a traveling nondissipative wave that
 * is neither preceded nor followed by another such disturbance."  I thought that these types of structures
 * may be the earliest and most primitive to arise in such a situation, and I also liked the name.
 * 
 * @author Sam Poulos
 * @version 1.0
 */
public class reactions
{
    // Dimensions of the simulation: nx = width, ny = height, nz = number of species/
    static private int nx = 200, ny = 200, nz = 50;
    // Self-explanatory.
    static private int nSteps = 1000000;
    // Diffusion is implemented in a very crude manner, with a von Neumann neighborhood of 1.
    static private float diffusionFactor = (float) 0.5;
    // Concentrations are clamped to 0.0 to 1.0, with float precision.
    static private float backgroundConcentration = (float) 0.01;
    // Program will write individual simulation frames to a "frames" folder in the current
    // directory.
    static private String outputDirectory = "./frames/";
    
    
    static private Random ran = new Random();
    static private int[] primes = { 2,3,5,7,11,13,17,19,23,29,31,37 };
    // The reaction order is the number of reactions for which an individual chemical species
    // serves as a reactant.
    static private int reactionOrder = 6;

    public static void main(String args[])
    {
        new reactions();
    }
    
    public reactions() {
        // Introductory message and parameters.
        System.out.println("ii: Beginning simulation with the following parameters:");
        System.out.println("ii:  width = "+nx+", height = "+ny+", #species = "+nz);
        System.out.println("ii:  diffusion factor = "+diffusionFactor);
        System.out.println("ii:  background concentration = "+backgroundConcentration);
        
        // Grid objects are essentially the chemical vats.  Create one and initialize it.
        Grid grid = new Grid(nx, ny, nz);        
        grid.setDiffusion(diffusionFactor);     
        grid.setBackgroundConcentration(backgroundConcentration);
        grid.setReactionOrder(reactionOrder);
        grid.setOutputDirectory(outputDirectory);
        
        // Create individual reactions for our chemical species.
        generateReactions(grid);

        grid.init();
        
        // Perform the simulation.
        for(int s=0;s<nSteps;s++) {
            grid.step();
            
            if(s%100==0) { System.out.println("Completed step "+s+"."); }
        }
        
        System.out.println("Simulation complete.");
    }
    
    // Used to tinker around.  Will return a random "rate" between 0.01 and 0.1.
    private float getRate() {
        return (float) Math.pow(10, ran.nextFloat()-2);
    }
    
    // Where most of my experimentation has been carried out.
    private void generateReactions(Grid grid) {
        int ip;
        int iz;
        int prime;
        int flip=1;
        
        // I wanted reactions to encompass all of the chemical species, but not set up "loops" where species would react
        // to form themselves again in a small number of reactions.  For instance, I didn't want s01 + s02 -> s11 + s12 -> s01 + s02.
        // 
        // At the moment, I am experimenting with s(i) + s(i+prime) -(s(i+2*prime))-> s(i+3*prime) + s(i+4*prime), where s(i+2*prime) is the enzyme.
        //
        // The reaction rate is equal to RATE * [reactant 1] * [reactant 2] * [product 1]^2 * [product 2]^2.
        // The volume of the reaction is limited to the amount of the limiting reagent, and to the lesser of (1-[product 1]) and
        // (1-[product 2]), as I like to keep the concentrations between 0.0 and 1.0 for simplicity.  Furthermore, the volume is limited
        // by another factor of reactionOrder, so that multiple reactions do not deplete all of a single reactant in one step.
        for(ip=0;ip<reactionOrder/2;ip++) {
            // As each reaction below names two reactants, iterating through the outer loop reactionOrder/2 times will still have each
            // species participate as a reactant in reactionOrder reactions.
            prime = primes[ip];
            for(iz=0;iz<nz;iz++) {
                
                // I have experimented with reaction constants differing wildly within orders of magnitude.  I usually
                // try to find a value between 1e-2 and 1e5 which allows the species to react, and hopefully produces
                // more than just oscillating spirals.
                grid.registerReaction(iz,iz+prime*flip,iz+prime*flip*3,iz+prime*flip*4,iz+prime*flip*2, (float) 1e5);
            }
            
            flip *= -1;
        }
    }

}
