import java.util.*;
import java.io.*;
import javax.imageio.*;
import java.awt.*;
import java.awt.image.BufferedImage;
/**
 * A "grid" is loosely a vat where chemical reactions, diffusions and ancillary
 * operations take place.
 * 
 * @author Sam Poulos
 * @version 1.0
 */
public class Grid
{
    private int nx, ny, nz;
    // swamp is the main three-dimensional array which holds the concentration of
    // each chemical species at each grid point.
    private float[][][] swamp;
    // reactions holds the specifications for the reactions as generated in the
    // reactions class.
    private ArrayList<reactionParameters> reactions = new ArrayList<reactionParameters>();
    
    private int reactionOrder;
    private float diffusionFactor;
    private float backgroundConcentration;
    private String outputDirectory;
    
    private int stepNum = 0;
    private Random rand = new Random();
    
    public Grid(int ix, int iy, int iz)
    {
        nx = ix; ny = iy; nz = iz;
        swamp = new float[nx][ny][nz];
    }
    
    public void setBackgroundConcentration(float ic) {
        backgroundConcentration = ic;
    }
    
    public void setReactionOrder(int ro) {
        reactionOrder = ro;
    }
    
    public void setDiffusion(float id) {
        diffusionFactor = id;
    }
    
    public void setOutputDirectory(String id) {
        outputDirectory = id;
    }
    
    public void init() {
        int ix; int iy; int iz; int ib;
        
        // Initialize the concentration of each species at each point to
        // backgroundConcentration.
        for(ix=0;ix<nx;ix++) {
            for(iy=0;iy<ny;iy++) {
                for(iz=0;iz<nz;iz++) {
                    swamp[ix][iy][iz] = backgroundConcentration;
                }
            }
        }
        
        // Seed the vat with random noise.  After some experimentation, I found that
        // blocks of random concentrations worked best.
        int rx; int ry;        
        for(ib=0;ib<100;ib++) {
            rx=rand.nextInt(nx);
            ry=rand.nextInt(ny);
            for(ix=rx-2;ix<rx+3;ix++) {
                for(iy=ry-2;iy<ry+3;iy++) {
                    for(iz=0;iz<nz;iz++) {
                        swamp[wrapX(ix)][wrapY(iy)][iz] = (float) rand.nextFloat();
                    }
                }
            }
        }
    }
    
    // To simplify the initial code, species numbers can be outside [0, nz-1].  If they are, then they are
    // "wrapped" around to the other end of the range.
    public  void registerReaction(int substrate, int s2, int product, int p2, int enzyme, float rate)
    {
        reactions.add(new reactionParameters(wrapZ(substrate), wrapZ(s2), wrapZ(product), wrapZ(p2), wrapZ(enzyme), rate));
    }
    
    // Seems simple enough...
    public void step() {
        reactAndDiffuse();
        drawFrame();
        stepNum++;
    }
    
    
    // Where the magic happens.  React AND diffuse.
    public void reactAndDiffuse()
    {
        int ix; int iy; int iz; int ir; float rate;
        reactionParameters curReact;
        int enzyme; int substrate; int product;
        int s2; int p2;
        float volume;
        
        int size = reactions.size();
        float[] curColumn = new float[nz];
        float[][] curSlice = new float[nx][ny];
        
        // Feed in species 0 by setting its concentration to 1.0 in a square in the middle of the
        // vat.  Also allow this area to be a sink by setting the concentration of every other
        // species to 0.0.  Concentration will be drawn in via diffusion.
        for(ix=nx/2-10;ix<=nx/2+10;ix++) {
            for(iy=ny/2-10;iy<ny/2+10;iy++) {
                swamp[ix][iy][0] = (float) 1.0;
                for(iz=1;iz<nz;iz++) {
                    swamp[ix][iy][iz] = (float) 0.0;
                }
            }
        }

        // Iterate over every position in the vat.
        for(ix=0;ix<nx;ix++) {
            for(iy=0;iy<ny;iy++) {
                
                // Iterate over every reaction
                for(ir=0;ir<size;ir++) {
                    // Get the current reaction from the ArrayList.
                    curReact = reactions.get(ir);
                    
                    // Populate reaction-specific variables.
                    substrate = curReact.substrate;
                    product = curReact.product;
                    enzyme = curReact.enzyme;
                    rate = curReact.rate;
                    s2 = curReact.s2;
                    p2 = curReact.p2;
             
                    // The heart of the reaction: specifying that the rate is proportional
                    // to the concentration of the reactant and proportional to the square of
                    // the concentration of the product.  'rate' is also a proportionality constant
                    // which is specified when setting up the reactions.
                    volume = swamp[ix][iy][product] * swamp[ix][iy][product];
                    volume *= swamp[ix][iy][enzyme] * rate;
                    volume *= swamp[ix][iy][substrate];
                    
                    
                    // This is one of my recent areas of exploration.  I generalized the reactions
                    // two have two reactants and two products.  If you prefer vanilla GS, you can
                    // always specify reactions with only one reactant and one product.
                    if(s2 != -1) {
                        volume *= swamp[ix][iy][p2] * swamp[ix][iy][p2];
                        volume *= swamp[ix][iy][s2];
                    }
                    
                    // All of the reactions at a single location are run in parallel, and so for
                    // this crude implementation, I wanted to make sure that the reactants weren't
                    // depleted from all the reactions.  Hence the "1/reactionOrder" factor.
                    volume = (float) Math.min(volume, swamp[ix][iy][substrate]/reactionOrder);
                    volume = (float) Math.min(volume, (1-swamp[ix][iy][product])/reactionOrder);
                    
                    if(s2 != -1) {
                        volume = (float) Math.min(volume, swamp[ix][iy][s2]/reactionOrder);
                        volume = (float) Math.min(volume, (1-swamp[ix][iy][p2])/reactionOrder);
                    }
                    
                    // Shouldn't make a difference, but just in case of problems with limited
                    // precision.
                    volume = (float) Math.max(0, volume);
                    
                    // Record the volume of the reaction in a temporary buffer.
                    curColumn[substrate] -= volume;
                    curColumn[product] += volume;
                    
                    if(s2 != -1) {
                        curColumn[s2] -= volume;
                        curColumn[p2] += volume;
                    }
                }
                
                for(iz=0;iz<nz;iz++) {
                    // After all of the reactions are run for a single location, update the
                    // concentrations of the species and clear the buffer.
                    swamp[ix][iy][iz] += curColumn[iz];
                    curColumn[iz] = 0;
                }
            }
        }
        
        
        // Perform the diffusion, iterating by chemical species.
        for(iz=0;iz<nz;iz++) {
            for(ix=0;ix<nx;ix++) {
                for(iy=0;iy<ny;iy++) {
                    // Crude implementation of diffusion.  Take off a certain fraction
                    // of the concentration, and then distribute it equally to the
                    // 4 neighbors.  Record the diffusion into a temporary buffer.
                    volume = swamp[ix][iy][iz] * diffusionFactor;
                    curSlice[ix][iy] -= volume;
                    volume /= 4;
                    
                    curSlice[wrapX(ix)][wrapY(iy+1)] += volume;
                    curSlice[wrapX(ix)][wrapY(iy-1)] += volume;
                    curSlice[wrapX(ix+1)][wrapY(iy)] += volume;
                    curSlice[wrapX(ix-1)][wrapY(iy)] += volume;
                }
            }
            
            for(ix=0;ix<nx;ix++) {
                for(iy=0;iy<ny;iy++) {
                    // Update concentrations and clear the buffer.
                    swamp[ix][iy][iz] += curSlice[ix][iy];
                    curSlice[ix][iy] = 0;
                }
            }
        }
    }
    

    private void drawFrame() {
        BufferedImage bi = new BufferedImage(nx, ny, BufferedImage.TYPE_3BYTE_BGR);
        
        int ix; int iy; int iz; int ic; float curVal;
        float[][][] buff = new float[nx][ny][3];
        float[] maxes = new float[3];
        float[] mins = new float[3];
        
        float[] curCol = new float[3];
        float[] ranges = new float[3];
        
        // We want to map the concentrations of each species into an aesthetically pleasing
        // display of red, green and blue.  So at each point, the sum of the concentrations of
        // every third species will be placed into a buffer for further processing.  For instance,
        // red = [s0] + [s3] + [s6] + ..., green = [s1] + [s4] + [s7] + ..., blue = [s2] + [s5] + [s8] + ...
        for(ix=0;ix<nx;ix++) {
            for(iy=0;iy<ny;iy++) {
                for(iz=0;iz<nz;iz++) {
                        buff[ix][iy][iz%3] += swamp[ix][iy][iz];
                }
            }
        }
        
 
        // Now we want to find the global minimum and maximum for each set of concentrations, over the
        // entire vat.
        for(ic=0;ic<3;ic++) {
            mins[ic] = buff[0][0][ic];
        }
        
        
        for(ix=0;ix<nx;ix++) {
            for(iy=0;iy<ny;iy++) {
                // Skip an area containing the location of reactant source and sink.
                //
                // Note that if you want to modify the program so this area no longer exists, you
                // can delete the line below and another line further down.
                if(iy>=ny/2-10&&iy<=ny/2+10&&ix>=nx/2-10&&ix<=nx/2+10) { continue; }
                
                // Crude code for determining maxima and minima.
                for(ic=0;ic<3;ic++) {
                    curVal = buff[ix][iy][ic];
                    if(curVal < mins[ic]) {
                        mins[ic] = curVal;
                    }
                    if(curVal > maxes[ic]) {
                        maxes[ic] = curVal;
                    }
                }
            }
        }
        
        // Calculate the dynamic ranges.
        for(ic=0;ic<3;ic++) {
            ranges[ic] = maxes[ic] - mins[ic];
        }
        
        
        // Map the sums of concentrations to actual RGB values.
        for(ix=0;ix<nx;ix++) {
            for(iy=0;iy<ny;iy++) {
                for(ic=0;ic<3;ic++) {
                    // Again, the area of sources and sinks may skew the aesthetics, so just render
                    // it as black.
                    if((ix>=nx/2-10&&ix<=nx/2+10)&&(iy>=ny/2-10&&iy<=ny/2+10)) { curCol[ic] = 0; continue; }
                    
                    // Originally I wanted to fully utilize the dynamic range of each color intensity.  Then
                    // I found that raising things to the third power looked cooler sometimes...
                    curCol[ic] = (float) Math.pow((buff[ix][iy][ic] - mins[ic])/ranges[ic],3);
                }
                bi.setRGB(ix,iy,toRGB(curCol[0], curCol[1], curCol[2]));
            }
        }
        
        try {
            String fileName = outputDirectory+String.format("step%1$05d.png",stepNum);
            File file = new File(fileName);
            ImageIO.write(bi, "png", file);
        } catch (IOException e) {
        }
    }
    
    // From example code I found online.  Basically to fit the data from floats into an int for the
    // setRGB function of BufferedImage.
    private int toRGB(float ir, float ig, float ib) {
        int r = toInt(ir); int g = toInt(ig); int b = toInt(ib);
        int ret = (r << 16) | (g << 8) | b;
        return ret;
    }
    
    private int toInt(float f) {
        int ret = (int) Math.floor(f * 256);
        if(ret == 256) { ret = 255; }
        
        return ret;
    }
    
    private int wrapX(int x) {
        if(x < 0) {
            x += nx;
        } else if(x >= nx) {
            x -= nx;
        }
    
    return x;
    }
    
    
    private int wrapY(int y) {
        if(y < 0) {
            y += ny;
        } else if(y >= ny) {
            y -= ny;
        }
        
        return y;
    }
    
    private int wrapZ(int z) {
        if(z < 0) {
            z += nz;
        } else if(z >= nz) {
            z -= nz;
        }
        
        return z;
    }
    
}

// Specify a reaction, complete with reactant(s), product(s), enzyme and an additional rate constant.
class reactionParameters {
    int substrate; int product; int enzyme; float rate;
    
    int s2 = -1; int p2 = -1;
    
    // Vanilla GS with one reactant and one product.
    public reactionParameters(int is, int ip, int ien, float ir) {
        substrate=is; product=ip; rate=ir; enzyme=ien;
    }
    
    // GS with two reactants and two products (with what I am experimenting at the moment).
    public reactionParameters(int is, int is2, int ip, int ip2, int ien, float ir) {
        substrate=is; product=ip; rate=ir; enzyme=ien;
        
        s2 = is2; p2 = ip2;
    }
}
