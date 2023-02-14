/******************************************************************************
 Define the Model class. This class implements a single simulation run.
******************************************************************************/

//*****************************************************************************
// HEADER FILES
#include "header.hpp"

//*****************************************************************************
// CLASS DEFINITION

class Model {
    // Female and male parent and offspring generations.
    Population f0, f1, m0, m1;

    // Current replicate, generation, and number of living males.
    int replicate, generation, livingFemales, livingMales;

    // Cost parameters, cumulative sums, and summary statistics.
    double c, k, meanPfix, meanPflex, meanV, meanAlpha, meanBeta,
        meanGamma, meanTrait, meanFlexDisplay;

    // Writer to write to the displays output file.
    Writer displayWriter;

    public:


    //*****************************************************************************
    // CONSTRUCTOR
    Model(double Replicate, double C, double K) : displayWriter(displayOutput) {
        // Simulation parameters.
        replicate = Replicate;
        c = C;
        k = K;
        generation = 0;

        // Initialize the genetic values of each individual in the population.
        for (int i = 0; i < N / 2; ++i) {

            // Preference for fixed display
            f0[i].Pfix0 = m0[i].Pfix0 = initPfix;
            f0[i].Pfix1 = m0[i].Pfix1 = initPfix;

            // Preference for flexible display
            f0[i].Pflex0 = m0[i].Pflex0 = initPflex;
            f0[i].Pflex1 = m0[i].Pflex1 = initPflex;

            // Fixed display
            f0[i].T0 = m0[i].T0 = initTrait;
            f0[i].T1 = m0[i].T1 = initTrait;

            // Flexible display intercept
            f0[i].A0 = m0[i].A0 = initAlpha;
            f0[i].A1 = m0[i].A1 = initAlpha;

            // Flexible display slope
            f0[i].B0 = m0[i].B0 = initBeta;
            f0[i].B1 = m0[i].B1 = initBeta;

            // Flexible display when alone
            f0[i].G0 = m0[i].G0 = initGamma;
            f0[i].G1 = m0[i].G1 = initGamma;

            // Viability
            f0[i].V0 = m0[i].V0 = initV;
            f0[i].V1 = m0[i].V1 = initV;
        }

        // Calculate phenotypes
        phenotypes();
    }

    //*****************************************************************************
    // CALCULATE PHENOTYPES AND RESCALE FITNESS

    void phenotypes() {

        // Initialize the number of living females and males at 0.
        livingFemales = 0;
        livingMales = 0;

        // Calculate female trait values.
        for (int i = 0; i < N/2; ++i)  { // for all females

            // Take the mean of the parental genes.
            f0[i].pFix = 0.5 * (f0[i].Pfix0 + f0[i].Pfix1);    // preference for the fixed display
            f0[i].pFlex = 0.5 * (f0[i].Pflex0 + f0[i].Pflex1); // preference for the flexible display
            f0[i].v = 0.5 * (f0[i].V0 + f0[i].V1);             // genetic viability

            // Female viability decreases with expressed preference and distance
            // from genetic viability optimum.
            double viab = exp(
                - b * f0[i].pFix * f0[i].pFix
                - b * f0[i].pFlex * f0[i].pFlex
                - abs(Vopt - f0[i].v));

            // The individual is alive with probability = viab.
            if(uniform() < viab) {
                ++livingFemales;
            }
            f0[i].livingIndex = livingFemales;
        }

        // Calculate male trait values.
        for (int i = 0; i < N/2; ++i) {  // for all males

            // Take the mean of the parental genes.
            m0[i].t = max(0.5 * (m0[i].T0 + m0[i].T1), 0.0);  // investment in fixed display
            m0[i].alpha = 0.5 * (m0[i].A0 + m0[i].A1);        // intercept of flexible display
            m0[i].beta = 0.5 * (m0[i].B0 + m0[i].B1);         // slope of flexible display
            m0[i].gamma = 0.5 * (m0[i].G0 + m0[i].G1);        // flexible display when alone
            m0[i].v = 0.5 * (m0[i].V0 + m0[i].V1);            // genetic viability

            // Male viability depends on expressed trait and distance from Vopt. The
            // expressed trait cannot be negative.
            double viab = exp(
                - k * m0[i].t * m0[i].t
                - abs(Vopt - m0[i].v));

            // The individual is alive with probability = viab.
            if(uniform() < viab) {
                m0[i].alive = true;
                ++livingMales;
            } else {
                m0[i].alive = false;
            }
            m0[i].livingIndex = livingMales;
        }
    }


    //*****************************************************************************
    // MUTATE GENOTYPE

    // Helper function. Mutates a single gene. The bias parameter is the
    // probability that the mutation is downward (negative).
    inline double mutateGene(double rate, double stdev, double bias = 0.5) {
        double mutation = 0.0;
        if (uniform() < rate)
            uniform() < bias ? mutation = - abs(normal(0.0, stdev)) : mutation = abs(normal(0.0, stdev));
        return(mutation);
    }

    void mutate(Individual &kid) {
        // Mutate pFix
        if (generation > mutGenP) {  // fixed from generation 0 to mutGenP, evolves thereafter
            // Mutate fixed preference
            kid.Pfix0 += mutateGene(mutPfix, stdevPfix);
            kid.Pfix1 += mutateGene(mutPfix, stdevPfix);

            // Mutate flexible preference
            kid.Pflex0 += mutateGene(mutPflex, stdevPflex);
            kid.Pflex1 += mutateGene(mutPflex, stdevPflex);
        }

        // Mutate viability
        kid.V0 += mutateGene(mutV, stdevV, biasV);
        kid.V1 += mutateGene(mutV, stdevV, biasV);

        // Mutate fixed trait
        kid.T0 += mutateGene(mutTrait, stdevTrait);
        kid.T1 += mutateGene(mutTrait, stdevTrait);

        // Mutate Alpha
        kid.A0 += mutateGene(mutAlpha, stdevAlpha);
        kid.A1 += mutateGene(mutAlpha, stdevAlpha);
        
        // Fixed from generation 0 to mutGenSlope, evolves thereafter
        if (generation > mutGenSlope) {  
            // Mutate Beta
            kid.B0 += mutateGene(mutBeta, stdevBeta);
            kid.B1 += mutateGene(mutBeta, stdevBeta);

            // Mutate Gamma
            kid.G0 += mutateGene(mutGamma, stdevGamma);
            kid.G1 += mutateGene(mutGamma, stdevGamma);
        }
        
        // If the display when alone is turned off, then set equal to the
        // intercept.
        if(generation <= mutGenSlope || !displaySlope) {
            kid.G0 = kid.A0;
            kid.G1 = kid.A1;
        }
    }


    //*****************************************************************************
    // PRODUCE OFFSPRING

    void createkid(int mother, int father, Individual &kid) {
        // Preference for fixed display
        uniform() < 0.5 ? kid.Pfix0 = f0[mother].Pfix0
                        : kid.Pfix0 = f0[mother].Pfix1;
        uniform() < 0.5 ? kid.Pfix1 = m0[father].Pfix0
                        : kid.Pfix1 = m0[father].Pfix1;

        // Preference for flexible display
        uniform() < 0.5 ? kid.Pflex0 = f0[mother].Pflex0
                        : kid.Pflex0 = f0[mother].Pflex1;
        uniform() < 0.5 ? kid.Pflex1 = m0[father].Pflex0
                        : kid.Pflex1 = m0[father].Pflex1;

        // Fixed display
        uniform() < 0.5 ? kid.T0 = f0[mother].T0 : kid.T0 = f0[mother].T1;
        uniform() < 0.5 ? kid.T1 = m0[father].T0 : kid.T1 = m0[father].T1;

        // Flexible display intercept
        uniform() < 0.5 ? kid.A0 = f0[mother].A0 : kid.A0 = f0[mother].A1;
        uniform() < 0.5 ? kid.A1 = m0[father].A0 : kid.A1 = m0[father].A1;

        // Flexible display slope
        uniform() < 0.5 ? kid.B0 = f0[mother].B0 : kid.B0 = f0[mother].B1;
        uniform() < 0.5 ? kid.B1 = m0[father].B0 : kid.B1 = m0[father].B1;

        // Flexible display when alone
        uniform() < 0.5 ? kid.G0 = f0[mother].G0 : kid.G0 = f0[mother].G1;
        uniform() < 0.5 ? kid.G1 = m0[father].G0 : kid.G1 = m0[father].G1;
        
        // Viability
        uniform() < 0.5 ? kid.V0 = f0[mother].V0 : kid.V0 = f0[mother].V1;
        uniform() < 0.5 ? kid.V1 = m0[father].V0 : kid.V1 = m0[father].V1;
    }


    //*****************************************************************************
    // FEMALE CHOOSES MATE

    int choose(double pFix, double pFlex) {

        // Determine the number of competing males.
        int numComp = poisson(meanPoisson) + 1;
        
        int Candidates[numComp], father;
        double DisplayTotal[numComp], DisplayFix[numComp],
            DisplayFlex[numComp], Attractiveness[numComp];

        // Select random males
        for (int i = 0; i < numComp; i++) {
            // Select from among living males with equal probability.
            int s = uniform() * livingMales;
            for(int j = 0; j < N/2; j++) {
                if(s <= m0[j].livingIndex) {
                    Candidates[i] = j;
                    break;
                }
            }
        }

        // Males choose their level of courtship investment.
        double sumAttractiveness = 0.0;
        for (int i = 0; i < numComp; ++i) {
            // Fixed part of the display.
            DisplayFix[i] = m0[Candidates[i]].t;

            // Flexible part of the display. Minimum 0.
            if(numComp > 1) {
                // If more than one male displays, each one's flexible display investment is
                // a linear function of alpha and beta centered at meanPoisson
                DisplayFlex[i] = max(m0[Candidates[i]].alpha + (numComp - 1 - meanPoisson) * m0[Candidates[i]].beta, 0.0);
            } else {
                // If one male displays alone, his flexible display investigate is gamma.
                DisplayFlex[i] = max(m0[Candidates[i]].gamma, 0.0);
            }

            // To implement a revealing handicap, multiply each component of display by
            // male's genetic viability (males with higher viability benefit more)
            Attractiveness[i] = sumAttractiveness + exp(exp(-abs(Vopt - m0[Candidates[i]].v)) *
                    (pFix * DisplayFix[i] + pFlex * DisplayFlex[i]));

            m0[i].attract += Attractiveness[i] - sumAttractiveness;
            ++m0[i].samplings;

            // Keep a cumulative total for later use
            sumAttractiveness = Attractiveness[i];

            // Males die proportional to the square of their flexible display.
            double probSurvive = exp(- c * DisplayFlex[i] * DisplayFlex[i]);
            if(uniform() > probSurvive) {
                m0[Candidates[i]].alive = false;
            }
        }

        // On the final generation, record the displays that take place.
        if (generation == numGen) {
            for (int i = 0; i < numComp; ++i) {
                displayWriter.writeMsg(displayRow(m0[Candidates[i]].v,
                    DisplayFix[i],
                    m0[Candidates[i]].T0,
                    m0[Candidates[i]].T1,
                    m0[Candidates[i]].alpha,
                    m0[Candidates[i]].beta,
                    m0[Candidates[i]].A0,
                    m0[Candidates[i]].A1,
                    m0[Candidates[i]].B0,
                    m0[Candidates[i]].B1,
                    m0[Candidates[i]].G0,
                    m0[Candidates[i]].G1,
                    DisplayFlex[i],
                    m0[Candidates[i]].V0,
                    m0[Candidates[i]].V1,
                    numComp));
            }
        }

        double r = uniform();
        for (int i = 0; i < numComp; ++i) {
            if (r * sumAttractiveness <= Attractiveness[i]) {
                father = Candidates[i];  // choose candidate male with chance
                                         // proportional to attractiveness of his display
                                         // Note: males that die from displaying may still be
                                         // selected here.
                break;
            }
        }

        // Recalculate number of living males and their indices.
        livingMales = 0;
        for(int i = 0; i < N/2; i++) {
            livingMales += m0[i].alive;
            m0[i].livingIndex = livingMales;
        }

        return(father);
    }


    //*****************************************************************************
    // PAIRING AND REPRODUCTION TO PRODUCE NEXT GENERATION

    void nextGen() {
        int mother, NumSons = 0, NumDaughters = 0;

        for(int i = 0; i < N; ++i) {

            // Randomly select a mother from the living females.
            double r = uniform() * livingFemales;
            for (int j = 0; j < N/2; ++j) {
                if (r <= f0[j].livingIndex) {
                    mother = j;
                    break;
                }
            }

            double pFix = f0[mother].pFix;    // female's preference for fixed display
            double pFlex = f0[mother].pFlex;  // female's preference for flexible

            // Pair female with male
            int father = choose(pFix, pFlex);  // select winner

            // Produce offspring
            Individual kid;
            createkid(mother, father, kid);

            // Mutate offspring
            mutate(kid);

            // Populate offspring generation
            if (i % 2 == 0) {  // determine sex (alternate to maintain 1:1 sex ratio)
                f1[NumDaughters] = kid;
                ++NumDaughters;
            } else {
                m1[NumSons] = kid;
                ++NumSons;
            }
        }
    }


    //*****************************************************************************
    // CALCULATE STATISTICS

    void statistics() {
        double sumPfix = 0.0, sumPflex = 0.0, sumAlpha = 0.0, sumBeta = 0.0,
        sumGamma = 0.0, sumTrait = 0.0, sumFlexDisplay = 0.0, sumV = 0.0;

        for (int i = 0; i < N/2; ++i) {
            // Sums
            sumPfix += 0.5 * (f0[i].Pfix0 + f0[i].Pfix1);
            sumPflex += 0.5 * (f0[i].Pflex0 + f0[i].Pflex1);
            sumAlpha += 0.5 * (f0[i].A0 + f0[i].A1);
            sumBeta += 0.5 * (f0[i].B0 + f0[i].B1);
            sumGamma += 0.5 * (f0[i].G0 + f0[i].G1);
            sumTrait += max(0.5 * (f0[i].T0 + f0[i].T1), 0.0);
            sumV +=  f0[i].v;
            sumFlexDisplay += max(0.5 * (f0[i].A0 + f0[i].A1), 0.0);
        }

        // Means
        meanPfix = sumPfix / double(N/2);
        meanPflex = sumPflex / double(N/2);
        meanV = sumV / double(N/2);
        meanAlpha = sumAlpha / double(N/2);
        meanBeta = sumBeta / double(N/2);
        meanGamma = sumGamma / double(N/2);
        meanTrait = sumTrait / double(N/2);
        meanFlexDisplay = sumFlexDisplay / double(N/2);
    }


    //*****************************************************************************
    // COMPOSE ONE ROW OF SUMMARY STATISTICS

    string dataRow() {
        string row =
            to_string(replicate).append(",").
            append(to_string(generation).append(",")).
            append(to_string(meanPfix)).append(",").
            append(to_string(meanPflex)).append(",").
            append(to_string(meanTrait)).append(",").
            append(to_string(meanAlpha)).append(",").
            append(to_string(meanBeta)).append(",").
            append(to_string(meanGamma)).append(",").
            append(to_string(meanFlexDisplay)).append(",").
            append(to_string(meanV)).append(",").
            append(to_string(double(livingFemales)/double(N/2))).append(",").
            append(to_string(double(livingMales)/double(N/2))).append("\n");
        return row;
    }


    //*****************************************************************************
    // CENSUS OF MALE VALUES FOR DISPLAY ELEMENTS AND GENETIC QUALITY

    string census() {
        double displayFlex;  // store the flexible display level of the current male

        // Start an empty string. We'll append data to this.
        string censusValues = "";

        // Add one row for each male in the population
        for (int i = 0; i < N/2; ++i) {
            // Calculate mean flexible display level.
            displayFlex = max(m0[i].alpha, 0.0);

            censusValues.append(
                to_string(replicate)).append(",").
                append(to_string(generation)).append(",").
                append(to_string(m0[i].v)).append(",").
                append(to_string(m0[i].t)).append(",").
                append(to_string(m0[i].alpha)).append(",").
                append(to_string(m0[i].beta)).append(",").
                append(to_string(m0[i].gamma)).append(",").
                append(to_string(displayFlex)).append("\n");
        }
        return censusValues;
    }


    //*****************************************************************************
    // WRITE PARAMETERS TO OUTPUT FILE

    string params() {
        string paramVals =
            to_string(replicate).append(",").
            append(to_string(N)).append(",").
            append(to_string(mutGenP)).append(",").
            append(to_string(mutGenSlope)).append(",").
            append(to_string(numGen)).append(",").
            append(to_string(skipSummary)).append(",").
            append(to_string(skipCensus)).append(",").
            append(to_string(fixedDisplay)).append(",").
            append(to_string(flexibleDisplay)).append(",").
            append(to_string(displaySlope)).append(",").
            append(to_string(initDisplay)).append(",").
            append(to_string(initP)).append(",").
            append(to_string(meanPoisson)).append(",").
            append(to_string(mutDisplay)).append(",").
            append(to_string(mutP)).append(",").
            append(to_string(stdevDisplay)).append(",").
            append(to_string(stdevP)).append(",").
            append(to_string(mutV)).append(",").
            append(to_string(stdevV)).append(",").
            append(to_string(biasV)).append(",").
            append(to_string(b)).append(",").
            append(to_string(k)).append(",").
            append(to_string(c)).append(",").
            append(to_string(Vopt)).append("\n");
        return paramVals;
    }


    //*****************************************************************************
    // RECORD ONE DISPLAY

    string displayRow(double viab, double fix, double t0, double t1, double a, double b,
        double a0, double a1, double b0, double b1, double g0, double g1, double flex,
        double v0, double v1, int competitors) {
        // Start an empty string. We'll append data to this.
        string censusValues =
            to_string(replicate).append(",").
            append(to_string(generation).append(",")).
            append(to_string(viab).append(",")).
            append(to_string(fix)).append(",").
            append(to_string(t0)).append(",").
            append(to_string(t1)).append(",").
            append(to_string(a)).append(",").
            append(to_string(b)).append(",").
            append(to_string(a0)).append(",").
            append(to_string(a1)).append(",").
            append(to_string(b0)).append(",").
            append(to_string(b1)).append(",").
            append(to_string(g0)).append(",").
            append(to_string(g1)).append(",").
            append(to_string(flex)).append(",").
            append(to_string(v0)).append(",").
            append(to_string(v1)).append(",").
            append(to_string(competitors)).append("\n");
        return censusValues;
    }


    //*****************************************************************************
    // RUN THE MODEL
    
    int run() {
        // Create writers to write to the output files
        Writer threadWriter1(summaryOutput);
        Writer threadWriter2(censusOutput);
        Writer threadWriter3(paramOutput);

        // Calculate summary statistics and write genetic trait values to
        // the main output file
        statistics();

        // Write the summary statistics to the output file
        threadWriter1.writeMsg(dataRow());

        // Simulate each generation
        for (generation = 1; generation <= numGen; ++generation) {
            // Create the next generation.
            nextGen();

            // Write output every "skip" generations
            if (generation % skipSummary == 0) {
                // Recalculate summary statistics and write genetic trait values to
                // the main output file
                statistics();
                
                // Write the summary statistics to the output file
                threadWriter1.writeMsg(dataRow());
            }

            // Take a census just before we allow mutations in
            // preferences and again on the final generation.
            if (generation % skipCensus == 0) {
                threadWriter2.writeMsg(census());
            }

            // Overwrite parental generations with offspring (female)
            for (int i = 0; i < N/2; ++i)
                f0[i] = f1[i];
            
            // Overwrite parental generations with offspring (male)
            for (int i = 0; i < N/2; ++i)
                m0[i] = m1[i];

            // Recalculate phenotypes
            phenotypes();
        }

        // Add paramater values to the census output file
        threadWriter3.writeMsg(params());

        // Determine the number of simulations that have finished based
        // on the number of lines in paramOutput.csv
        int runsComplete = threadWriter3.readRows() - 1;

        // Return runsComplete
        return runsComplete;
    }
};