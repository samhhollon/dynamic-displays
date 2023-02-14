/******************************************************************************
 Good-genes model of sexual selection on male courtship effort, based on
 'revealing handicap' mechanism as described by van Doorn & Weissing (2006).
 Males with high genetic viability enjoy greater marginal fecundity benefits
 of displaying than do males with low genetic viability. Females have evolvable
 preferences for fixed and flexible display elements. Males can its display
 investment depending on the number of rivals displaying.
******************************************************************************/


//*****************************************************************************
// HEADER FILES
#include "model.hpp"


//*****************************************************************************
// MAIN FUNCTION

int main() {
    // Get the current time.
    auto timeStart = chrono::system_clock::now();

    // Set trait parameters based on which displays (fixed and/or flexible)
    // are turned on.
    setTraitParams();

    // Create writer for the summary statistics output file.
    Writer mainWriter(summaryOutput);
    mainWriter.writeMsg("replicate,generation,meanPfix,meanPflex,"
        "meanTrait,meanAlpha,meanBeta,meanGamma,meanFlexDisplay,"
        "meanV,livingFemales,livingMales\n");

    // Create writer for the column headers in the census output file.
    Writer mainWriter2(censusOutput);
    mainWriter2.writeMsg("replicate,generation,v,t,alpha,beta,gamma,meanFlex\n");

    // Create writer for the column headers in the parameter output file.
    Writer mainWriter3(paramOutput);
    mainWriter3.writeMsg("replicate,N,mutGenP,mutGenSlope,numGen,"
        "skipSummary,skipCensus,fixedDisplay,flexibleDisplay,"
        "displaySlope,initDisplay,initP,meanPoisson,mutDisplay,"
        "mutP,stdevDisplay,stdevP,mutV,stdevV,biasV,b,k,c,Vopt\n");

    // Create writer for the column headers in the display output file.
    Writer mainWriter4(displayOutput);
    mainWriter4.writeMsg("replicate,generation,v,t,t0,t1,alpha,beta,a0,a1,b0,b1,"
        "g0,g1,flexibleDisplay,v0,v1,numComp\n");

    // Create a vector of promises to await the return value of each simulation run.
    // The length of this vector is the product of the lengths of each of the vector
    // of parameter values iterated through during the simulation batch run.
    int cArrSize = *(&cArr + 1) - cArr;
    int kArrSize = *(&kArr + 1) - kArr;
    int numRuns = cArrSize * kArrSize * numRep;
    std::vector<std::future<void>> results(numRuns);

    // Create a thread pool.
    ctpl::thread_pool pool(num_threads);

    // Track how many runs have started with a counter.
    int counter = 0;

    // Iterate through each value of stdevFecArr.
    // Iterate through each value of cArr.
    for (const double &cVal : cArr) {

        // Iterate through each value of kArr.
        for (const double &kVal : kArr) {
                
            // Iterate through each replicate.
            for (int rep = 0; rep < numRep; ++rep) {

                // Push a simulation run to the thread pool.
                results[counter] = pool.push([counter, cVal, kVal, numRuns](int) {
                    // Create a model.
                    Model m(counter, cVal, kVal);

                    // Run the model and get its return value to calculate the number of runs
                    // completed so far.
                    int runsComplete = m.run();

                    // Calculate the percentage of runs complete.
                    int percentComplete = round(100 * runsComplete / numRuns);

                    // Print the number of runs completed, the total number to be run, and the
                    // percentage completed.
                    cout << runsComplete << " of " << numRuns
                        << " runs (" << percentComplete << "%) complete\n";
                });

                // Increment the counter
                ++counter;
            }
        }
    }

    // // Get each simulation's return value.
    for (int i = 0; i < counter; ++i) {
        results[i].get();
    }

    // Get the current time.
    auto timeEnd = chrono::system_clock::now();

    // Calculate the time the program took to run.
    chrono::duration<double> timeRun = timeEnd - timeStart;

    // Convert the runtime to hours.
    double runTimeHours = round(100 * timeRun.count() / 3600) / 100;

    // Get the number of minutes per model run (replicate).
    double replicateMinutes = round(100 * timeRun.count() / (60 * counter)) / 100;

    // Output the time the program took to run.
    cout << "The computation finished in "
         << runTimeHours << " hours "
         << "(" << replicateMinutes
         << " minutes per run).\n";
}