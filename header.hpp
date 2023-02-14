/******************************************************************************
 Store constants, types, and variables for use in other files.
******************************************************************************/

//*****************************************************************************
// HEADER FILES
#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
#include <mutex>
#include <random>
#include <memory>
#include <thread>
#include <chrono>
#include "ctpl_stl.h"


//*****************************************************************************
// NAMESPACE
using namespace std;


//*****************************************************************************
// CONSTANTS

// Constants specifying the number and structure of runs.
const int numRep = 20;          // number of replicates
const int N = 5000;             // population size
const int mutGenP = 5000;       // number of generations after which preference mutates from initial value
const int mutGenSlope = 10000;  // number of generations after which the flexible display slope/display when alone mutates if turned on
const int numGen = 40000;       // number of generations
const int skipSummary = 100;    // interval between print-outs
const int skipCensus = 20000;   // interval between census printouts

// Constants specifying the display traits to include and what variable
// the flexible display strategy is based on.
const bool fixedDisplay = false;    // fixed display trait
const bool flexibleDisplay = true;  // flexible display trait
const bool displaySlope = false;    // whether flexible display's slope and display when alone are turned on

// Constants specifying the initial values of the traits turned on
const double initDisplay = 0.0;  // starting value of all displays turned on
const double initP = 1.0;        // starting value of all preferences turned on

// Constants specifying the variation in the number of males that compete
// to mate with each female.
const int meanPoisson = 3;    // mean if Poisson distributed. Minimum 1.

// Constants specifying mutation rates and standard deviations for displays
// and preferences when a given display is turned on.
const double mutDisplay = 0.05;    // mutation rate of displays when turned on
const double mutP = 0.05;          // mutation rate of preferences when turned on
const double stdevDisplay = 0.1;   // mutation stdev of displays when turned on (0.01)
const double stdevP = 0.1;         // mutation stdev of preferences when turned on (0.01)

// Constants specifying the mutation rates and stdevs that don't depend on
// which traits are turned on.
const double stdevPfix = stdevP;          // mutation stdev in pref. for the fixed display
const double stdevPflex = stdevP;         // mutation stdev in pref. intercept for the flexible display
const double stdevAlpha = stdevDisplay;   // mutation stdev in flexible-display intercept
const double stdevBeta = stdevDisplay;    // mutation stdev in flexible-display slope
const double stdevGamma = stdevDisplay;   // mutation stdev in flexible display when alone
const double stdevTrait = stdevDisplay;   // mutation stdev in fixed display
const double mutV = 0.10;                 // mutation rate in viability v
const double stdevV = 0.5;                // mutation stdev in v
const double biasV = 0.95;                // probability that mutation in v is downwards (0.5 = unbiased)

// Constants specifying preference costs
const double b = 0.005;                  // cost of each preference component
const double kArr[] = {0.175};           // cost of fixed display
const double cArr[] = {0.032};           // cost of flexible display
const double Vopt = 0.0;                 // optimal value of v

// Constants specifying initial trait values that don't depend on which traits
// are turned on.
const double initV = Vopt;  // starting value of v

// Number of threads in the thread pool
const int num_threads = thread::hardware_concurrency();


//*****************************************************************************
// TRAIT PARAMETERS

// Variables for fixed display.
double mutPfix = 0.0;    // mutation rate in pref. for the fixed display
double mutTrait = 0.0;   // mutation rate in fixed trait
double initPfix = 0.0;   // starting value of female preference for fixed display
double initTrait = 0.0;  // starting value of fixed trait

// Variables for flexible display.
double mutPflex = 0.0;   // mutation rate in pref. for the flexible display
double mutAlpha = 0.0;   // mutation rate in alpha
double mutBeta = 0.0;    // mutation rate in beta
double mutGamma = 0.0;   // mutation rate in beta
double initPflex = 0.0;  // starting female preference for flexible display
double initAlpha = 0.0;  // starting value of alpha
double initBeta = 0.0;   // starting value of beta
double initGamma = 0.0;  // starting value for gamma

// Set the trait parameters below based on which display traits have been
// turned on using constants fixed and flexible.
void setTraitParams() {
    // Fixed display
    if (fixedDisplay) {
        // Mutation rates
        mutPfix = mutP;
        mutTrait = mutDisplay;

        // Initial values
        initPfix = initP;
        initTrait = initDisplay;
    }

    // Flexible display
    if (flexibleDisplay) {
        // Intercept
        mutPflex = mutP;
        mutAlpha = mutDisplay;
        initPflex = initP;
        initAlpha = initDisplay;

        // Display slope and display when alone
        if(displaySlope) {
            mutBeta = mutDisplay;
            initBeta = initDisplay;
            mutGamma = mutDisplay;
            initGamma = initDisplay;
        } else {
            mutBeta = 0.0;
            initBeta = 0.0;
            mutGamma = initAlpha;
            initGamma = initAlpha;
        }
    }
}


//*****************************************************************************
// RANDOM NUMBERS

// Thread-safe random number generator (continuous uniform, [0,1])
double uniform() {
    static thread_local mt19937* generator = nullptr;
    if (!generator) generator = new mt19937(clock() + hash<thread::id>()(this_thread::get_id()));
    uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(*generator);
}

// Thread-safe random number generator (Poisson)
double poisson(const int &mean) {
    static thread_local mt19937* generator = nullptr;
    if (!generator) generator = new mt19937(clock() + hash<thread::id>()(this_thread::get_id()));
    poisson_distribution<int> distribution(mean);
    return distribution(*generator);
}

// Thread-safe random number generator (normal)
double normal(const double &mean, const double &stdev) {
    static thread_local mt19937* generator = nullptr;
    if (!generator) generator = new mt19937(clock() + hash<thread::id>()(this_thread::get_id()));
    normal_distribution<double> distribution(mean, stdev);
    return distribution(*generator);
}


//*****************************************************************************
// TYPES AND CLASSES

// Define individual traits
typedef struct Individual {
    double Pfix0, Pfix1, Pflex0, Pflex1,         // genetic values for female preferences
        A0, A1, B0, B1, G0, G1, T0, T1, V0, V1,  // genetic values for male strategy and viability
        pFix, pFlex, alpha, beta, gamma, t, v,   // phenotypic values
        attract, samplings;                      // phenotypic values
    int livingIndex;                             // number of offspring and sampling events
    bool alive;                                  // whether male is alive
} Population[N/2];


// Define a synchronized file class. Each thread will write to the same
// output file, so we use a mutex to control access.
class SynchronizedFile {
    public:
        SynchronizedFile(const string& path) : _path(path) {
            // Open file for writing
            out.open(_path);

            // Initialize the number of rows at 0.
            rows = 0;
        }

        void addRow(const string& dataToWrite) {
            // Ensure that only one thread can execute at a time
            lock_guard<mutex> lock(_writerMutex);

            // Write to the file
            out << dataToWrite;

            // Increment rows
            ++rows;
        }

        int getRows() {
            // Ensure that only one thread can execute at a time
            lock_guard<mutex> lock(_writerMutex);

            return rows;
        }

    private:
        ofstream out;
        string _path;
        mutex _writerMutex;
        int rows;
};


// Create a writer class that writes to a synchronized file
class Writer {
    public:
        Writer(shared_ptr<SynchronizedFile> sf) : _sf(sf) {}

        void writeMsg(string msg) {
            _sf -> addRow(msg);
        }

        int readRows() {
            return _sf -> getRows();
        }
    private:
        shared_ptr<SynchronizedFile> _sf;
};


//*****************************************************************************
// OUTPUT FILES

// Create a synchronized file to store summary statistics
auto summaryOutput = make_shared<SynchronizedFile>("2.data/summaryOutput.csv");

// Create a synchronized file to store snapshots of male traits.
auto censusOutput = make_shared<SynchronizedFile>("2.data/censusOutput.csv");

// Create a synchronized file to store parameter values
auto paramOutput = make_shared<SynchronizedFile>("2.data/paramOutput.csv");

// Create a synchronized file to store display events.
auto displayOutput = make_shared<SynchronizedFile>("2.data/displayOutput.csv");