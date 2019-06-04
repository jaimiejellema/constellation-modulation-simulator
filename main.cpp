#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <random>
#include <cstdio>
#include <string>

#include <fstream>
#include <climits>

const uint32_t CONST_SIZE_POWER = 8;
const uint32_t CONST_SIZE = 1 << CONST_SIZE_POWER;

const uint64_t SIMULATION_BATCH_SIZE_POWER = 8;
const uint64_t SIMULATION_BATCH_SIZE = (1 << SIMULATION_BATCH_SIZE_POWER);
const uint64_t SIMULATION_BATCH_COUNT = 1'000;

#define CALC_TRANSITION_PROB
#define ITERATE_SNR

typedef double coord_t;

struct treeNode {
    coord_t I;
    coord_t Q;
    uint32_t index;
    treeNode* leave[2];
};

treeNode* rootNode;
treeNode constellation[CONST_SIZE];
treeNode constellationJellema[CONST_SIZE];
treeNode constellationLarsson[CONST_SIZE];

uint32_t lookupTableJellema[CONST_SIZE];

inline double dist(struct treeNode *a, struct treeNode *b) {
    return (a->I - b->I) * (a->I - b->I) + (a->Q - b->Q) * (a->Q - b->Q);
}

bool compareI (treeNode const& lhs, treeNode const& rhs) {
    return (lhs.I < rhs.I);
}


bool compareQ (treeNode const& lhs, treeNode const& rhs) {
    return (lhs.Q < rhs.Q);
}

bool compareRadial(treeNode const& lhs, treeNode const& rhs) {
    return ( (lhs.I * lhs.I) + (lhs.Q * lhs.Q) ) < ( (rhs.I * rhs.I) + (rhs.Q * rhs.Q) );
}

bool compareAngle(treeNode const& lhs, treeNode const& rhs) {
    return atan2(lhs.Q, lhs.I)  < atan2(rhs.Q, rhs.I);
}

bool compareCross(treeNode const& lhs, treeNode const& rhs) {
    return (abs(lhs.I) +  abs(lhs.Q))  < (abs(rhs.I) +  abs(rhs.Q));
}

bool compareCrossInverse(treeNode const& lhs, treeNode const& rhs) {
    //TODO: pretty sure this is not how I intended it to work
    return (abs(lhs.Q) -  abs(lhs.I))  < (abs(rhs.Q) -  abs(rhs.I));
}


/*
 * Will generate the constellation and put the corresponding constellation points in the constellation array
 */

void createConstellation() {
    double power = 1.0;

    double golden_angle = M_PI * (3 - sqrt(5));
    double scalingFactor = sqrt( (2 * power) / (CONST_SIZE + 1)); //Yes that +1 is intentional don't ask me why, ask Larsson

    //Note 1 up to and INCLUDING CONST_SIZE to prevent multiplication with 0
    for (uint32_t i = 1; i <= CONST_SIZE; i++) {
        constellation[i - 1].I = (scalingFactor * sqrt(i) * cos(golden_angle * i));
        constellation[i - 1].Q = (scalingFactor * sqrt(i) * sin(golden_angle * i));
        constellation[i - 1].index = i - 1;
    }
}

/*
 * Since I will probably forget this.
 * While building tree, alternate axis. I and Q
 *
 * Also remember that the following diagram needs to hold
 *
 * 0 |___|___|  0
 * 1 |___|___|  1
 *     0   1
 *
 * NOT like this
 * 0 |___|___|  1
 * 1 |___|___|  0
 *     0   1
 *
 * Depending on if two parents higher is either 0 or 1, division is either 01 or 10 respectively.
 */

treeNode* buildTree(treeNode nodes[], int depth, bool invertI, bool invertQ) {
    size_t size = CONST_SIZE >> depth;
    size_t halfWayIndex = (size - 1) / 2;

    treeNode *node = new treeNode();


    if (depth % 2 == 0) {
        std::sort(nodes, nodes + size, &compareI);
    } else {
        std::sort(nodes, nodes + size, &compareQ);
    }


    /*
    if (depth < (CONST_SIZE_POWER/2) ) {
        std::sort(nodes, nodes + size, &compareI);
    } else {
        std::sort(nodes, nodes + size, &compareQ);
    }
    */

    /*
    if (depth < 2) {
        std::sort(nodes, nodes + size, &compareRadial);
    } else {
        std::sort(nodes, nodes + size, &compareAngle);
    }
    */


    /*
    if (depth % 4 == 0 || depth % 4 == 1) {
        if (depth % 2 == 0) {
            std::sort(nodes, nodes + size, &compareI);
        } else {
            std::sort(nodes, nodes + size, &compareQ);
        }
    } else {
        if (depth % 2 == 0) {
            std::sort(nodes, nodes + size, &compareCross);
        } else {
            std::sort(nodes, nodes + size, &compareCrossInverse);
        }
    }
    */



    //when the tree is build, the leaves are the actual constellation points.
    if (size == 1) {
        node->leave[0] = nullptr;
        node->leave[1] = nullptr;
        return &nodes[0];
    }

    //Will probably mess up my pointers....
    treeNode* leftNodes  = nodes;
    treeNode* rightNodes = nodes + halfWayIndex + 1;


    if (depth % 2 == 0) {
        if (invertI) {
            node->leave[1] = buildTree(leftNodes,  depth + 1, invertI, invertQ);
            node->leave[0] = buildTree(rightNodes, depth + 1, !invertI, invertQ);
        } else {
            node->leave[0] = buildTree(leftNodes,  depth + 1, !invertI, invertQ);
            node->leave[1] = buildTree(rightNodes, depth + 1, invertI, invertQ);
        }
    } else {
        if (invertQ) {
            node->leave[1] = buildTree(leftNodes,  depth + 1, invertI, invertQ);
            node->leave[0] = buildTree(rightNodes, depth + 1, invertI, !invertQ);
        } else {
            node->leave[0] = buildTree(leftNodes,  depth + 1, invertI, !invertQ);
            node->leave[1] = buildTree(rightNodes, depth + 1, invertI, invertQ);
        }
    }

    return node;
}


void numberTree(treeNode* node) {

    //std::cout << "0 index: " << ((node->index << 1) | 0) << std::endl;
    //std::cout << "1 index: " << ((node->index << 1) | 1) << std::endl;


    if (node->leave[0]) {
        node->leave[0]->index = ((node->index << 1) | 0);
        numberTree(node->leave[0]);
    }

    //std::cout << "Node index: " << node->index << std::endl;

    if (node->leave[1]) {
        node->leave[1]->index = ((node->index << 1) | 1);
        numberTree(node->leave[1]);
    }
}

/*
 * Perform a nearestNeighbourhood search using the constellation tree that was build.
 */
uint32_t nearestNeightbour(treeNode nodes[], coord_t I, coord_t Q) {
    double currentBestDistance = 10;
    uint32_t indexBest = 0;

    for (uint32_t i = 0; i < CONST_SIZE; i++) {
        double diffI = nodes[i].I - I;
        double diffQ = nodes[i].Q - Q;
        double distance =  diffI * diffI + diffQ * diffQ;
        if (distance < currentBestDistance) {
            indexBest = i;
            currentBestDistance = distance;
        }
    }

    return indexBest;
}

void generateLookup(treeNode nodes[]) {
    for (uint32_t i = 0; i < CONST_SIZE; i++) {
        lookupTableJellema[i] = nearestNeightbour(nodes, constellation[i].I, constellation[i].Q);
    }
}


//Thanks Brian!
int hammingDist(unsigned int x, unsigned int y) {
    int dist = 0;
    unsigned int val = x^y;
    while(val) {
        ++dist;
        val &= val - 1;
    }
    return dist;
}


/*
 * Be carefull this is NOT the mapping from Larsson to me, this was just a testing thing
 */

void printTree(treeNode* node) {
    if (node->leave[0]) {
        printTree(node->leave[0]);
    }

    if (node->leave[1]) {
        printTree(node->leave[1]);
    }

    if (node->leave[0] == nullptr) {
        std::cout << std::bitset<CONST_SIZE_POWER>(node->index) << " Original " << constellation[node->index].I << " " << constellation[node->index].Q << std::endl;
        std::cout << std::bitset<CONST_SIZE_POWER>(node->index) << " Tree     " << node->I << " " << node->Q << std::endl;
    }
}

void printConstellation(treeNode nodes[]) {
    for (uint32_t i = 0; i < CONST_SIZE; i++) {
        std::cout << nodes[i].I << ","<<nodes[i].Q << "," << nodes[i].index << std::endl;
    }
}

double calculateStdDev(std::vector<uint64_t> data) {
    double mean = 0;
    double meanSquareError = 0.0;

    mean = accumulate(data.begin(), data.end(), (uint64_t) 0) / (double) data.size();

    for (uint32_t i = 0; i < data.size(); i++) {
        meanSquareError += pow(data[i] - mean, 2);
    }

    return sqrt(meanSquareError / data.size());
}

int main() {
    createConstellation();
    std::cout << "creating mapping" << std::endl;
    std::copy(std::begin(constellation), std::end(constellation), std::begin(constellationJellema));
    std::copy(std::begin(constellation), std::end(constellation), std::begin(constellationLarsson));

    std::cout << "building tree" << std::endl;
    rootNode = buildTree(constellationJellema, 0, false, false);
    rootNode->index = 0;
    std::cout << "numbering nodes" << std::endl;
    numberTree(rootNode);

    //printConstellation(constellationJellema);

    generateLookup(constellationJellema);

    //std::cout << "printing nodes" << std::endl;
    //printTree(rootNode);

#ifdef ITERATE_SNR
    int begin = -5;
    int end = 35;
    double interval = 0.1;

    int count = (end - begin) / interval;

    for (int i = 0; i <= count; i++) {
        double SNR_dB = begin + interval * i;
        double SNR = pow(10, SNR_dB / 10);

#else
        double SNR_dB = 5;
        double SNR = pow(10, SNR_dB / 10);
#endif
        std::cout << "simulating SNR[dB]: " << SNR_dB << std::endl;

        double mean = 0.0;
        double stddev = 1.0 / SNR;
        std::mt19937_64 generator(std::random_device{}());
        std::normal_distribution<double> distributionNormal(mean, stddev);
        std::uniform_int_distribution<uint32_t> distributionUniform(0, CONST_SIZE - 1);

        //Simulate a ton of random data
        std::vector<uint64_t> errorCountSymbol;
        std::vector<uint64_t> errorCountBitLarsson;
        std::vector<uint64_t> errorCountBitJellema;

        errorCountSymbol.reserve(SIMULATION_BATCH_COUNT);
        errorCountBitLarsson.reserve(SIMULATION_BATCH_COUNT);
        errorCountBitJellema.reserve(SIMULATION_BATCH_COUNT);


#ifndef CALC_TRANSITION_PROB

        for (uint64_t batch = 0; batch < SIMULATION_BATCH_COUNT; batch++) {

            //Include progress print statements if it's going to take a while
            if ((SIMULATION_BATCH_SIZE * SIMULATION_BATCH_COUNT) > 1'000'000) {
                //std::cout << "Running Batch: " << batch << " out of " << SIMULATION_BATCH_COUNT << std::endl;
            }

            uint64_t symbolError = 0;
            uint64_t bitErrorLarsson = 0;
            uint64_t bitErrorJellema = 0;
            for (uint64_t i = 0; i < SIMULATION_BATCH_SIZE; i++) {
                uint32_t nodeIndex = distributionUniform(generator);

                treeNode p = constellation[nodeIndex];

                double receivedI = p.I + distributionNormal(generator);
                double receivedQ = p.Q + distributionNormal(generator);

                uint32_t closests = nearestNeightbour(constellation, receivedI, receivedQ);

                if (closests != nodeIndex) {
                    symbolError++;
                    bitErrorLarsson += hammingDist(closests, nodeIndex);
                    bitErrorJellema += hammingDist(lookupTableJellema[closests], lookupTableJellema[nodeIndex]);
                }
            }

            errorCountSymbol.push_back(symbolError);
            errorCountBitLarsson.push_back(bitErrorLarsson);
            errorCountBitJellema.push_back(bitErrorJellema);
        }

        double symbolErrorMean =
                accumulate(errorCountSymbol.begin(), errorCountSymbol.end(), 0.0) / errorCountSymbol.size();
        double bitErrorLarssonMean =
                accumulate(errorCountBitLarsson.begin(), errorCountBitLarsson.end(), 0.0) / errorCountBitLarsson.size();
        double bitErrorJellemaMean =
                accumulate(errorCountBitJellema.begin(), errorCountBitJellema.end(), 0.0) / errorCountBitJellema.size();

        double symbolErrorStdDev = calculateStdDev(errorCountSymbol);
        double bitErrorLarssonStdDev = calculateStdDev(errorCountBitLarsson);
        double bitErrorJellemaStdDev = calculateStdDev(errorCountBitJellema);

        /*
        printf("Symbol errors:     %20.3f %10.3f %10.3f%%\n", symbolErrorMean, symbolErrorStdDev,
               (symbolErrorMean * 100.0) / SIMULATION_BATCH_SIZE);
        printf("Bit error Larsson: %20.3f %10.3f %10.3f%% %10.3f\n", bitErrorLarssonMean, bitErrorLarssonStdDev,
               (bitErrorLarssonMean * 100.0) / (SIMULATION_BATCH_SIZE * CONST_SIZE_POWER),
               (bitErrorLarssonStdDev * 100.0) / (SIMULATION_BATCH_SIZE * CONST_SIZE));
        printf("Bit error Jellema: %20.3f %10.3f %10.3f%% %10.3f\n", bitErrorJellemaMean, bitErrorJellemaStdDev,
               (bitErrorJellemaMean * 100.0) / (SIMULATION_BATCH_SIZE * CONST_SIZE_POWER),
               (bitErrorJellemaStdDev * 100.0) / (SIMULATION_BATCH_SIZE * CONST_SIZE));
       */
        printf("%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f\n",
                12, SNR_dB,
                12, symbolErrorMean / SIMULATION_BATCH_SIZE,
                12, bitErrorLarssonMean / (SIMULATION_BATCH_SIZE * CONST_SIZE_POWER),
                12, bitErrorJellemaMean / (SIMULATION_BATCH_SIZE * CONST_SIZE_POWER),

                12, symbolErrorMean,
                12, symbolErrorStdDev,
                12, bitErrorLarssonMean,
                12, bitErrorLarssonStdDev,
                12, bitErrorJellemaMean,
                12, bitErrorJellemaStdDev
        );

#ifdef ITERATE_SNR
    }
#endif

#else

        std::ofstream probabilities;
        //std::ofstream mapping;

        char snr_string[256];

        sprintf(snr_string, "%.3f", SNR_dB);

        std::string filename = "test_" + std::to_string(CONST_SIZE_POWER) + std::string() + "__" + snr_string;

        probabilities.open(filename, std::ios::out | std::ios::binary);
        //mapping.open("mapping_8_test.bin", std::ios::out | std::ios::binary);

        std::vector<uint32_t> receivedSymbol;
        receivedSymbol.reserve(CONST_SIZE);

        for (uint32_t i = 0; i < CONST_SIZE; i++) {
            receivedSymbol.push_back(0);
            //mapping.write(reinterpret_cast<char*>(&lookupTableJellema[i]), sizeof(uint32_t));
        }



        for (uint64_t symbol = 0; symbol < CONST_SIZE; symbol++) {
            std::cout << "Running symbol: " << (symbol + 1) << " out of " << CONST_SIZE << std::endl;

            std::fill(receivedSymbol.begin(), receivedSymbol.end(), 0);

            for (uint64_t i = 0; i < SIMULATION_BATCH_SIZE; i++) {
                treeNode p = constellation[symbol];

                double receivedI = p.I + distributionNormal(generator);
                double receivedQ = p.Q + distributionNormal(generator);

                uint32_t closests = nearestNeightbour(constellation, receivedI, receivedQ);

                receivedSymbol[closests] += 1;
            }


            for (uint32_t j = 0; j < CONST_SIZE; j++) {
                receivedSymbol[j] = receivedSymbol[j] * (UINT32_MAX / SIMULATION_BATCH_SIZE);
            }

            probabilities.write(reinterpret_cast<char *>(receivedSymbol.data()),
                                receivedSymbol.size() * sizeof(uint32_t));
        }

        probabilities.close();

        //mapping.close();

#ifdef ITERATE_SNR
    }
#endif


#endif

    return 0;
}