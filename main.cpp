#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <random>


const unsigned int CONST_SIZE_POWER = 8;
const unsigned int CONST_SIZE = 1 << CONST_SIZE_POWER;
const unsigned int SIMULATION_BATCH_SIZE = 100'000;
const unsigned int SIMULATION_BATCH_NUMBER = 10;

typedef double coord_t;

struct treeNode {
    coord_t I;
    coord_t Q;
    unsigned int index;
    treeNode* leave[2];
    coord_t splitPlane;
};

treeNode* rootNode;
treeNode constellation[CONST_SIZE];
treeNode constellationJellema[CONST_SIZE];
treeNode constellationLarsson[CONST_SIZE];

unsigned int lookupTableJellema[CONST_SIZE];

inline double dist(struct treeNode *a, struct treeNode *b) {
    return (a->I - b->I) * (a->I - b->I) + (a->Q - b->Q) * (a->Q - b->Q);
}

bool compareI (treeNode const& lhs, treeNode const& rhs) {
    return (lhs.I < rhs.I);
}

bool compareIreverse (treeNode const& lhs, treeNode const& rhs) {
    return (lhs.I > rhs.I);
}


bool compareQ (treeNode const& lhs, treeNode const& rhs) {
    return (lhs.Q < rhs.Q);
}

bool compareQreverse (treeNode const& lhs, treeNode const& rhs) {
    return (lhs.Q > rhs.Q);
}

/*
 * Will generate the constellation and put the corresponding constellation points in the constellation array
 */

void createConstellation() {
    double golden_angle = M_PI * (3 - sqrt(5));
    double scalingFactor = 1 / sqrt(CONST_SIZE);

    //Note 1 up to and INCLUDING CONST_SIZE to prevent multiplication with 0

    for (int i = 1; i <= CONST_SIZE; i++) {
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
        node->splitPlane = (nodes[halfWayIndex].I + nodes[halfWayIndex + 1].I) / 2;
    } else {
        std::sort(nodes, nodes + size, &compareQ);
        node->splitPlane = (nodes[halfWayIndex].Q + nodes[halfWayIndex + 1].Q) / 2;
    }


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
            std::cout << "Test: 1"<< std::endl;
            node->leave[1] = buildTree(leftNodes,  depth + 1, invertI, invertQ);
            node->leave[0] = buildTree(rightNodes, depth + 1, !invertI, invertQ);
        } else {
            std::cout << "Test: 2"<< std::endl;
            node->leave[0] = buildTree(leftNodes,  depth + 1, !invertI, invertQ);
            node->leave[1] = buildTree(rightNodes, depth + 1, invertI, invertQ);
        }
    } else {
        if (invertQ) {
            std::cout << "Test: 3"<< std::endl;
            node->leave[1] = buildTree(leftNodes,  depth + 1, invertI, invertQ);
            node->leave[0] = buildTree(rightNodes, depth + 1, invertI, !invertQ);
        } else {
            std::cout << "Test: 4"<< std::endl;
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
    std::cout << "Node index: " << node->index << std::endl;

    if (node->leave[1]) {
        node->leave[1]->index = ((node->index << 1) | 1);
        numberTree(node->leave[1]);
    }
}

/*
 * Perform a nearestNeighbourhood search using the constellation tree that was build.
 */
int nearestNeightbour(treeNode nodes[], coord_t I, coord_t Q) {
    double currentBestDistance = 10;
    int indexBest = 0;

    for (int i = 0; i < CONST_SIZE; i++) {
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
    for (int i = 0; i < CONST_SIZE; i++) {
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

    std::cout << "building tree" << std::endl;
    rootNode = buildTree(constellationLarsson, 0, false, false);
    rootNode->index = 0;
    std::cout << "numbering nodes" << std::endl;
    numberTree(rootNode);

    generateLookup(constellationJellema);

    //std::cout << "printing nodes" << std::endl;
    //printTree(rootNode);


    std::cout << "simulating" << std::endl;

    const double mean = 0.0;
    const double stddev = 0.2;
    std::mt19937_64 generator(std::random_device{}());
    std::normal_distribution<double> distributionNormal(mean, stddev);
    std::uniform_int_distribution<int> distributionUniform(0, CONST_SIZE - 1);

    //Simulate a ton of random data
    int errorCountSymbol = 0;
    int errorCountBitLarsson = 0;
    int errorCountBitJellema = 0;
    for (int64_t i = 0; i < SIMULATION_SIZE ; i++) {
        int nodeIndex = distributionUniform(generator);

        treeNode p = constellation[nodeIndex];

        double receivedI = p.I + distributionNormal(generator);
        double receivedQ = p.Q + distributionNormal(generator);

        int closests = nearestNeightbour(constellation, receivedI, receivedQ);
        
        if (closests != nodeIndex) {
            errorCountSymbol++;
            errorCountBitLarsson += hammingDist(closests, nodeIndex);
            errorCountBitJellema += hammingDist(lookupTableJellema[closests], lookupTableJellema[nodeIndex]);

        }
    }

    std::cout << "Symbol errors: " << errorCountSymbol << " "<< (errorCountSymbol * 100) / SIMULATION_SIZE << "%" << std::endl;
    std::cout << "Bit error Larsson: " << errorCountBitLarsson << " " << (errorCountBitLarsson * 100.0) / (SIMULATION_SIZE * CONST_SIZE_POWER) << "%" << std::endl;
    std::cout << "Bit error Jellema: " << errorCountBitJellema << " " << (errorCountBitJellema * 100.0) / (SIMULATION_SIZE * CONST_SIZE_POWER) << "%" << std::endl;

    return 0;
}