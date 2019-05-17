#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <random>


const unsigned int CONST_SIZE_POWER = 8;
const unsigned int CONST_SIZE = 1 << CONST_SIZE_POWER;
const unsigned int SIMULATION_SIZE = 1'000'000;

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
treeNode constellationCopy[CONST_SIZE];

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


void createMappingSquareAA(std::vector<std::vector<bool>> splits) {
    for (int i = 0; i < CONST_SIZE_POWER; i++) {
        int partitions = 1 << i;
        for (int j = 0; j < partitions; j++) {
            int size_div = CONST_SIZE / partitions;
            int start_index = j*size_div;
            int end_index   = (j+1)*size_div;


            if (splits[i][j] == 1) {
                if (i%2) {
                    std::sort(constellation + start_index, constellation + end_index, &compareQ);
                } else {
                    std::sort(constellation + start_index, constellation + end_index, &compareI);
                }
            } else {
                if (i%2) {
                    std::sort(constellation + start_index, constellation + end_index, &compareQreverse);
                } else {
                    std::sort(constellation + start_index, constellation + end_index, &compareIreverse);
                }
            }
        }
    }
}

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
    std::cout << "Node index: " << node->index << std::endl;

    if (node->leave[1]) {
        node->leave[1]->index = ((node->index << 1) | 1);
        numberTree(node->leave[1]);
    }
}

/*
 * Perform a nearestNeighbourhood search using the constellation tree that was build.
 * TODO: implement proper nearest neighbour search
 */
int nearestNeightbour(coord_t I, coord_t Q) {
    double currentBestDistance = 10;
    int indexBest = 0;

    for (int i = 0; i < CONST_SIZE; i++) {
        double diffI = constellation[i].I - I;
        double diffQ = constellation[i].Q - Q;
        double distance =  diffI * diffI + diffQ * diffQ;
        if (distance < currentBestDistance) {
            indexBest = i;
            currentBestDistance = distance;
        }
    }

    return indexBest;
};

int visited;


//TODO: this function is broken as fuck due to not giving I and Q values to the tree intermediates
void nearest(struct treeNode *root, struct treeNode *nd, int i,
             struct treeNode **best, double *best_dist) {
    double d, dx, dx2;

    if (!root) {
        return;
    }
    d = dist(rootNode, nd);
    if (i == 0) {
        dx = root->I - nd->I;
    } else {
        dx = root->Q - nd->Q;
    }
    dx2 = dx * dx;

    visited++;

    if (!*best || d < *best_dist) {
        *best_dist = d;
        *best = root;
    }

    i ^= 0x01;

    nearest(dx > 0 ? root->leave[0] : root->leave[1], nd, i, best, best_dist);

    if (dx2 >= *best_dist) {
        return;
    }

    nearest(dx > 0 ? root->leave[1] : root->leave[0], nd, i, best, best_dist);
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


std::vector<std::vector<bool>> splitGenerator() {
    std::vector<std::vector<bool>> splits;

    for (int i = 0; i < CONST_SIZE_POWER; i++) {
        splits.emplace_back(std::vector<bool>());
    }

    splits[0].push_back(true);
    splits[1].push_back(true);
    splits[1].push_back(true);

    for (int i = 2; i < CONST_SIZE_POWER; i++) {
        int j = 0;
        while (j < splits[i-2].size()) {
            if (splits[i-2][j] == 1) {
                splits[i].push_back(false);
                splits[i].push_back(false);
                splits[i].push_back(true);
                splits[i].push_back(true);
            } else {
                splits[i].push_back(true);
                splits[i].push_back(true);
                splits[i].push_back(false);
                splits[i].push_back(false);
            }
            j = j + 1;
        }

    }
    return splits;
}

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
    std::copy(std::begin(constellation), std::end(constellation), std::begin(constellationCopy));

    createMappingSquareAA(splitGenerator());

    //std::copy(std::begin(constellation), std::end(constellation), std::begin(constellationCopy));

    std::cout << "building tree" << std::endl;
    rootNode = buildTree(constellationCopy, 0, false, false);
    rootNode->index = 0;
    std::cout << "numbering nodes" << std::endl;
    numberTree(rootNode);
    std::cout << "printing nodes" << std::endl;
    printTree(rootNode);


    std::cout << "simulating" << std::endl;

    const double mean = 0.0;
    const double stddev = 0.2;
    std::default_random_engine generator(std::random_device{}());
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

        int closests = nearestNeightbour(receivedI, receivedQ);

        //TODO: nearest neighbour using KD tree is broken see function description.
        /*
        treeNode test;
        test.I = receivedI;
        test.Q = receivedQ;
        double best_dist;
        struct treeNode *found;
        found = nullptr;

        nearest(rootNode, &test, 0, &found, &best_dist);
        int closestsFast = found->index;

        std::cout << found->I << " " << constellation[closests].I << std::endl;
        std::cout << found->Q << " " << constellation[closests].Q << std::endl;
        */

        if (closests != nodeIndex) {
            errorCountSymbol++;
            errorCountBitJellema += hammingDist(closests, nodeIndex);
            errorCountBitLarsson += hammingDist(constellation[closests].index, constellation[nodeIndex].index);
        }
    }

    std::cout << "Symbol errors: " << errorCountSymbol << " "<< (errorCountSymbol * 100) / SIMULATION_SIZE << "%" << std::endl;
    std::cout << "Bit error Larsson: " << errorCountBitLarsson << " " << (errorCountBitLarsson * 100.0) / (SIMULATION_SIZE * CONST_SIZE_POWER) << "%" << std::endl;
    std::cout << "Bit error Jellema: " << errorCountBitJellema << " " << (errorCountBitJellema * 100.0) / (SIMULATION_SIZE * CONST_SIZE_POWER) << "%" << std::endl;

    return 0;
}