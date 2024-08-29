//
// Created by nboardma on 7/18/2019.
//

#ifndef CPPFILES_NETWORKELEMENTS_H
#define CPPFILES_NETWORKELEMENTS_H
#include <vector>
#include <random>
using namespace std;

class Node {
// myID: refers to the node number, 1, 2, ...
// myValue: a value associated with the node. Ex: fail time, value, etc,

//TODO: maybe add in an node identifier (sink, sensor, target). If target add a weight value

public:
    Node();
    Node(int inputID);

    void SetID(int temp);
    int GetNodeID();

    int GetValue();
    void SetValue(int data);

    void UpdateCapacityLabel(int newCapacity);
    int GetMaxCapacityPath();

    Node *GetPrev();
    Node *GetNext();
    Node *GetFirstElement();

    void SetPrev(Node *newPrev);
    void SetNext(Node *newNext);
    void SetFirstItem(Node *firstElement);

    void UpdateReached();
    bool IsReached();
    void MarkPermanent();
    bool IsPermanent();

    void ResetNodeReached();
    void ResetNodePermanent();

    void AddReachableNode(int nodeR);
    int GetReachableNode(int i);
    int GetDegree();
    void ResetReachableNodes();

    int GetNumberInCell();
    void AddNodeToBucket(int nodeNew);
    int GetNodeInBucket(int i);
    void ResetNodesInBucket();

    void SetNewLoc(double xLoc, double yLoc);
    double GetXLoc();
    double GetYLoc();

    bool IsCovered();
    void SetCovered();
    void ResetCovered();

    bool IsFailed();
    void SetFailed();
    void ResetFailed();

protected:
    vector<int> forwardStar;
    int nodeDegree = 0;
    bool isReached = false;
    bool isCovered = false;

    bool isFailed = false;

    int myValue;
    int myID;
    int myMaxCapacityPath = 0;

    double myXLoc;
    double myYLoc;

    vector<int> nodesInBucket;
    int myNodesInBucket = 0;
    bool isPermanent = false;

private:
    Node *prev;
    Node *next;
    Node *firstInBucket;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class Edge {
// TODO: add an identifier for a directed edge or an undirected edge
public:
    Edge(int start, int end): myStartNode(start), myEndNode(end), myWeight(0){};
    Edge(int start, int end, double weight): myStartNode(start), myEndNode(end), myWeight(weight){};

public:
    int myStartNode;
    int myEndNode;
    double myWeight;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class NodeLinkedList {
    // numItems - number of nodes in the network
    // vector<int> time - a vector indicating order of failure

public:
    NodeLinkedList();
    NodeLinkedList(int numItems);
    NodeLinkedList(int numItems, vector<int> order);
    NodeLinkedList(int numItems, int numSens, int numTarg);
    NodeLinkedList(int numItems, int numSens, int numTarg, vector<int> order);

    void InitializeBuckets(int numItems);
    void InitializeBuckets(int numItems, vector<int> failOrder);

    void AddToList(int listNum, Node *newElement);
    Node *RemoveFirst(int listNum);
    Node *RemoveTargetElement(int listNum, Node *removeElement);
    bool IsListEmpty(int listNum);

    void BuildForwardStar(vector<Edge*> edgeList);

    void DialsAlgorithm(vector<Edge*> edgeList, int numNodes);
    void DialsAlgorithm(vector<Edge*> edgeList, int numNodes, int numTar, double &timeDuration);
    void DijkstrasAlgorithm(vector<Edge*> edgeList, int numNodes);
    void DijkstrasAlgorithm(vector<Edge*> edgeList, int numNodes, int numTar, double &timeDuration);
    void DialsAlgorithmUpdateRule(vector<Edge *> edgeList, int numNodes, int numTar, double &timeDuration);

    void UpdateDSpecFeq(int numFailed);
    void CalculateDSpec(int numReps);
    void OutputDSpec();

    void EfficientFrontier(double maxDelta, double costFixed, double costVar, double myShape, double myScale);

    void SignatureRelation(int newSize, double t, double maxDelta, double costFixed, double costVar, double myShape, double myScale);
    void SignatureRelationOutput(int sizeRange, double t, double maxDelta, double costFixed, double costVar, double myShape, double myScale);
    void SignatureRelationVariance(int sizeRange, int numReps, double t, double maxDelta, double costFixed, double costVar, double myShape, double myScale);

    void CalculateReliabilityWeibull(double t, double myShape, double myScale);
    void CalculateReliabilityVarianceWeibull(double t, int numReps, double myShape, double myScale);

    void CalculateStableMaintenanceReliabilityWeibull(double delta, double mShape, double mScale);
    void CalculateStableMaintenanceVarianceWeibull(double delta, int numReps , double mShape, double mScale);

    void CalculateMaintenancePolicyCostWeibull(double delta, double costFixed, double costVar, double mShape, double mScale);
    void CalculateMaintenancePolicyCostWeibull(int numSens, double delta, double costFixed, double costVar, double mShape, double mScale);

    void GenerateNewFailOrder(vector<int> &nodeTime);

    void GenerateNewRGG_naive();
    void GenerateNewRGG_Bucket();
    void CompareAdjacentCells(int cell1, int cell2);
    void UniformTargetLoc();

    void OutputCriticalTimes();
    void addFailTime(int failTime);
    void addFailOrder(int nodeNum);
    int GetFailTime(int num);
    int GetFailOrder(int num);

    void OutputTargetCriticalTimes(int numTargets);
    void OutputTargetFailOrder(int numTargets);
    void addTargetFailTime(int failTime);
    void addTargetFailOrder(int nodeNum);
    int GetSensorFailOrder(int senNum);
    int GetTargetFailTime(int failNum);
    int GetTargetFailOrder(int targNum);

    void OutputSensorFailOrder(int numSensors);
    void addSensorFailTime(int failTime);
    void addSensorFailOrder(int nodeNum);

    void SetRanges(double comRange, double sensingRange);

    int CoverageNumber(double coverage, int numT);

    void UpdateBucketTimes(vector<int> newOrder);
    void ClearOrderTimes();

protected:
    Node *firstItem;
    Node *lastItem;
    Node *myBuckets;
    int myNumElements;  // total number of nodes (includes sink, sensors, targets)
    int myNumSensors;   // number of sensor nodes, also includes the sink node
    int myNumTargets;   // number of targets

    random_device myDevX;       // used to generate x locations for a new RGG
    random_device myDevY;       // used to generate y locations for a new RGG
    random_device myDevOrder;   // used to generate a new order of sensor failures

    double mySensorRange;
    double myCommunicationRange;

    vector<int> orderedCriticalFailTime;
    vector<int> orderedFail;

    vector<int> orderedTargetCriticalFailTime;
    vector<int> orderedTargetFail;

    vector<int> orderedSensorCriticalFailTime;
    vector<int> orderedSensorFail;

    vector<int> DSpecFreq;
    vector<int> DSpecElements;

    vector<int> myDSpec;
    vector<double> myDSpecProb;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////


class SimulationNodeMC : public Node {

public:
    SimulationNodeMC();

    int GetNodeID();
    double GetFailTime();
    void SetFailTime(double t);

    void UpdateFailCapacity (double newTime);
    double GetFailCapacity();

    void SetNewLoc(double xLoc, double yLoc);
    void ResetNode();

    double GetXLoc();
    double GetYLoc();

    bool IsCovered();
    void SetCovered();
    void ResetCovered();

    bool IsFailed();
    void SetFailed();
    void ResetFailed();

    void AddNodeToBucket(int nodeNew);
    void ResetReachableNodes();
    void ResetNodesInBucket();

    void UpdateForwardStar(int removeNum);

    int GetNodeInBucket(int i);
    int GetNumberInCell();

private:
    double myFailTime;
    double myCapacityFailTime = 0;

};

class NodeLinkedListMC : public NodeLinkedList {

public:
    NodeLinkedListMC(int numItems, int numSens, int numTarg);
    NodeLinkedListMC(int numItems, int numSens, int numTarg, vector<double> failTimes);

    void InitializeBuckets(int numItems);
    void InitializeBuckets(int numItems, vector<double> failTimes);
    void DijkstrasAlgorithm(vector<Edge*> edgeList, int numNodes, int numTar, double &timeDuration);

    void addFailTime(double failTime);
    void addTargetFailTime(double failTime);
    void addSensorFailTime(double failTime);
    void addNetworkFailTime(double failTime);

    double GetTargetFailTime(int failNum);
    void UniformTargetLoc();
    void GenerateNewRGG_Bucket();
    void CompareAdjacentCells(int cell1, int cell2);

    void GenerateNewFailTime(vector<double> &nodeTime);
    void UpdateFailTimes(vector<double> newTimes);
//    void UpdateFailTimes(vector<double> newTimes, int epochNum, double delta);

    void CalculateReliability(double t, int numReps);

    double TargetCoveragePercentage();

    void DijkstrasAlgorithmMainten(int numNodes, int numTar, double &timeDuration, double delta, vector<int> &failCount);

    void BreadthFirstSearch(int numNodes, int numTar, double &timeDuration);
    double BreadthFirstSearchReturn(int numNodes, int numTar, double &timeDuration);

    void GenerateNewFailTimeMaintenance(vector<double> &nodeTime, double t, double curTime);
    void GenerateNewRGG_BucketMaint();
    void GenerateNewRGG_BucketMaint(double costFixed, double costVar, double &policyCost);
    void UpdateRGG(double time);

protected:
    SimulationNodeMC *myBuckets;

    vector<double> orderedCriticalFailTime;
    vector<double> orderedTargetCriticalFailTime;
    vector<double> orderedSensorCriticalFailTime;

    vector<double> networkFailTimes;

    random_device myDevTime;


};




/*
class BinaryHeap {

public:
    BinaryHeap();
    BinaryHeap(int numItems, int sourceNode);

    // heap operations
    void CreateHeap(int numItems, int sourceNode);
    int FindMin();
    void DeleteMin();
    void DecreaseKey(int nodeNum, double newVal);
    void IncreaseKey(int nodeNum, double newVal);
    void SiftUp(int nodeNum);
    void SiftDown(int nodeNum);
    void InsertElement(int nodeNum, double val);

protected:
    SimulationNodeMC *myHeap;
    int myNumElements = 0;
    vector<int> heapOrder;

};

class Dijkstras {

public:
    Dijkstras(int numNodes);
    Dijkstras(int numNodes, int sourceNode);

    void DijkstrasHeap(vector<Edge*> edgeList);

    void DijkstrasNewSourceHeap(int source);

    void DijkstrasMultipleHeap(vector<Edge*> edgeList, vector<int> sources);

    void GetPath(int endNode);

    void InitializeHeap(int numItems, int sourceNode);

    void ResetDijkstras(int newSource);
    void SetSource(int num);

protected:
    BinaryHeap dijkstrasHeap;
    int myNumNodes;
    int mySourceNode;

};*/

// reads in a txt file where the format of the first row is:
// number of nodes(sensors)          number of arcs              number of targets
void ReadNetwork(char* filename, vector<Edge*> &edgeVect, vector<int> &nodeVect, int &numTargets);

// reads in a txt file where the format of the first row is:
// number of nodes          number of arcs
void ReadNetwork(char* filename, vector<Edge*> &edgeVect, vector<int> &nodeVect);

// used to update the fail times for a new test instance. Assumes the remaining network topology is the same
void UpdateTimes(char* filename, vector<int> &nodeVect, int &numIterations);

#endif //CPPFILES_NETWORKELEMENTS_H
