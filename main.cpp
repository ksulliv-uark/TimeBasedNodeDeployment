#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include "NetworkElements.h"
#include "Reliability.h"

using namespace std;

// Run a single instance with a single input file
void RunAlgorithms(char *fileName){

    cout << " New Test " << endl;

    vector<Edge*> networkEdges;
    vector<int> nodeTime;
    int numTargets;

    ReadNetwork(fileName, networkEdges, nodeTime, numTargets);

    // number of nodes includes sink node, number of sensors, and number of targets
    int numNodes = nodeTime.size();

    double timeDials;
    double timeDijkstras;

    /*
    cout << "Num Nodes : " << numNodes << endl;
    cout << "Num Sensors : " << numNodes - numTargets << endl;
    cout << "Num Targets : " << numTargets << endl;
    */

    NodeLinkedList myDialsBuckets(numNodes, (numNodes - numTargets), numTargets, nodeTime);
    NodeLinkedList myDijkstrasBuckets(numNodes, (numNodes - numTargets), numTargets, nodeTime);

    myDialsBuckets.BuildForwardStar(networkEdges);
    myDijkstrasBuckets.BuildForwardStar(networkEdges);

    myDialsBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);
    myDijkstrasBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);

    cout << " Time Dials : " << timeDials << endl;
    cout << " Time Dijkstras : " << timeDijkstras << endl;

    myDialsBuckets.OutputCriticalTimes();
    cout << " -------------------- " << endl;
    myDialsBuckets.OutputTargetCriticalTimes(numTargets);
    cout << " -------------------- " << endl;
    myDialsBuckets.OutputTargetFailOrder(numTargets);
//    myDialsBuckets.OutputSensorFailOrder(numNodes - numTargets - 1);

//    myDijkstrasBuckets.OutputCriticalTimes();

//    int critNum = myDialsBuckets.CoverageNumber(0.8, numTargets);

}

// used to run when there is a single update file - update file contains new fail times
void UpdateInstance(char *fileName, char *updateFile){

    vector<Edge*> networkEdges;
    vector<int> nodeTime;
    int numTargets;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    ReadNetwork(fileName, networkEdges, nodeTime, numTargets);

    // number of nodes includes sink node, number of sensors, and number of targets
    int numNodes = nodeTime.size();

    double timeDials;
    double timeDijkstras;

    NodeLinkedList myNetworkBuckets(numNodes, nodeTime);

    myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);
//    myNetworkBuckets.DijkstrasAlgorithm(networkEdges,numNodes,numTargets, timeDijkstras);

    myNetworkBuckets.OutputTargetCriticalTimes(numTargets);

    // Start of updating node failure times
    int tempLoc = 0;
    int numIterations;
    vector<int> updateTimes; // stores all update times in a single vector

    // read in new times from an input file
    UpdateTimes(updateFile, updateTimes, numIterations);

    for (long j = 1; j <= numIterations ; ++j) {
        // clear out current fail times
        nodeTime.clear();
        // read in new times from those stored in vector
        for (int i = 0; i < numNodes; ++i){
            nodeTime.push_back(updateTimes.at(tempLoc));
            tempLoc = tempLoc + 1;
        }

        myNetworkBuckets.UpdateBucketTimes(nodeTime);
        myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);

//        myNetworkBuckets.UpdateBucketTimes(nodeTime);
//        myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);

        myNetworkBuckets.OutputTargetCriticalTimes(numTargets);

    }

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finishN - startN;

    cout << "Time To Run : " << elapsed.count() << endl;

}

// used to run when there are multiple update files - update file contains new fail times
void MultipleInstancesDials(int numFiles, char *fileNames[]){

    double coverageLevel = 0.8;

    vector<Edge*> networkEdges;
    vector<int> nodeTime;
    int numTargets;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    ReadNetwork(fileNames[1], networkEdges, nodeTime, numTargets);

    // number of nodes includes sink node, number of sensors, and number of targets
    int numNodes = nodeTime.size();
    double timeDials;

    NodeLinkedList myNetworkBuckets(numNodes, nodeTime);
    myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);

    myNetworkBuckets.OutputTargetCriticalTimes(numTargets);
//    int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
//    cout << "Coverage drops below " << coverageLevel*100 << "% with sensor failure : " << myNetworkBuckets.GetTargetFailTime(critNum) << endl;

    // Start of updating node failure times
    for (int k = 2; k < numFiles ; ++k) {

        int tempLoc = 0;
        int numIterations;
        vector<int> updateTimes; // stores all update times in a single vector

        // read in new times from an input file
        UpdateTimes(fileNames[k], updateTimes, numIterations);

        for (long j = 1; j <= numIterations ; ++j) {
            // clear out current fail times
            nodeTime.clear();
            // read in new times from those stored in vector
            for (int i = 0; i < numNodes; ++i){
                nodeTime.push_back(updateTimes.at(tempLoc));
                tempLoc = tempLoc + 1;
            }

            myNetworkBuckets.UpdateBucketTimes(nodeTime);
            myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);

            myNetworkBuckets.OutputTargetCriticalTimes(numTargets);
//            int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
//            cout << "Coverage drops below " << coverageLevel*100 << "% with sensor failure : " << myNetworkBuckets.GetTargetFailTime(critNum) << endl;

        }

    }

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finishN - startN;

    cout << "TimeDials " << elapsed.count() << endl;

}

// used to run when there are multiple update files - update file contains new fail times
void MultipleInstancesDijkstras(int numFiles, char *fileNames[]){

    double coverageLevel = 0.7;

    vector<Edge*> networkEdges;
    vector<int> nodeTime;
    int numTargets;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    ReadNetwork(fileNames[1], networkEdges, nodeTime, numTargets);

    // number of nodes includes sink node, number of sensors, and number of targets
    int numNodes = nodeTime.size();
    double timeDijkstras;

    NodeLinkedList myNetworkBuckets(numNodes, nodeTime);
    myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);

    myNetworkBuckets.OutputTargetCriticalTimes(numTargets);
//    int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
//    cout << "Coverage drops below " << coverageLevel*100 << "% with sensor failure : " << myNetworkBuckets.GetTargetFailTime(critNum) << endl;

    // Start of updating node failure times
    for (int k = 2; k < numFiles ; ++k) {
//        cout << "File Name : " << fileNames[k] << endl;

        int tempLoc = 0;
        int numIterations;
        vector<int> updateTimes; // stores all update times in a single vector

        // read in new times from an input file
        UpdateTimes(fileNames[k], updateTimes, numIterations);

        for (long j = 1; j <= numIterations ; ++j) {
            // clear out current fail times
            nodeTime.clear();
            // read in new times from those stored in vector
            for (int i = 0; i < numNodes; ++i){
                nodeTime.push_back(updateTimes.at(tempLoc));
                tempLoc = tempLoc + 1;
            }

            myNetworkBuckets.UpdateBucketTimes(nodeTime);
            myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);

            myNetworkBuckets.OutputTargetCriticalTimes(numTargets);
//            int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
//            cout << "Coverage drops below " << coverageLevel*100 << "% with sensor failure : " << myNetworkBuckets.GetTargetFailTime(critNum) << endl;

        }

    }

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finishN - startN;

    cout << "TimeDijkstra " << elapsed.count() << endl;

}

// used to generate new fail times
void NewFailOrderInstance(int numSensors, int numTargets, vector<int> &nodeTime){

// TODO change default random engine to a different generator
//    default_random_engine myGen;
//    mersenne_twister_engine myMersenne();

    random_device myDev;

    nodeTime.clear();
    nodeTime.reserve(numSensors + numTargets);

    for (int j = 0; j < numSensors ; ++j) {
        nodeTime.push_back(j + 1);
    }

    for (int j = 0; j < numSensors ; ++j) {
        // generate random number
        uniform_int_distribution<int> myDist(0, numSensors - j - 1);
        int myNum = myDist(myDev);

        // store the last number that is not marked
        int tempFirst = nodeTime.at(myNum);
        // element to swap
        int tempLast = nodeTime.at(numSensors - j - 1);

        // swap elements
        nodeTime[myNum] = tempLast;
        nodeTime[numSensors - j - 1] = tempFirst;

    }

    // set fail time for sink
    nodeTime.insert(nodeTime.begin(), numSensors + 1);

    // initial time for all targets
    for (int j = 0; j < numTargets ; ++j) {
        nodeTime.push_back(numSensors + 1);
    }
}


// generate a new test instance within c++
void GenerateInstancesDials(char *fileName, int numReps, double time, double coverageLevel, double cFixed, double cVar){

    vector<Edge*> networkEdges;
    vector<int> nodeTime;
    int numTargets;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    ReadNetwork(fileName, networkEdges, nodeTime, numTargets);

    // number of nodes includes sink node, number of sensors, and number of targets
    int numNodes = nodeTime.size();
    double timeDials;

    NodeLinkedList myNetworkBuckets(numNodes, (numNodes - numTargets), numTargets, nodeTime);

    myNetworkBuckets.BuildForwardStar(networkEdges);
    myNetworkBuckets.SetRanges(0.0750000, 0.0750000);
    myNetworkBuckets.UniformTargetLoc();

    myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);

 //   myNetworkBuckets.OutputTargetCriticalTimes(numTargets);

    int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
    myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));


    // Start of updating node failure times
    for (int k = 1; k < numReps ; ++k) {

//        NewFailOrderInstance(numNodes - numTargets - 1, numTargets, nodeTime);

        myNetworkBuckets.GenerateNewFailOrder(nodeTime);

//        myNetworkBuckets.GenerateNewRGG_naive();
//        cout << " ---------------- New RGG ----------------  " << endl;
        myNetworkBuckets.GenerateNewRGG_Bucket();

        myNetworkBuckets.UpdateBucketTimes(nodeTime);
        myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);

//        myNetworkBuckets.OutputTargetCriticalTimes(numTargets);
//        cout << myNetworkBuckets.GetTargetFailTime(critNum) << endl;


        myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));

    }

    myNetworkBuckets.CalculateDSpec(numReps);
    myNetworkBuckets.CalculateReliabilityWeibull(time, 1.5, 10.0);

//    myNetworkBuckets.CalculateReliabilityVarianceWeibull(time, numReps+1, 1.5, 10.0);

    double i = 1.0;
    double curShape = 1.5;
    double curScale = 10.0;

    ///////////////////////////////////////////////

    cout << " Shape = " << curShape << " , Scale = " << curScale << " F Cost = " << cFixed << " v cost = " << cVar << endl;
    cout << " Stable Maintenance Reliability for changing delta " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateStableMaintenanceReliabilityWeibull(i, curShape, curScale);
        i = i + 0.1;

    }

    i = 1.0;
    cout << "Stable maintenance cost " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateMaintenancePolicyCostWeibull(i, cFixed, cVar, curShape, curScale);
        i = i + 0.1;
    }


//    myNetworkBuckets.OutputDSpec();

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finishN - startN;
    cout << numReps << " - " << (numNodes-441) <<  " TimeDials " << elapsed.count() << endl;

}

// generate a new test instance within c++
void GenerateInstancesDijkstras(char *fileName, int numReps, double time, double coverageLevel, double cFixed, double cVar){

    vector<Edge*> networkEdges;
    vector<int> nodeTime;
    int numTargets;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    ReadNetwork(fileName, networkEdges, nodeTime, numTargets);

    // number of nodes includes sink node, number of sensors, and number of targets
    int numNodes = nodeTime.size();
    double timeDijkstras;

    NodeLinkedList myNetworkBuckets(numNodes, (numNodes - numTargets), numTargets, nodeTime);

    myNetworkBuckets.BuildForwardStar(networkEdges);
    myNetworkBuckets.SetRanges(0.0750000, 0.0750000);
    myNetworkBuckets.UniformTargetLoc();

    myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);

//    myNetworkBuckets.OutputTargetCriticalTimes(numTargets);
    int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
    myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));

    // Start of updating node failure times
    for (int k = 1; k < numReps ; ++k) {

//        NewFailOrderInstance(numNodes - numTargets - 1, numTargets, nodeTime);

        myNetworkBuckets.GenerateNewFailOrder(nodeTime);

//        myNetworkBuckets.GenerateNewRGG_naive();
        myNetworkBuckets.GenerateNewRGG_Bucket();

        myNetworkBuckets.UpdateBucketTimes(nodeTime);
        myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);

//        myNetworkBuckets.OutputTargetCriticalTimes(numTargets);
        myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));

    }

    myNetworkBuckets.CalculateDSpec(numReps);
    myNetworkBuckets.CalculateReliabilityWeibull(time, 1.5, 10.0);
//    myNetworkBuckets.CalculateReliabilityVarianceWeibull(time, numReps + 1, 1.5, 10.0);


    double i = 1.0;
    double curShape = 1.5;
    double curScale = 10.0;

    cout << " Shape = " << curShape << " , Scale = " << curScale << " F Cost = " << cFixed << " v cost = " << cVar << endl;
    cout << " Stable Maintenance Reliability for changing delta " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateStableMaintenanceReliabilityWeibull(i, curShape, curScale);
        i = i + 0.1;
    }

    i = 1.0;
    cout << "Stable maintenance cost " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateMaintenancePolicyCostWeibull(i, cFixed, cVar, curShape, curScale);
        i = i + 0.1;
    }

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finishN - startN;

    cout << numReps << " - " << (numNodes-441) << " TimeDijkstras " << elapsed.count() << endl;
}

// creates a random network to initialize, no input file
void DialsRandomNetwork(int nSens, int numReps, double time, double coverageLevel, double cFixed, double cVar){

    vector<Edge*> networkEdges;
    vector<int> nodeTime;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    // number of nodes includes sink node, number of sensors, and number of targets
    int numTargets = 441;
    int numNodes = nSens + numTargets + 1;
    double timeDials;
    double totalTimeDials = 0;
    double networkGenerationTime = 0;
    double simTimeTotal = 0;

    NodeLinkedList myNetworkBuckets(numNodes, (numNodes - numTargets), numTargets);
    myNetworkBuckets.SetRanges(0.0750000, 0.0750000);
    myNetworkBuckets.UniformTargetLoc();

    auto startSimTime = chrono::high_resolution_clock::now();

    myNetworkBuckets.GenerateNewFailOrder(nodeTime);

    auto endSimeTime = chrono::high_resolution_clock::now();
    chrono::duration<double> simTime = endSimeTime - startSimTime;

    simTimeTotal = simTimeTotal + simTime.count();


    auto startNetwork = chrono::high_resolution_clock::now();

//    myNetworkBuckets.GenerateNewRGG_naive();
    myNetworkBuckets.GenerateNewRGG_Bucket();

    auto finishNetwork = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedNetwork = finishNetwork - startNetwork;
    networkGenerationTime = networkGenerationTime + elapsedNetwork.count();

    myNetworkBuckets.UpdateBucketTimes(nodeTime);

    myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);
    totalTimeDials = totalTimeDials + timeDials;

    int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
    myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));

    // Start of updating node failure times
    for (int k = 1; k < numReps ; ++k) {

        startSimTime = chrono::high_resolution_clock::now();

        myNetworkBuckets.GenerateNewFailOrder(nodeTime);

        endSimeTime = chrono::high_resolution_clock::now();
        simTime = endSimeTime - startSimTime;
        simTimeTotal = simTimeTotal + simTime.count();

        auto startNetworkLoop = chrono::high_resolution_clock::now();

//        myNetworkBuckets.GenerateNewRGG_naive();
        myNetworkBuckets.GenerateNewRGG_Bucket();

        finishNetwork = chrono::high_resolution_clock::now();
        elapsedNetwork = finishNetwork - startNetworkLoop;
        networkGenerationTime = networkGenerationTime + elapsedNetwork.count();

        myNetworkBuckets.UpdateBucketTimes(nodeTime);

        myNetworkBuckets.DialsAlgorithm(networkEdges, numNodes, numTargets, timeDials);
        totalTimeDials = totalTimeDials + timeDials;
        myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));
    }

    auto finishLoop = chrono::high_resolution_clock::now();

    myNetworkBuckets.CalculateDSpec(numReps);
    myNetworkBuckets.CalculateReliabilityWeibull(time, 1.5, 10.0);

//    myNetworkBuckets.CalculateReliabilityVarianceWeibull(time, numReps, 1.5, 10.0);

    double i = 1.0;
    double curShape = 1.5;
    double curScale = 10.0;
    double maxDelta = 10.0;

    ///////////////////////////////////////////////

    auto startPolicy = chrono::high_resolution_clock::now();
//    cout << " Shape = " << curShape << " , Scale = " << curScale << " F Cost = " << cFixed << " v cost = " << cVar << endl;
    cout << " Stable Maintenance Reliability for changing delta " << endl;
    while (i < maxDelta){
        myNetworkBuckets.CalculateStableMaintenanceReliabilityWeibull(i, curShape, curScale);
        i = i + 0.1;

    }

    i = 1.0;
    cout << " Stable Maintenance Variance for changing delta " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateStableMaintenanceVarianceWeibull(i, numReps, curShape, curScale);
        i = i + 0.5;
    }


    i = 1.0;
    cout << "Stable maintenance cost " << endl;
    while (i < maxDelta){
        myNetworkBuckets.CalculateMaintenancePolicyCostWeibull(i, cFixed, cVar, curShape, curScale);
        i = i + 0.1;
    }

    myNetworkBuckets.OutputDSpec();

    ////////////////////////////////////////////////

//        myNetworkBuckets.SignatureRelationOutput(nSens, time, maxDelta, cFixed, cVar, curShape, curScale);

    auto endPolicy = chrono::high_resolution_clock::now();
    chrono::duration<double> timeEvaluatingMaintPolicies = endPolicy - startPolicy;

//    myNetworkBuckets.SignatureRelationVariance(nSens, numReps, time, maxDelta, cFixed, cVar, curShape, curScale);

    auto endPolicyVar = chrono::high_resolution_clock::now();
    chrono::duration<double> timeEvaluatingMaintVar = endPolicyVar - endPolicy;

//    myNetworkBuckets.SignatureRelation(851, time, maxDelta, cFixed, cVar, curShape, curScale);

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finishN - startN;

    chrono::duration<double> dTime = finishLoop - startNetwork;

    cout << " Time Spent Generating Network : " << networkGenerationTime << endl;
    cout << " Time Spent implementing Dials : " << totalTimeDials << endl;
    cout << " Time Generating Fail Time : " << simTimeTotal << endl;
    cout << " Total Time for destruction algorithm : " << dTime.count() << endl;
    cout << " Time evaluating TBM policies : " << timeEvaluatingMaintPolicies.count() << endl;
    cout << " Time estimating Variance : " << timeEvaluatingMaintVar.count() << endl;
    cout << numReps << " - " << " TotalTimeDials " << elapsed.count() << endl;

}

void DijkstrasRandomNetwork(int nSens, int numReps, double time, double coverageLevel, double cFixed, double cVar){

    vector<Edge*> networkEdges;
    vector<int> nodeTime;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    // number of nodes includes sink node, number of sensors, and number of targets
    int numTargets = 441;
    int numNodes = nSens + numTargets + 1;
    double timeDijkstras;
    double totalTimeDijkstras = 0;
    double networkGenerationTime = 0;
    double simTimeTotal = 0;

    NodeLinkedList myNetworkBuckets(numNodes, (numNodes - numTargets), numTargets);
    myNetworkBuckets.SetRanges(0.0750000, 0.0750000);
    myNetworkBuckets.UniformTargetLoc();

    auto startSimTime = chrono::high_resolution_clock::now();

    myNetworkBuckets.GenerateNewFailOrder(nodeTime);

    auto endSimeTime = chrono::high_resolution_clock::now();
    chrono::duration<double> simTime = endSimeTime - startSimTime;

    simTimeTotal = simTimeTotal + simTime.count();

    auto startNetwork = chrono::high_resolution_clock::now();

//    myNetworkBuckets.GenerateNewRGG_naive();
        myNetworkBuckets.GenerateNewRGG_Bucket();

    auto finishNetwork = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedNetwork = finishNetwork - startNetwork;
    networkGenerationTime = networkGenerationTime + elapsedNetwork.count();

    myNetworkBuckets.UpdateBucketTimes(nodeTime);

    myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);
    totalTimeDijkstras = totalTimeDijkstras + timeDijkstras;

    int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
    myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));


    // Start of updating node failure times
    for (int k = 1; k < numReps ; ++k) {

        startSimTime = chrono::high_resolution_clock::now();

        myNetworkBuckets.GenerateNewFailOrder(nodeTime);

        endSimeTime = chrono::high_resolution_clock::now();
        simTime = endSimeTime - startSimTime;
        simTimeTotal = simTimeTotal + simTime.count();

        auto startNetworkLoop = chrono::high_resolution_clock::now();

//        myNetworkBuckets.GenerateNewRGG_naive();
            myNetworkBuckets.GenerateNewRGG_Bucket();

        finishNetwork = chrono::high_resolution_clock::now();
        elapsedNetwork = finishNetwork - startNetworkLoop;
        networkGenerationTime = networkGenerationTime + elapsedNetwork.count();

        myNetworkBuckets.UpdateBucketTimes(nodeTime);

        myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);
        totalTimeDijkstras = totalTimeDijkstras + timeDijkstras;
        myNetworkBuckets.UpdateDSpecFeq(myNetworkBuckets.GetTargetFailTime(critNum));
    }

    auto finishLoop = chrono::high_resolution_clock::now();

    myNetworkBuckets.CalculateDSpec(numReps);
    myNetworkBuckets.CalculateReliabilityWeibull(time, 1.5, 10.0);

    double i = 1.0;
    double curShape = 1.5;
    double curScale = 10.0;
    double maxDelta = 10.0;

    ///////////////////////////////////////////////

    auto startPolicy = chrono::high_resolution_clock::now();

//    cout << " Shape = " << curShape << " , Scale = " << curScale << " F Cost = " << cFixed << " v cost = " << cVar << endl;
    cout << " Stable Maintenance Reliability for changing delta " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateStableMaintenanceReliabilityWeibull(i, curShape, curScale);
        i = i + 0.1;
    }


    i = 1.0;
    cout << " Stable Maintenance Variance for changing delta " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateStableMaintenanceVarianceWeibull(i, numReps, curShape, curScale);
        i = i + 0.5;
    }

    i = 1.0;
    cout << "Stable maintenance cost " << endl;
    while (i < 10.0){
        myNetworkBuckets.CalculateMaintenancePolicyCostWeibull(i, cFixed, cVar, curShape, curScale);
        i = i + 0.1;
    }

    myNetworkBuckets.OutputDSpec();

    ////////////////////////////////////////////////

/*    int currNum = nSens - 1;
    while (currNum >= 450) {
        // + 1 to include the sink node
        myNetworkBuckets.SignatureRelation(currNum + 1, time, maxDelta, cFixed, cVar, curShape, curScale);
        currNum = currNum - 1;
    }*/

//        myNetworkBuckets.SignatureRelationOutput(nSens, time, maxDelta, cFixed, cVar, curShape, curScale);

    auto endPolicy = chrono::high_resolution_clock::now();
    chrono::duration<double> timeEvaluatingMaintPolicies = endPolicy - startPolicy;

//        myNetworkBuckets.SignatureRelationVariance(nSens, numReps, time, maxDelta, cFixed, cVar, curShape, curScale);

    auto endPolicyVar = chrono::high_resolution_clock::now();
    chrono::duration<double> timeEvaluatingMaintVar = endPolicyVar - endPolicy;

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finishN - startN;

    chrono::duration<double> dTime = finishLoop - startNetwork;

    cout << " Time Spent Generating Network : " << networkGenerationTime << endl;
    cout << " Time Spent implementing Dijkstras : " << totalTimeDijkstras << endl;
    cout << " Time Generating Fail Time : " << simTimeTotal << endl;
    cout << " Total Time for destruction algorithm : " << dTime.count() << endl;
    cout << " Time evaluating TBM policies : " << timeEvaluatingMaintPolicies.count() << endl;
    cout << " Time estimating Variance : " << timeEvaluatingMaintVar.count() << endl;
    cout << numReps << " - " << " TimeDijkstras " << elapsed.count() << endl;

}

// Monte Carlo method to estimate network reliability
void MonteCarloReliability(int nSens, int numReps, double time, double coverageLevel, double cFixed, double cVar){

    vector<Edge*> networkEdges;
    vector<double > nodeTime;

    // start time record
    auto startN = chrono::high_resolution_clock::now();

    // number of nodes includes sink node, number of sensors, and number of targets
    int numTargets = 441;
    int numNodes = nSens + numTargets + 1;
    double timeDijkstras;
    double totalTimeDijkstras = 0;
    double networkGenerationTime = 0;

    NodeLinkedListMC myNetworkBuckets(numNodes, (numNodes - numTargets), numTargets);
    myNetworkBuckets.SetRanges(0.0750000, 0.0750000);
    myNetworkBuckets.UniformTargetLoc();

    myNetworkBuckets.GenerateNewFailTime(nodeTime);

    auto startNetwork = chrono::high_resolution_clock::now();
    myNetworkBuckets.GenerateNewRGG_Bucket();

    auto finishNetwork = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedNetwork = finishNetwork - startNetwork;
    networkGenerationTime = networkGenerationTime + elapsedNetwork.count();

    myNetworkBuckets.UpdateFailTimes(nodeTime);

    myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);
    totalTimeDijkstras = totalTimeDijkstras + timeDijkstras;

    int critNum = myNetworkBuckets.CoverageNumber(coverageLevel, numTargets);
    myNetworkBuckets.addNetworkFailTime(myNetworkBuckets.GetTargetFailTime(critNum));

    // Start of updating node failure times
    for (int k = 1; k < numReps ; ++k) {
        myNetworkBuckets.GenerateNewFailTime(nodeTime);

        startNetwork = chrono::high_resolution_clock::now();
        myNetworkBuckets.GenerateNewRGG_Bucket();

        finishNetwork = chrono::high_resolution_clock::now();
        elapsedNetwork = finishNetwork - startNetwork;
        networkGenerationTime = networkGenerationTime + elapsedNetwork.count();

        myNetworkBuckets.UpdateFailTimes(nodeTime);

        myNetworkBuckets.DijkstrasAlgorithm(networkEdges, numNodes, numTargets, timeDijkstras);
        myNetworkBuckets.addNetworkFailTime(myNetworkBuckets.GetTargetFailTime(critNum));

        totalTimeDijkstras = totalTimeDijkstras + timeDijkstras;
    }

    myNetworkBuckets.CalculateReliability(time, numReps);

    double i = 1.0;
    double curShape = 1.5;
    double curScale = 10.0;
    double maxDelta = 10.0;

    ///////////////////////////////////////////////

    auto startPolicy = chrono::high_resolution_clock::now();

    auto endPolicy = chrono::high_resolution_clock::now();
    chrono::duration<double> timeEvaluatingMaintPolicies = endPolicy - startPolicy;

    // end time record
    auto finishN = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finishN - startN;

    cout << " Time Spent Generating Network : " << networkGenerationTime << endl;
    cout << " Time Spent implementing Dijkstras : " << totalTimeDijkstras << endl;
    cout << " Time evaluating TBM policies : " << timeEvaluatingMaintPolicies.count() << endl;
    cout << numReps << " - " << " TimeDijkstras " << elapsed.count() << endl;

}

// Used to plot the transient reliability of a given (n1, delta) TBM policy
void MonteCarloTBM(int nSens, int numReps, double time, double coverageLevel, double timeInterval, double delta,
                   double cFixed, double cVar){

    vector<Edge*> networkEdges;
    vector<double> nodeTime;
    vector<int> failCount(4001,0);
    vector<double> maintCost;

//    vector<int> tempRep;

    // number of nodes includes sink node, number of sensors, and number of targets
    int numTargets = 441;
    int numNodes = nSens + numTargets + 1;
    double policyCost = 0;
    double timeDijkstras;
    double totalTimeDijkstras = 0;

    NodeLinkedListMC myNetworkBuckets(numNodes, (numNodes - numTargets), numTargets);

    auto startTime = chrono::high_resolution_clock::now();

    for (int k = 1; k < numReps ; ++k) {

        myNetworkBuckets.SetRanges(0.0750000, 0.0750000);
        myNetworkBuckets.UniformTargetLoc();

        // Generate initial fail times
        myNetworkBuckets.GenerateNewFailTime(nodeTime);

        // generate initial RGG
        myNetworkBuckets.GenerateNewRGG_Bucket();
        myNetworkBuckets.UpdateFailTimes(nodeTime);

        // implement BFS to determine the number of targets covered
        double tarCovered = myNetworkBuckets.BreadthFirstSearchReturn(numNodes, numTargets, timeDijkstras);

        if (tarCovered < coverageLevel){
            failCount[0] = failCount[0] + 1;
        }

        double timeToMaint = delta;

        // update rgg (shouldn't change at this point)
        myNetworkBuckets.UpdateRGG(0);

        double curTime = timeInterval;
        int curIndex = 1;

        timeToMaint = timeToMaint - timeInterval;

        while (curTime < time) {

            // Time for maintenance
            if (timeToMaint <= 0.0001){

                // redeploy sensors, and generate a new fail time for those added
                myNetworkBuckets.GenerateNewFailTimeMaintenance(nodeTime, time, curTime);

                // randomly replace failed nodes
                myNetworkBuckets.GenerateNewRGG_BucketMaint(cFixed, cVar, policyCost);
//                myNetworkBuckets.GenerateNewRGG_BucketMaint();
                maintCost.push_back(policyCost / delta);
//                tempRep.push_back(k);

                myNetworkBuckets.UpdateFailTimes(nodeTime);
                tarCovered = myNetworkBuckets.BreadthFirstSearchReturn(numNodes, numTargets, timeDijkstras);

                timeToMaint = delta - timeInterval;

                if (tarCovered < coverageLevel){
                    failCount[curIndex] = failCount[curIndex] + 1;
                }

            }
            else{

                myNetworkBuckets.UpdateRGG(curTime);
                myNetworkBuckets.UpdateFailTimes(nodeTime);
                tarCovered = myNetworkBuckets.BreadthFirstSearchReturn(numNodes, numTargets, timeDijkstras);

                if (tarCovered < coverageLevel){
                    failCount[curIndex] = failCount[curIndex] + 1;
                }

                timeToMaint = timeToMaint - timeInterval;

            }

            curTime = curTime + timeInterval;
            curIndex = curIndex + 1;
        }

    }

    ///////////////////////////////////////////////

    double totalCost = 0;
    for (int j = 0; j < maintCost.size(); ++j) {
        totalCost = totalCost + maintCost.at(j);
    }

    double averageCost = totalCost / (double) maintCost.size();

    double timeIndex = 0;

    for (int i = 0; i < failCount.size() ; ++i) {
        cout << " Time " << timeIndex << " -- NumRepsFailed " << failCount.at(i) << " outOf " << numReps << endl;
        timeIndex = timeIndex + timeInterval;
    }
//
//    ofstream myFile;
//    myFile.open ("MaintReliability.txt");
//
//    cout << " -------- Actual Maintenance Costs ------------------" << endl;
//    for (int l = 0; l <tempRep.size() ; ++l) {
//
//        myFile << " Rep : " << tempRep.at(l) << " -- cost : " << maintCost.at(l) << endl;
//
//    }
//
//    myFile.close();

    auto finishTime = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedNetwork = finishTime - startTime;

    cout << " Time Spent transient reliability : " << elapsedNetwork.count() << endl;
    cout << " Long run average cost : " << averageCost << endl;
    cout << " Time Spent implementing Dijkstras : " << totalTimeDijkstras << endl;
}

void TestShuffle(char *fileName){

    vector<Edge*> networkEdges;
    vector<int> nodeTime;
    int numTargets;

    ReadNetwork(fileName, networkEdges, nodeTime, numTargets);

    int numNodes = nodeTime.size();

    NodeLinkedList myNetworkBuckets(numNodes, (numNodes - numTargets), numTargets, nodeTime);

    for (int i = 1; i < numNodes - numTargets; ++i) {
        cout << nodeTime.at(i) << " " ;
    }
    cout << endl;

    for (int j = 0; j < 10000 ; ++j) {
        myNetworkBuckets.GenerateNewFailOrder(nodeTime);
        for (int i = 1; i < numNodes - numTargets; ++i) {
            cout << nodeTime.at(i) << " " ;
        }
        cout << endl;
    }
}

// TODO: Make sure the sensor life distributions in Network Reliability and ReliabilityVariance match
int main(int argc, char* argv[]) {

    int numReps = 50000;
    double t = 40.0;
    double coverage = 0.80;

    double cfixed = 100;
    double cVar = 1.0;


//    RunAlgorithms("C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\TestFinal.txt");
//    UpdateInstance("C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\TestFinal.txt","C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\TestFinalUpdate.txt");
//    UpdateInstance("C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\TestInstance550.txt","C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\TestInstanceUpdate550p1.txt");

//    GenerateInstancesDials("C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\test750.txt", numReps, t, coverage, cfixed, cVar);
//    GenerateInstancesDijkstras("C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\test750.txt", numReps, t, coverage, cfixed, cVar);


//        DialsRandomNetwork(450, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(500, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(550, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(600, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(650, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(700, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(750, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(800, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(850, numReps, t, coverage, cfixed, cVar);
//        DialsRandomNetwork(900, numReps, t, coverage, cfixed, cVar);

//    DialsRandomNetwork(850, numReps, t, coverage, cfixed, cVar);
//   DijkstrasRandomNetwork(850, 20000, t, coverage, cfixed, cVar);



//        MonteCarloReliability(600, 50000, t, coverage, cfixed, cVar);

    MonteCarloTBM(650, 1000, 400, coverage, 0.1, 5.6, cfixed, cVar);



//TestDSpec("C:\\Users\\nboardma\\Desktop\\test550.txt", 550, t, 1.5,  10);


//    GenerateInstancesDials("C:\\Users\\nboardma\\Documents\\Personal\\CPPFiles\\test450.txt", numReps, t, coverage, cfixed, cVar);


/////////////////////////////////////////////////////////
// Use this format to run from the command line, where the file name is input

//    UpdateInstance(argv[1], argv[2]);

//    GenerateInstancesDials(argv[1], numReps, t, coverage, cfixed, cVar);
//cout <<" Coverage level 95" << endl;
//    GenerateInstancesDials(argv[1], numReps, t, 0.95);
//
//    cout <<" Coverage level 90" << endl;
//    GenerateInstancesDials(argv[1], numReps, t, 0.90);
//
//    cout <<" Coverage level 85" << endl;
//    GenerateInstancesDials(argv[1], numReps, t, 0.85);
//
//    cout <<" Coverage level 80" << endl;
//    GenerateInstancesDials(argv[1], numReps, t, 0.80);


//    GenerateInstancesDijkstras(argv[1], numReps, t, coverage);
//    cout <<" Coverage level 95" << endl;
//    GenerateInstancesDijkstras(argv[1], numReps, t, 0.95);
//
//    cout <<" Coverage level 90" << endl;
//    GenerateInstancesDijkstras(argv[1], numReps, t, 0.90);
//
//    cout <<" Coverage level 85" << endl;
//    GenerateInstancesDijkstras(argv[1], numReps, t, 0.85);
//
//    cout <<" Coverage level 80" << endl;
//    GenerateInstancesDijkstras(argv[1], numReps, t, 0.80);



//    MultipleInstancesDials(argc, argv);
//    MultipleInstancesDijkstras(argc, argv);


    return 0;
}