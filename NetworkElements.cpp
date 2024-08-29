//
// Created by nboardma on 7/18/2019.
//

#include "NetworkElements.h"
#include "Reliability.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>

using namespace std;

// constructor for NodeLinkedList
NodeLinkedList::NodeLinkedList(){
    firstItem = nullptr;
    lastItem = nullptr;
}
NodeLinkedList::NodeLinkedList(int numItems){
    firstItem = nullptr;
    lastItem = nullptr;
    myNumElements = numItems;
    InitializeBuckets(numItems);
}
NodeLinkedList::NodeLinkedList(int numItems, vector<int> order){
    firstItem = nullptr;
    lastItem = nullptr;
    myNumElements = numItems;
    InitializeBuckets(numItems, order);
}

// use when no input file is to be provided, RGG and fail order can be generated directly in c++
NodeLinkedList::NodeLinkedList(int numItems, int numSens, int numTarg){
    firstItem = nullptr;
    lastItem = nullptr;
    myNumElements = numItems;   // total number of nodes
    myNumSensors = numSens;     // number of sensors (including sink)
    myNumTargets = numTarg;     // number of targets
    InitializeBuckets(numItems);
}

NodeLinkedList::NodeLinkedList(int numItems, int numSens, int numTarg, vector<int> order){
    firstItem = nullptr;
    lastItem = nullptr;
    myNumElements = numItems;   // total number of nodes
    myNumSensors = numSens;     // number of sensors (including sink)
    myNumTargets = numTarg;     // number of targets
    InitializeBuckets(numItems, order);
}

// initializes buckets. To begin, all nodes are in bucket zero and point to each other in sequential order
void NodeLinkedList::InitializeBuckets(int numItems){

    DSpecFreq.assign(myNumSensors, 0);

    myBuckets = new Node[numItems];

    // update all values for the sink node
    Node *ptr = myBuckets;
    ptr->SetID(0);
    ptr->UpdateCapacityLabel(numItems);
    ptr->UpdateReached();
    ptr->MarkPermanent();
    ptr->SetFirstItem((ptr + 1));
    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    for (int i = 1; i < numItems ; ++i) {
        myBuckets[i].SetID(i);
        myBuckets[i].SetNext((ptr + 1));
        myBuckets[i].SetPrev((ptr - 1));
        myBuckets[i].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }
    myBuckets[numItems-1].SetNext(nullptr);

}

// initializes buckets. To begin, all nodes are in bucket zero and point to each other in sequential order
// also sets the value for each node according to the input failure order in vector<int> failTime
void NodeLinkedList::InitializeBuckets(int numItems, vector<int> failOrder){

    DSpecFreq.assign(myNumSensors, 0);

    myBuckets = new Node[numItems];

    Node *ptr = myBuckets;
    // update all values for the sink node
    ptr->SetID(0);
    ptr->UpdateCapacityLabel(numItems);
    ptr->UpdateReached();
    ptr->MarkPermanent();
    ptr->SetValue(failOrder[0]);
    ptr->SetFirstItem((ptr + 1));

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    // values for sensor nodes
    for (int i = 1; i < myNumSensors ; ++i) {
        myBuckets[i].SetID(i);
        myBuckets[i].SetNext((ptr + 1));
        myBuckets[i].SetPrev((ptr - 1));

        myBuckets[failOrder[i]].SetValue(i);

        myBuckets[i].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    // values for target nodes
    for (int j = myNumSensors ; j < numItems ; ++j) {

        myBuckets[j].SetID(j);
        myBuckets[j].SetNext((ptr + 1));
        myBuckets[j].SetPrev((ptr - 1));

        myBuckets[j].SetValue(myNumElements);

        myBuckets[j].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    myBuckets[numItems-1].SetNext(nullptr);

}

// removes the first element from the list input by listNum
Node *NodeLinkedList::RemoveFirst(int listNum){

    Node *tempFirst = myBuckets[listNum].GetFirstElement();
    Node *tempNext = tempFirst->GetNext();

    // update the pointer to the first element on the list
    myBuckets[listNum].SetFirstItem(tempNext);

    // update the pointer of the new first element in the list
    if (tempNext != nullptr){
        tempNext->SetPrev(nullptr);
    }

    // update the Next pointer for element removed
    tempFirst->SetNext(nullptr);

    return tempFirst;

}

// removes a specific element from an input list
Node *NodeLinkedList::RemoveTargetElement(int listNum, Node *removeElement){

    // temporarily store Prev and Next pointer for the element removed
    Node *tempPrev = removeElement->GetPrev();
    Node *tempNext = removeElement->GetNext();

    // checks to see if the element removed is the first element in the list
    if (tempPrev == nullptr){
        myBuckets[listNum].SetFirstItem(tempNext);
    }

    // update Previous and Next pointers for nodes on either side of the removed node
    if (tempPrev != nullptr){
        tempPrev->SetNext(tempNext);
    }

    if (tempNext != nullptr){
        tempNext->SetPrev(tempPrev);
    }

    // remove pointers from the element removed
    removeElement->SetPrev(nullptr);
    removeElement->SetNext(nullptr);

    return removeElement;

}

// Adds a new element newElement to the front of an input list
void NodeLinkedList::AddToList(int listNum, Node *newElement){

    // store the original first item in the list
    Node *temp = myBuckets[listNum].GetFirstElement();

    // if the current first item in the list is null, place item in the front of the list
    if (temp == nullptr) {
        newElement->SetNext(nullptr);
    }
    // otherwise update the pointer for the original first item to have pointer prev equal to added element
    // update pointer next for the newly added element to be the original element
    else {
        newElement->SetNext(temp);
        temp->SetPrev(newElement);
    }

    // update pointer to the first element element in the list
    myBuckets[listNum].SetFirstItem(newElement);

    newElement->SetPrev(nullptr);

}

// returns true if list is empty, or false if there are elements in the list
bool NodeLinkedList::IsListEmpty(int listNum){
    Node *tempFirst = myBuckets[listNum].GetFirstElement();
    return (tempFirst == nullptr);
}

// Build the adjacency list for each node
void NodeLinkedList::BuildForwardStar(vector<Edge *> edgeList) {
    for (int i = 0; i < edgeList.size() ; ++i) {
        myBuckets[edgeList[i]->myStartNode].AddReachableNode(edgeList[i]->myEndNode);
    }
}

// performs Dials algorithm on a network. The Sink node is Node 0, sensor nodes 1,2,..., n
void NodeLinkedList::DialsAlgorithm(vector<Edge*> edgeList, int numNodes) {

    int numPermanent = 1;

    // Build the adjacency list for each node
    for (int i = 0; i < edgeList.size() ; ++i) {
        myBuckets[edgeList[i]->myStartNode].AddReachableNode(edgeList[i]->myEndNode);
    }

    // start time record
    auto start = chrono::high_resolution_clock::now();

    // iterates through edges in the network that leave the sink node, and update the bucket and max capacity path label
    // of the newly reached node
    Node sinkNode = myBuckets[0];

    for (int j = 0; j < sinkNode.GetDegree(); ++j) {

        int updateNode = sinkNode.GetReachableNode(j);
        int updateLabel = min(sinkNode.GetValue(), myBuckets[updateNode].GetValue());

        myBuckets[updateNode].UpdateReached();
        myBuckets[updateNode].UpdateCapacityLabel(updateLabel);

        Node *ptr = myBuckets;
        ptr = ptr + updateNode;

        this->RemoveTargetElement(0, ptr);
        this->AddToList(updateLabel, ptr);

    }

    int currBucket = numNodes - 1;

    while (numPermanent < numNodes && currBucket > 0){

        // if the current bucket is not empty
        while (!(this->IsListEmpty(currBucket))){
            // remove the next node from the bucket

            Node *permNode = myBuckets;
            permNode = permNode + (this->RemoveFirst(currBucket))->GetNodeID();
            permNode->MarkPermanent();
            this->addFailOrder(permNode->GetNodeID());
            this->addFailTime(permNode->GetMaxCapacityPath());
            numPermanent = numPermanent + 1;

            // iterate through edges
            for (int i = 0; i < permNode->GetDegree(); ++i) {

                int updateNode1 = permNode->GetReachableNode(i);

                if (!myBuckets[updateNode1].IsPermanent()){
                    myBuckets[updateNode1].UpdateReached();

                    // determine the current label and the new label for the reachable node

                    int currLabel =  myBuckets[updateNode1].GetMaxCapacityPath();
                    int updateLabel = min(permNode->GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());

                    //TODO : think about adding an IF statement to see if updateLabel == node current max capacity
                    // if it is then the AddToList does not need to be reaccomplished. Right now it places the node
                    // at the beginning of the list anyway

                    Node *ptr = myBuckets;
                    ptr = ptr + updateNode1;

                    this->RemoveTargetElement(currLabel, ptr);

                    myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);

                    this->AddToList(updateLabel, ptr);
                }
            }
        }
        currBucket = currBucket - 1;
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;
}

// performs Dials algorithm on a network. The Sink node is Node 0, sensor nodes 1,2,..., n
// and target nodes n+1, n+2, ... , n + numTargets
void NodeLinkedList::DialsAlgorithm(vector<Edge*> edgeList, int numNodes, int numTar, double &timeDuration) {

    // numNodes - Includes every node in the network (sink, sensors, and targets)
    int numPermanent = 1;

    // start time record
    auto start = chrono::high_resolution_clock::now();

    // iterates through edges in the network that leave the sink node, and update the bucket and max capacity path label
    // of the newly reached node
    Node sinkNode = myBuckets[0];

    for (int j = 0; j < sinkNode.GetDegree(); ++j) {

        int updateNode = sinkNode.GetReachableNode(j);
        int updateLabel = min(sinkNode.GetValue(), myBuckets[updateNode].GetValue());

        myBuckets[updateNode].UpdateReached();
        myBuckets[updateNode].UpdateCapacityLabel(updateLabel);

        Node *ptr = myBuckets;
        ptr = ptr + updateNode;

        this->RemoveTargetElement(0, ptr);
        this->AddToList(updateLabel, ptr);

    }

    int currBucket = numNodes - numTar;

    while (numPermanent < numNodes && currBucket > 0){

        // if the current bucket is not empty
        while (!(this->IsListEmpty(currBucket))){
            // remove the next node from the bucket

            Node *permNode = myBuckets;
            permNode = permNode + (this->RemoveFirst(currBucket))->GetNodeID();
            permNode->MarkPermanent();
            // stores the failure order and time at which node is disconnected
            this->addFailOrder(permNode->GetNodeID());
            this->addFailTime(permNode->GetMaxCapacityPath());
            numPermanent = numPermanent + 1;

            // stores the failure order and time at which Targets are disconnected
            if ( permNode->GetNodeID() >= numNodes - numTar ) {
                this->addTargetFailOrder(permNode->GetNodeID());
                this->addTargetFailTime(permNode->GetMaxCapacityPath());
            }
            else {
                this->addSensorFailOrder(permNode->GetNodeID());
                this->addSensorFailTime(permNode->GetMaxCapacityPath());
            }

            // iterate through edges
            for (int i = 0; i < permNode->GetDegree(); ++i) {

                int updateNode1 = permNode->GetReachableNode(i);

                if (!myBuckets[updateNode1].IsPermanent()){
                    myBuckets[updateNode1].UpdateReached();

                    // determine the current label and the new label for the reachable node

                    int currLabel =  myBuckets[updateNode1].GetMaxCapacityPath();
                    int updateLabel = min(permNode->GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());

                    //TODO : think about adding an IF statement to see if updateLabel == node current max capacity
                    // if it is then the AddToList does not need to be reaccomplished. Right now it places the node
                    // at the beginning of the list anyway

                    Node *ptr = myBuckets;
                    ptr = ptr + updateNode1;

                    this->RemoveTargetElement(currLabel, ptr);

                    myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);

                    this->AddToList(updateLabel, ptr);
                }
            }
        }
        currBucket = currBucket - 1;
    }

    // If a target is not connected to the sink then it will never be marked permanent, so it will never be added to the
    // order of target failures. This for loop checks for such a case, and adds the targets to ths list with a fail time of zero
    if (orderedTargetFail.size() < numTar){
        Node *ptr = myBuckets;
        ptr = ptr + (numNodes - numTar);
        for (int i = 0; i < numTar ; ++i) {
            if (!ptr->IsPermanent()){
                this->addTargetFailOrder(ptr->GetNodeID());
                this->addTargetFailTime(ptr->GetMaxCapacityPath());
            }
            ptr = ptr + 1;
        }
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;
    timeDuration = elapsed.count();

}

// This function performs Dials Algorithm the same as above, but uses the max(min()) rule to update
// The results are exactly the same
void NodeLinkedList::DialsAlgorithmUpdateRule(vector<Edge *> edgeList, int numNodes, int numTar, double &timeDuration) {

    // numNodes - Includes every node in the network (sink, sensors, and targets)
    int numPermanent = 1;

    // start time record
    auto start = chrono::high_resolution_clock::now();

    // iterates through edges in the network that leave the sink node, and update the bucket and max capacity path label
    // of the newly reached node
    Node sinkNode = myBuckets[0];

    for (int j = 0; j < sinkNode.GetDegree(); ++j) {

        int updateNode = sinkNode.GetReachableNode(j);
        int updateLabel = max(myBuckets[updateNode].GetMaxCapacityPath(), min(sinkNode.GetMaxCapacityPath(), myBuckets[updateNode].GetValue()));

        myBuckets[updateNode].UpdateReached();
        myBuckets[updateNode].UpdateCapacityLabel(updateLabel);

        Node *ptr = myBuckets;
        ptr = ptr + updateNode;

        this->RemoveTargetElement(0, ptr);
        this->AddToList(updateLabel, ptr);

    }

    int currBucket = numNodes - numTar;

    while (numPermanent < numNodes && currBucket > 0){

        // if the current bucket is not empty
        while (!(this->IsListEmpty(currBucket))){
            // remove the next node from the bucket

            Node *permNode = myBuckets;
            permNode = permNode + (this->RemoveFirst(currBucket))->GetNodeID();

            permNode->MarkPermanent();
            // stores the failure order and time at which node is disconnected
            this->addFailOrder(permNode->GetNodeID());
            this->addFailTime(permNode->GetMaxCapacityPath());
            numPermanent = numPermanent + 1;

            // stores the failure order and time at which Targets are disconnected
            if ( permNode->GetNodeID() >= numNodes - numTar ) {
                this->addTargetFailOrder(permNode->GetNodeID());
                this->addTargetFailTime(permNode->GetMaxCapacityPath());
            }
            else {
                this->addSensorFailOrder(permNode->GetNodeID());
                this->addSensorFailTime(permNode->GetMaxCapacityPath());
            }

            // iterate through edges
            for (int i = 0; i < permNode->GetDegree(); ++i) {

                int updateNode1 = permNode->GetReachableNode(i);

                int currLabel =  myBuckets[updateNode1].GetMaxCapacityPath();
                int updateLabel = max(myBuckets[updateNode1].GetMaxCapacityPath(), min(permNode->GetMaxCapacityPath(), myBuckets[updateNode1].GetValue()));

                Node *ptr = myBuckets;
                ptr = ptr + updateNode1;

                if (currLabel == updateLabel){
                    // prevents the situation in which a node is placed back in the same bucket (particularly if the node has already been permanent)
                }
                else {
                    this->RemoveTargetElement(currLabel, ptr);
                    myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);
                    this->AddToList(updateLabel, ptr);

                }
            }
        }
        currBucket = currBucket - 1;
    }

    // If a target is not connected to the sink then it will never be marked permanent, so it will never be added to the
    // order of target failures. This for loop checks for such a case, and adds the targets to ths list with a fail time of zero
    if (orderedTargetFail.size() < numTar){
        Node *ptr = myBuckets;
        ptr = ptr + (numNodes - numTar);
        for (int i = 0; i < numTar ; ++i) {
            if (!ptr->IsPermanent()){
                this->addTargetFailOrder(ptr->GetNodeID());
                this->addTargetFailTime(ptr->GetMaxCapacityPath());
            }
            ptr = ptr + 1;
        }
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;
    timeDuration = elapsed.count();

}


// performs Dijkstra's algorithm on a network. The Sink node is Node 0, sensor nodes 1,2,..., n
void NodeLinkedList::DijkstrasAlgorithm(vector<Edge *> edgeList, int numNodes) {

    // sink node is already permanent
    int numPermanent = 1;

    // vector to store ID's of temporary nodes
    vector<int> tempNodes;
    for (int i = 1; i <= (numNodes-1) ; ++i) {
        tempNodes.push_back(i);
    }

    // Build the adjacency list for each node
    for (int i = 0; i < edgeList.size() ; ++i) {
        myBuckets[edgeList[i]->myStartNode].AddReachableNode(edgeList[i]->myEndNode);
    }

    // start time record
    auto start = chrono::high_resolution_clock::now();

    // update all nodes connected to sink
    for (int i = 0; i < myBuckets[0].GetDegree(); ++i) {
        int updateNode1 = myBuckets[0].GetReachableNode(i);
        int updateLabel = min( myBuckets[0].GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());
        myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);
    }

    while (numPermanent < numNodes){

        int tempNode;
        int tempIndex = -1;
        int currMax = 0;
        int permNode = 0;

        // find the node with the largest distance label that is not permanent
        for (int j = 0; j < tempNodes.size(); ++j) {
            tempNode = tempNodes.at(j);
            if (myBuckets[tempNode].GetMaxCapacityPath() > currMax) {
                // store node ID and max label

                currMax = myBuckets[tempNode].GetMaxCapacityPath();
                permNode = myBuckets[tempNode].GetNodeID();
                tempIndex = j;
            }
        }

        if (tempIndex == -1) {
            // node not connected
        } else {
            // remove from temporary nodes vector
            tempNodes.erase(tempNodes.begin() + tempIndex);

            // mark node permanent
            myBuckets[permNode].MarkPermanent();

            this->addFailOrder(permNode);
            this->addFailTime(currMax);

            // iterate through edges
            for (int k = 0; k < myBuckets[permNode].GetDegree(); ++k) {

                int updateNode1 = myBuckets[permNode].GetReachableNode(k);
                if (!myBuckets[updateNode1].IsPermanent()) {
                    int updateLabel = min(myBuckets[permNode].GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());
                    myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);
                }
            }
        }
        numPermanent = numPermanent + 1;
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;

/*
    // sink node is already permanent
    int numPermanent = 1;

    // Build the adjacency list for each node
    for (int i = 0; i < edgeList.size() ; ++i) {
        myBuckets[edgeList[i]->myStartNode].AddReachableNode(edgeList[i]->myEndNode);
    }

    // update all nodes connected to sink
    for (int i = 0; i < myBuckets[0].GetDegree(); ++i) {
        int updateNode1 = myBuckets[0].GetReachableNode(i);
        int updateLabel = min( myBuckets[0].GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());
        myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);
    }

    while (numPermanent < numNodes){
        int currMax = 0;
        int permNode = 0;
        Node *ptr = myBuckets;
        // find the node with the largest distance label that is not permanent
        for (int j = 0; j < numNodes; ++j) {
            if ((ptr->GetMaxCapacityPath() > currMax) && !ptr->IsPermanent()) {
                // store node ID and max label
                currMax = ptr->GetMaxCapacityPath();
                permNode = ptr->GetNodeID();
            }

            ptr = ptr + 1;
        }

        // mark node permanent
        myBuckets[permNode].MarkPermanent();

        this->addFailOrder(permNode);
        this->addFailTime(currMax);

        // iterate through edges
        for (int k = 0; k < myBuckets[permNode].GetDegree(); ++k) {

            int updateNode1 = myBuckets[permNode].GetReachableNode(k);
            if (!myBuckets[updateNode1].IsPermanent()){
                int updateLabel = min( myBuckets[permNode].GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());
                myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);
            }
        }

        numPermanent = numPermanent + 1;
    }

*/

}

// performs Dials algorithm on a network. The Sink node is Node 0, sensor nodes 1,2,..., n
// and target nodes n+1, n+2, ... , n + numTargets
void NodeLinkedList::DijkstrasAlgorithm(vector<Edge *> edgeList, int numNodes, int numTar, double &timeDuration) {

    // sink node is already permanent
    int numPermanent = 1;

    // vector to store ID's of temporary nodes
    vector<int> tempNodes;
    for (int i = 1; i <= (numNodes-1) ; ++i) {
        tempNodes.push_back(i);
    }

    // start time record
    auto start = chrono::high_resolution_clock::now();

    // update all nodes connected to sink
    for (int i = 0; i < myBuckets[0].GetDegree(); ++i) {
        int updateNode1 = myBuckets[0].GetReachableNode(i);
        int updateLabel = min( myBuckets[0].GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());
        myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);
    }

    while (numPermanent < numNodes) {

        int tempNode;
        int tempIndex = -1;
        int currMax = 0;
        int permNode = 0;

        // find the node with the largest distance label that is not permanent
        for (int j = 0; j < tempNodes.size(); ++j) {
            tempNode = tempNodes.at(j);
            if (myBuckets[tempNode].GetMaxCapacityPath() > currMax) {
                // store node ID and max label

                currMax = myBuckets[tempNode].GetMaxCapacityPath();
                permNode = myBuckets[tempNode].GetNodeID();
                tempIndex = j;
            }
        }

        if (tempIndex == -1) {
            // node not connected - should be able to break here
        } else {

            // remove from temporary nodes vector
            tempNodes.erase(tempNodes.begin() + tempIndex);

            // mark node permanent
            myBuckets[permNode].MarkPermanent();

            this->addFailOrder(permNode);
            this->addFailTime(currMax);

            // stores the failure order and time at which Targets are disconnected
            if (myBuckets[permNode].GetNodeID() >= numNodes - numTar) {
                this->addTargetFailOrder(myBuckets[permNode].GetNodeID());
                this->addTargetFailTime(myBuckets[permNode].GetMaxCapacityPath());
            }
            else {
                this->addSensorFailOrder(myBuckets[permNode].GetNodeID());
                this->addSensorFailTime(myBuckets[permNode].GetMaxCapacityPath());
            }

            // iterate through edges
            for (int k = 0; k < myBuckets[permNode].GetDegree(); ++k) {

                int updateNode1 = myBuckets[permNode].GetReachableNode(k);
                if (!myBuckets[updateNode1].IsPermanent()) {
                    int updateLabel = min(myBuckets[permNode].GetMaxCapacityPath(), myBuckets[updateNode1].GetValue());
                    myBuckets[updateNode1].UpdateCapacityLabel(updateLabel);
                }
            }
        }
        numPermanent = numPermanent + 1;
    }

    // If a target is not connected to the sink then it will never be marked permanent, so it will never be added to the
    // order of target failures. This for loop checks for such a case, and adds the targets to ths list with a fail time of zero
    if (orderedTargetFail.size() < numTar){
        Node *ptr = myBuckets;
        ptr = ptr + (numNodes - numTar);
        for (int i = 0; i < numTar ; ++i) {
            if (!ptr->IsPermanent()){
                this->addTargetFailOrder(ptr->GetNodeID());
                this->addTargetFailTime(ptr->GetMaxCapacityPath());
            }
            ptr = ptr + 1;
        }
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;
    timeDuration = elapsed.count();
}

// within an iteration, add the number of sensors failed that cause network failure
// DSpecFreq - vector of size (numSensors x 1) with frequency for ith sensor failure causing network failure
// DSpecElements - vector containing unique elements of DSpecFreq with elements greater than zero
void NodeLinkedList::UpdateDSpecFeq(int numFailed) {

    if (DSpecFreq.at(numFailed) == 0) {
        DSpecFreq[numFailed] = 1;
        DSpecElements.push_back(numFailed);
    }
    else {
        DSpecFreq[numFailed] = DSpecFreq[numFailed] + 1;
    }
}

// Once all replications are complete, calculate the D Spectrum. The number of failures is also in order
// myDSpec - number of sensors failed that cause network failure
// myDSpecProb - probability that the corresponding entry in myDSpec causes network failure
void NodeLinkedList::CalculateDSpec(int numReps){
    for (int i = 0; i < DSpecFreq.size() ; ++i) {
        if (DSpecFreq.at(i) == 0){
            // skip
        }
        else {
            myDSpec.push_back(i);
            double myProb = (double) DSpecFreq.at(i) / (double) numReps;
            myDSpecProb.push_back(myProb);
        }
    }
}

void NodeLinkedList::OutputDSpec() {
    for (int i = 0; i < myDSpec.size() ; ++i) {
        cout << myDSpec.at(i) << "     ----    " << myDSpecProb.at(i) << endl;
    }
}

// outputs critical times for every node in the network
void NodeLinkedList::OutputCriticalTimes(){
    for (int i = 0; i < myNumElements; ++i) {
        cout << " Node Number " << myBuckets[i].GetNodeID() << " --- Critical Time " << myBuckets[i].GetMaxCapacityPath() << "\n";
    }
}

// outputs critical times for targets only, in order
void NodeLinkedList::OutputTargetCriticalTimes(int numTargets) {
    cout << "TargetFail ";
    for (int i = (numTargets - 1); i >= 0; --i) {
        cout << orderedTargetCriticalFailTime.at(i) << " ";
    }
    cout << endl;
}

// outputs the target ID in order for which targets are disconnected
void NodeLinkedList::OutputTargetFailOrder(int numTargets){
    cout << "TargetFailOrder ";
    for (int i = (numTargets - 1); i >= 0; --i) {
        cout << orderedTargetFail.at(i) << " ";
    }
    cout << endl;
}

// outputs the sensor ID in order for which sensors are disconnected
void NodeLinkedList::OutputSensorFailOrder(int numSensors){
    cout << "SensorFailOrder ";
    for (int i = (numSensors - 1); i >= 0; --i) {
        cout << orderedSensorFail.at(i) << " " ;
    }
    cout << endl;
}

// vector that stores the ordered failed times
void NodeLinkedList::addFailTime(int failTime){
    orderedCriticalFailTime.push_back(failTime);
}

// vector that stores the order in which nodes are disconnected from the sink (node ID) - for all nodes
void NodeLinkedList::addFailOrder(int nodeNum){
    orderedFail.push_back(nodeNum);
}

// returns the time from an input location from the ordered list of fail times
// NOTE: the list is in reverse order (largest to smallest), does not include sink
// myNumElements-1 because myNumElements includes the sink node
int NodeLinkedList::GetFailTime(int num){
    return orderedCriticalFailTime.at(myNumElements - 1 - num);
}

// returns a number of the ordered list of failed nodes (node ID)
// NOTE: the list is in reverse order (largest to smallest), does not include sink
// myNumElements-1 because myNumElements includes the sink node
int NodeLinkedList::GetFailOrder(int num){
    return orderedFail.at(myNumElements - 1 - num);
}

// vector that stores the ordered failed times for targets
void NodeLinkedList::addTargetFailTime(int failTime){
    orderedTargetCriticalFailTime.push_back(failTime);
}

// vector that stores the order in which targets are disconnected from the sink (target ID)
void NodeLinkedList::addTargetFailOrder(int nodeNum){
    orderedTargetFail.push_back(nodeNum);
}

// vector that stores the ordered fail time for sensors
void NodeLinkedList::addSensorFailTime(int failTime){
    orderedSensorCriticalFailTime.push_back(failTime);
}

// vector that stores the order in which sensors are disconnected from the sink (sensors ID)
void NodeLinkedList::addSensorFailOrder(int nodeNum){
    orderedSensorFail.push_back(nodeNum);
}

int NodeLinkedList::GetSensorFailOrder(int senNum){
    return orderedSensorFail.at(orderedSensorFail.size() - senNum);
}

// returns a number of the ordered list of disconnected targets times
// NOTE: the list is in reverse order (largest to smallest)
int NodeLinkedList::GetTargetFailTime(int failNum){
    return orderedTargetCriticalFailTime.at(orderedTargetCriticalFailTime.size() - failNum);
}

// returns an element from the ordered list of targets that are disconnected
// NOTE: the list is in reverse order (largest to smallest)
int NodeLinkedList::GetTargetFailOrder(int targNum){
    return orderedTargetFail.at(orderedTargetFail.size() - targNum);
}

// determines how many TARGETS need to be disconnected/failed before coverage drops below a desired threshold
// returns this number. Coverage is between 0 and 1
int NodeLinkedList::CoverageNumber(double coverage, int numT) {
    // myNumElements includes the sink node, so subtract 1 to get number of sensor nodes
    int coverNum;
    double exactNum = numT * (1-coverage);
    int wholeNum = (int) ceil(numT * (1-coverage));

    if (abs(wholeNum - exactNum) <= 0.0001){
        coverNum = wholeNum + 1;
    }
    else {
        coverNum = wholeNum;
    }
    return coverNum;
}

// used to reinitialize buckets for a new test instance
void NodeLinkedList::UpdateBucketTimes(vector<int> newOrder){

    this->ClearOrderTimes();

    Node *ptr = myBuckets;
    // update all values for the sink node
    ptr->SetID(0);
    ptr->UpdateCapacityLabel(myNumElements);
    ptr->UpdateReached();
    ptr->MarkPermanent();
    ptr->SetValue(newOrder[0]);
    ptr->SetFirstItem((ptr + 1));

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    for (int i = 1; i < myNumSensors ; ++i) {

        // reset initial status
        myBuckets[i].ResetNodeReached();
        myBuckets[i].ResetNodePermanent();
        myBuckets[i].UpdateCapacityLabel(0);

        // re-initialize
        myBuckets[i].SetNext((ptr + 1));
        myBuckets[i].SetPrev((ptr - 1));
        myBuckets[newOrder[i]].SetValue(i);

        myBuckets[i].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    // Update value for target nodes
    for (int j = myNumSensors ; j < myNumElements ; ++j) {

        // reset initial status
        myBuckets[j].ResetNodeReached();
        myBuckets[j].ResetNodePermanent();
        myBuckets[j].UpdateCapacityLabel(0);

        myBuckets[j].SetID(j);
        myBuckets[j].SetNext((ptr + 1));
        myBuckets[j].SetPrev((ptr - 1));
        myBuckets[j].SetValue(myNumElements);

        myBuckets[j].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    myBuckets[myNumElements-1].SetNext(nullptr);
}

// clears out all information from a previous test instance
void NodeLinkedList::ClearOrderTimes() {
    orderedTargetFail.clear();
    orderedFail.clear();
    orderedCriticalFailTime.clear();
    orderedTargetCriticalFailTime.clear();
    orderedSensorCriticalFailTime.clear();
    orderedSensorFail.clear();

}

// Input the overall time interval, and calculate reliability for each unit of time
// myDSpec - number of sensors failed that cause network failure
// myDSpecProb - probability that the corresponding entry in myDSpec causes network failure
void NodeLinkedList::CalculateReliabilityWeibull(double t, double myShape, double myScale) {
/*    for (int i = 0; i <= t ; ++i) {
        cout << " Reliability : " << NetworkReliabilityWeibull(i, myNumSensors - 1, myDSpec, myDSpecProb, myShape, myScale) << endl;
//        NetworkReliabilityWeibull(i, myNumSensors - 1, myDSpec, myDSpecProb, myShape, myScale);

    }*/

    double i = 0;
    while (i < t) {
        cout << "Reliability " << myNumSensors << " t = " << i << " -- "
             << NetworkReliabilityWeibull(i, myNumSensors - 1, myDSpec, myDSpecProb, myShape, myScale) << endl;
        i = i + 0.2;
    }
}


// Input the overall time interval, and calculate variance for each unit of time
// myDSpec - number of sensors failed that cause network failure
// myDSpecProb - probability that the corresponding entry in myDSpec causes network failure
void NodeLinkedList::CalculateReliabilityVarianceWeibull(double t, int numReps, double myShape, double myScale) {
    for (int i = 0; i <= t ; ++i) {
        cout << "Variance t =  " << i << " -- " << NetworkReliabilityVarianceWeibull(i, myNumSensors - 1, myDSpec, myDSpecProb, numReps, myShape, myScale) << endl;
    }
}

// Calculates the stable network reliability in the presence of a delta TBM policy
void NodeLinkedList::CalculateStableMaintenanceReliabilityWeibull(double delta, double mShape, double mScale) {

//    cout << "Stable Network Reliability delta = " << delta << " : " << StableMaintenanceReliabilityWeibull(delta, myNumSensors - 1, myDSpec, myDSpecProb) << endl;
    cout << "MaintenanceRel " << myNumSensors << " delta = " << delta << " -- " <<  StableMaintenanceReliabilityWeibull(delta, myNumSensors - 1, myDSpec, myDSpecProb, mShape, mScale) << endl;
//    StableMaintenanceReliabilityWeibull(delta, myNumSensors - 1, myDSpec, myDSpecProb, mShape, mScale);

}

// Calculates the stable network reliability in the presence of a delta TBM policy
void NodeLinkedList::CalculateStableMaintenanceVarianceWeibull(double delta, int numReps, double mShape, double mScale) {

    cout << "MaintenanceVar " << myNumSensors << " delta = " << delta << " -- " <<  StableMaintenanceVarianceWeibull(delta, myNumSensors - 1, myDSpec, myDSpecProb, numReps, mShape, mScale) << endl;

}

// Calculates the cost of stable network reliability in the presence of a delta TBM policy
void NodeLinkedList::CalculateMaintenancePolicyCostWeibull(double delta, double costFixed, double costVar, double mShape,
                                                           double mScale) {

//    cout << "Stable Network Reliability Cost with delta " << delta << " : " << MaintenancePolicyCostWeibull(delta, costFixed, costVar, myNumSensors - 1) << endl;
    cout << MaintenancePolicyCostWeibull(delta, costFixed, costVar, myNumSensors - 1, mShape, mScale) << endl;
//    MaintenancePolicyCostWeibull(delta, costFixed, costVar, myNumSensors - 1, mShape, mScale);

}

// Calculates the cost of stable network reliability in the presence of a delta TBM policy
void NodeLinkedList::CalculateMaintenancePolicyCostWeibull(int numSens, double delta, double costFixed, double costVar, double mShape,
                                                           double mScale) {

//    cout << "Stable Network Reliability Cost with delta " << delta << " : " << MaintenancePolicyCostWeibull(delta, costFixed, costVar, myNumSensors - 1) << endl;
    cout << MaintenancePolicyCostWeibull(delta, costFixed, costVar, numSens, mShape, mScale) << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Generates a new RGG, with a new x,y location for every sensor
void NodeLinkedList::GenerateNewRGG_naive() {

    Node *ptr = myBuckets;
    //Sink Location
    ptr->SetNewLoc(0.5,0.5);
    ptr->ResetReachableNodes();

    ptr = ptr + 1;

//    random_device myDev;
    uniform_real_distribution<> myUnif(0,1);

    double myXNum;
    double myYNum;

    // generate new location for every sensor
    for (int i = 1; i < myNumSensors ; ++i) {
        myXNum = myUnif(myDevX);
        myYNum = myUnif(myDevY);

        ptr->ResetReachableNodes();
        ptr->SetNewLoc(myXNum, myYNum);
        ptr = ptr + 1;
    }

    // build forward star for sink node separately, as we do not have to examine target nodes
    int startNode = 0;
    for (int l = startNode + 1; l < myNumSensors ; ++l) {
        double dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[l].GetXLoc(), 2) + pow(myBuckets[startNode].GetYLoc() - myBuckets[l].GetYLoc(), 2));
        if (dist <= myCommunicationRange){
            myBuckets[startNode].AddReachableNode(myBuckets[l].GetNodeID());
        }
    }

    // rebuild forward star for every node
    for (int j = 1; j < myNumSensors ; ++j) {
        for (int i = (j+1); i <myNumSensors ; ++i) {
            double dist = sqrt(pow(myBuckets[j].GetXLoc() - myBuckets[i].GetXLoc(), 2) + pow(myBuckets[j].GetYLoc() - myBuckets[i].GetYLoc(), 2));
            if (dist <= myCommunicationRange){
                myBuckets[j].AddReachableNode(myBuckets[i].GetNodeID());
                myBuckets[i].AddReachableNode(myBuckets[j].GetNodeID());
            }
        }


        // Add Target nodes to forward star
        for (int k = myNumSensors; k < (myNumSensors + myNumTargets) ; ++k) {
            double dist = sqrt(pow(myBuckets[j].GetXLoc() - myBuckets[k].GetXLoc(), 2) + pow(myBuckets[j].GetYLoc() - myBuckets[k].GetYLoc(), 2));
            if (dist <= mySensorRange){
                myBuckets[j].AddReachableNode(myBuckets[k].GetNodeID());
            }
        }
    }
}

// TODO: NOTE - if the number of grid cells is larger than the total number of nodes in the network, this process will have issues!!!
void NodeLinkedList::GenerateNewRGG_Bucket(){

//    random_device myDev;
    uniform_real_distribution<> myUnif(0,1);

    double myXNum;
    double myYNum;
    int myTempBin;

    int myBucketsPerRow = (int) ceil(1 / myCommunicationRange);

    Node *ptr = myBuckets;
    //Sink Location
    ptr->SetNewLoc(0.5,0.5);
    // clear the forward start of the node, and nodes in bin zero
    ptr->ResetReachableNodes();
    ptr->ResetNodesInBucket();

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    for (int i = 1; i < myNumSensors ; ++i) {
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
        ptr = ptr + 1;

    }

    // Update value for target nodes
    for (int j = myNumSensors ; j < myNumElements ; ++j) {
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
        ptr = ptr + 1;
    }

    // start placing sensor in bins, starting with sink node
    myXNum = myBuckets[0].GetXLoc();
    myYNum = myBuckets[0].GetYLoc();
    // determine which grid number the new sensor location falls in
    myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

    myBuckets[myTempBin].AddNodeToBucket(0);

    ptr = myBuckets;
    ptr = ptr + 1;

    // generate new location for every sensor
    for (int i = 1; i < myNumSensors ; ++i) {
        myXNum = myUnif(myDevX);
        myYNum = myUnif(myDevY);

        ptr->SetNewLoc(myXNum, myYNum);

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow + 1);

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;
    }

    // place target nodes in new bin
    for (int k = myNumSensors; k < myNumElements ; ++k) {

        myXNum = ptr->GetXLoc() + 0.0001;
        myYNum = ptr->GetYLoc() + 0.0001;

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

        //cout << " Node " << ptr->GetNodeID() << " Added to bucket " << myTempBin << endl;

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;

    }

    int myXBucket;
    int myYBucket;
    double dist;
    vector<int> myNodesInCurBucket;

    int curBucket = myBucketsPerRow + 1;
    for (int l = 0; l < pow(myBucketsPerRow,2) ; ++l) {
        // determine which grid cell we are currently in
        myXBucket = curBucket % myBucketsPerRow;
        myYBucket = (curBucket / myBucketsPerRow) % myBucketsPerRow ;

        if (myXBucket == 0){
            myXBucket = myBucketsPerRow;
            if (curBucket / myBucketsPerRow == myBucketsPerRow){
                myYBucket = myBucketsPerRow - 1;
            }
            else {
                myYBucket = myYBucket - 1;
            }
        }

        if (myYBucket == 0){
            myYBucket = myBucketsPerRow;
        }

        // all within bin comparison
        for (int i = 0; i < (myBuckets[l].GetNumberInCell()) ; ++i) {

            int startNode = myBuckets[l].GetNodeInBucket(i);

            // if the starting node is the sink node, then we only need to worry about connecting to sensor nodes, which is based on communication range
            if (startNode == 0) {
                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);
                    // determine if end node is a sensor node or a target node
                    if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
            // if starting node is a sensor node
            else if (startNode > 0 && startNode < myNumSensors){

                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);

                    // start node is a sensor node, end node is sink node, check for reverse arc
                    if (endNode == 0){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sink comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[endNode].AddReachableNode(startNode);
                        }

                    }
                    else if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sensor comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                            myBuckets[endNode].AddReachableNode(startNode);
                        }
                    }
                    else {
                        // end node must be a target node
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= mySensorRange){
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
        }

        int myAdjacentCell;

        // look to the cell immediately above current cell
        int myTopBucket = myYBucket + 1;
        myAdjacentCell = myXBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow +1);
        if (l < myAdjacentCell && myTopBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately to the right
        int myRightBucket = myXBucket + 1;
        myAdjacentCell = myRightBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myRightBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately below
        int myBottomBucket = myYBucket - 1;
        myAdjacentCell = myXBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myBottomBucket < 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to left of current cell
        int myLeftBucket = myXBucket - 1;
        myAdjacentCell = myLeftBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myLeftBucket > 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top right
        myAdjacentCell = myRightBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom right
        myAdjacentCell = myRightBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top left
        myAdjacentCell = myLeftBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l< myAdjacentCell) && (myLeftBucket > 0 ) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom left
        myAdjacentCell = myLeftBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myLeftBucket > 0 ) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        curBucket = curBucket + 1;
    }
}

void NodeLinkedList::CompareAdjacentCells(int cell1, int cell2) {
    double dist;
    int startNode;
    int endNode;

    for (int i = 0; i < (myBuckets[cell1].GetNumberInCell()) ; ++i) {

        startNode = myBuckets[cell1].GetNodeInBucket(i);

        // if the starting node is the sink node, then we only need to worry about connecting to sensor nodes, which is based on communication range
        if (startNode == 0) {
            for (int j = 0; j < myBuckets[cell2].GetNumberInCell() ; ++j) {
                endNode = myBuckets[cell2].GetNodeInBucket(j);
                // determine if end node is a sensor node or a target node
                if (endNode < myNumSensors){
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    if (dist <= myCommunicationRange) {
                        myBuckets[startNode].AddReachableNode(endNode);
                    }
                }
            }
        }
        // if starting node is a sensor node
        else if (startNode > 0 && startNode < myNumSensors){

            for (int j = 0; j < myBuckets[cell2].GetNumberInCell() ; ++j) {
                endNode = myBuckets[cell2].GetNodeInBucket(j);

                if (endNode == 0){
                    // start node is a sensor node, end node is sink node, we check to add arc since nodes are in different grid cell
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    // sensor to sink comparison, we can arc from sink to sensor node
                    if (dist <= myCommunicationRange) {
                        myBuckets[endNode].AddReachableNode(startNode);
                    }
                }
                else if (endNode < myNumSensors){
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    // sensor to sensor comparison, we can add both arcs
                    if (dist <= myCommunicationRange) {
                        myBuckets[startNode].AddReachableNode(endNode);
                        myBuckets[endNode].AddReachableNode(startNode);
                    }
                }
                else {
                    // end node must be a target node
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    if (dist <= mySensorRange){
                        myBuckets[startNode].AddReachableNode(endNode);
                    }
                }
            }
        }
        // starting node is a target node, but we still need to check for sensor to target arcs
        else {
            for (int j = 0; j < myBuckets[cell2].GetNumberInCell() ; ++j) {
                endNode = myBuckets[cell2].GetNodeInBucket(j);

                if ((endNode == 0) || (endNode >= myNumSensors)){
                    // start node is target, end node is either a target or sink node, we skip
                }
                else{
                    // end node must be a sensor node
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    if (dist <= mySensorRange){
                        myBuckets[endNode].AddReachableNode(startNode);
                    }
                }
            }
        }
    }
}


// used to generate new fail times
void NodeLinkedList::GenerateNewFailOrder(vector<int> &nodeTime){

    // TODO change default random engine to a different generator
    //    default_random_engine myGen;
    //    mersenne_twister_engine myMersenne();

    nodeTime.clear();
    nodeTime.reserve(myNumElements - 1);

    for (int j = 0; j < (myNumSensors - 1) ; ++j) {
        nodeTime.push_back(j + 1);
    }

    for (int j = 0; j < (myNumSensors - 1) ; ++j) {
        // generate random number
        uniform_int_distribution<int> myDist(0, myNumSensors - j - 2);
        int myNum = myDist(myDevOrder);

        // store the last number that is not marked
        int tempFirst = nodeTime.at(myNum);
        // element to swap
        int tempLast = nodeTime.at(myNumSensors - j - 2);

        // swap elements
        nodeTime[myNum] = tempLast;
        nodeTime[myNumSensors - j - 2] = tempFirst;

    }

    // set fail time for sink
    nodeTime.insert(nodeTime.begin(), myNumSensors);

    // initial time for all targets
    for (int j = 0; j < myNumTargets ; ++j) {
        nodeTime.push_back(myNumSensors);
    }
}

// NOTE : newSize includes the sink node!!!! This outputs to the command line
void NodeLinkedList::SignatureRelation(int newSize, double t, double maxDelta, double costFixed, double costVar, double myShape, double myScale) {

    vector<int> newFailed;
    vector<double> newProb;

    int sizeDiff = myNumSensors - newSize;

    if (sizeDiff < 0) {
        cout << " The new network size is larger than initial network" << endl;
    }
    else {

        int firstFail = myDSpec.at(0);

        if ((firstFail - sizeDiff) <= 0){
            newFailed.push_back(0);
            newProb.push_back(myDSpecProb.at(0));
        }

        for (int i = 1; i < myDSpec.size(); ++i) {
            if ((myDSpec.at(i) - sizeDiff) <= 0){
                newProb[0] = newProb[0] + myDSpecProb.at(i);
            }
            else{
                newFailed.push_back(myDSpec.at(i) - sizeDiff);
                newProb.push_back(myDSpecProb.at(i));
            }
        }
    }

    double i = 1.0;
    while (i < maxDelta){
        cout << "MaintenanceRel " << newSize << " delta = " << i << " -- " << StableMaintenanceReliabilityWeibull(i, newSize - 1, newFailed, newProb, myShape, myScale) << " -- " << MaintenancePolicyCostWeibull(i, costFixed, costVar, newSize - 1, myShape, myScale) << endl;
        i = i + 0.1;
    }
}

// outputs the maintenance reliability and cost to a txt file
void NodeLinkedList::SignatureRelationOutput(int sizeRange, double t, double maxDelta, double costFixed, double costVar,
                                             double myShape, double myScale){

    ofstream myFile;
    myFile.open ("MaintReliability.txt");

    int currNum = myNumSensors - 1;
    while (currNum >= 500) {
        // + 1 to include the sink node
        vector<int> newFailed;
        vector<double> newProb;

        int sizeDiff = myNumSensors - currNum;

        if (sizeDiff < 0) {
            cout << " The new network size is larger than initial network" << endl;
        }
        else {

            int firstFail = myDSpec.at(0);

            if ((firstFail - sizeDiff) <= 0){
                newFailed.push_back(0);
                newProb.push_back(myDSpecProb.at(0));
            }

            for (int i = 1; i < myDSpec.size(); ++i) {
                if ((myDSpec.at(i) - sizeDiff) <= 0){
                    newProb[0] = newProb[0] + myDSpecProb.at(i);
                }
                else{
                    newFailed.push_back(myDSpec.at(i) - sizeDiff);
                    newProb.push_back(myDSpecProb.at(i));
                }
            }
        }

        double i = 1.0;
        while (i < maxDelta){
            myFile << "MaintenanceRel " << currNum << " delta = " << i << " -- " << StableMaintenanceReliabilityWeibull(i, currNum - 1, newFailed, newProb, myShape, myScale) << " -- " << MaintenancePolicyCostWeibull(i, costFixed, costVar, currNum - 1, myShape, myScale) << endl;
            i = i + 0.1;
        }

        currNum = currNum - 1;
    }

    myFile.close();
}

// outputs the maintenance reliability and cost to a txt file
void NodeLinkedList::SignatureRelationVariance(int sizeRange, int numReps, double t, double maxDelta, double costFixed, double costVar,
                                             double myShape, double myScale){

    ofstream myFile;
    myFile.open ("MaintReliabilityVar.txt");

    int currNum = myNumSensors;

    while (currNum >=500){
        // + 1 to include the sink node
        vector<int> newFailed;
        vector<double> newProb;

        int sizeDiff = myNumSensors - currNum;

        if (sizeDiff < 0) {
            cout << " The new network size is larger than initial network" << endl;
        }
        else {

            int firstFail = myDSpec.at(0);

            if ((firstFail - sizeDiff) <= 0){
                newFailed.push_back(0);
                newProb.push_back(myDSpecProb.at(0));
            }

            for (int i = 1; i < myDSpec.size(); ++i) {
                if ((myDSpec.at(i) - sizeDiff) <= 0){
                    newProb[0] = newProb[0] + myDSpecProb.at(i);
                }
                else{
                    newFailed.push_back(myDSpec.at(i) - sizeDiff);
                    newProb.push_back(myDSpecProb.at(i));
                }
            }
        }

        double i = 1.0;
        while (i < maxDelta){
            myFile << "MaintenanceVar " << currNum << " delta = " << i << " -- " << StableMaintenanceVarianceWeibull(i, currNum - 1, newFailed, newProb, numReps, myShape, myScale) << endl;
            i = i + 0.1;
        }

        currNum = currNum - 100;
    }

    myFile.close();
}

// something is wrong with how solutions are added to the efficient frontier, so this is not exactly correct
void NodeLinkedList::EfficientFrontier(double maxDelta, double costFixed, double costVar, double myShape, double myScale) {

    cout << " calculating Efficient Frontier " << endl;

    vector<double> policyReliability;
    vector<double> policyCost;
    vector<int> networkSize;
    vector<double> networkDelta;

    double tempRel;
    double tempCost;

    double i = 1.0;
    while (i < 10.0){
        tempRel  = StableMaintenanceReliabilityWeibull(i, myNumSensors - 1, myDSpec, myDSpecProb, myShape, myScale);
        tempCost = MaintenancePolicyCostWeibull(i, costFixed, costVar, myNumSensors - 1, myShape, myScale);

        networkSize.push_back(myNumSensors - 1);
        networkDelta.push_back(i);
        policyReliability.push_back(tempRel);
        policyCost.push_back(tempCost);

        i = i + 0.1;
    }


    int currNum = myNumSensors - 1;
    while (currNum >= 501) {
        // + 1 to include the sink node
        vector<int> newFailed;
        vector<double> newProb;

        int sizeDiff = myNumSensors - currNum;

        if (sizeDiff < 0) {
            cout << " The new network size is larger than initial network" << endl;
        }
        else {

            int firstFail = myDSpec.at(0);

            if ((firstFail - sizeDiff) <= 0){
                newFailed.push_back(0);
                newProb.push_back(myDSpecProb.at(0));
            }

            for (int j = 1; j < myDSpec.size(); ++j) {
                if ((myDSpec.at(j) - sizeDiff) <= 0){
                    newProb[0] = newProb[0] + myDSpecProb.at(j);
                }
                else{
                    newFailed.push_back(myDSpec.at(j) - sizeDiff);
                    newProb.push_back(myDSpecProb.at(j));
                }
            }
        }


        i = 1.0;
        while (i < maxDelta){
            tempRel = StableMaintenanceReliabilityWeibull(i, currNum - 1, newFailed, newProb, myShape, myScale);
            tempCost = MaintenancePolicyCostWeibull(i, costFixed, costVar, currNum - 1, myShape, myScale);

            networkSize.push_back(currNum - 1);
            networkDelta.push_back(i);
            policyReliability.push_back(tempRel);
            policyCost.push_back(tempCost);

            i = i + 0.1;
        }

        currNum = currNum - 1;
    }

    vector<double> efficientRel;
    vector<double> efficientCost;
    vector<int> efficientSize;
    vector<double> efficientDelta;

    for (int k = 0; k < policyCost.size() ; ++k) {
        double currRel = policyReliability.at(k);
        double currCost = policyCost.at(k);

        bool foundBetter = false;

        for (int j = 0; j < policyCost.size() ; ++j) {

            if (k != j){
                if ((policyReliability.at(j) > currRel) && (policyCost.at(j) < currCost)){
                    foundBetter = true;
                }
            }
        }

        if (!foundBetter){
            efficientSize.push_back(networkSize.at(k));
            efficientCost.push_back(currCost);
            efficientRel.push_back(currRel);
            efficientDelta.push_back(networkDelta.at(k));
        }
    }

    ofstream myFile;
    myFile.open ("MaintReliability.txt");

    for (int l = 0; l < policyCost.size() ; ++l) {
        myFile << "MaintenanceRel " << networkSize.at(l) << " delta = " << networkDelta.at(l) << " -- " << policyReliability.at(l) << " -- " << policyCost.at(l) << endl;
    }

    myFile << " Efficient Frontier Policies " << endl;
    for (int m = 0; m < efficientCost.size() ; ++m) {
        myFile << "MaintenanceRel " << efficientSize.at(m) << " delta = " << efficientDelta.at(m) << " -- " << efficientRel.at(m) << " -- " << efficientCost.at(m) << endl;
    }

    myFile.close();

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Uniformly locate target nodes over network
void NodeLinkedList::UniformTargetLoc() {

    double n = sqrt(myNumTargets);

    double dist = 1 / (n - 1);

    Node *ptr = myBuckets;
    ptr = ptr + myNumSensors;

    for (int i = 0; i < n ; ++i) {
        for (int j = 0; j < n ; ++j) {
            ptr->SetNewLoc(dist * (double) i, dist * (double) j);
            ptr = ptr + 1;
        }
    }
}

void NodeLinkedList::SetRanges(double comRange, double sensingRange) {
    myCommunicationRange = comRange;
    mySensorRange = sensingRange;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                         Read Network Data
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Read in Network Data from .txt file
// reads in a txt file where the format of the first row is:
// number of nodes(sensors)     |     number of arcs        |      number of targets
void ReadNetwork(char* filename, vector<Edge*> &edgeVect, vector<int> &nodeVect, int &numT){

    ifstream file;
    int n;          //number of nodes
    int m;          //number of edges
    int numTar;     //number of targets
    int f;          //variable to store edge's "from node"
    int t;          //variable to store edge's "to node"
    int bval;       //variable to store node's value

    file.open(filename);    //open network file

    file >> n;      //read number of nodes
    file >> m;      //read number of edge

    file >> numTar; //read number of targets
    numT = numTar;

    //read node data, populate supply/demand vector
    for (int i = 0; i < (n + numTar); i++){
        file >> bval;
        nodeVect.push_back(bval);
    }

    //read edge data, populate edge list
    for (int e = 0; e < m; e++){
        file >> f;
        file >> t;
       // Edge *newEdge(f,t);
        Edge *a = new Edge(f,t);
        edgeVect.push_back(a);

    }
    file.close();
}

// Read in Network Data from .txt file
// reads in a txt file where the format of the first row is:
// number of nodes     |     number of arcs
void ReadNetwork(char* filename, vector<Edge*> &edgeVect, vector<int> &nodeVect){

    ifstream file;
    int n;      //number of nodes
    int m;      //number of edges
    int f;      //variable to store edge's "from node"
    int t;      //variable to store edge's "to node"
    int bval;   //variable to store node's value

    file.open(filename);    //open network file

    file >> n;     //read number of nodes
    file >> m;     //read number of edge

    //read node data, populate supply/demand vector
    for (int i = 0; i < n ; i++){
        file >> bval;
        nodeVect.push_back(bval);
    }

    //read edge data, populate edge list
    for (int e = 0; e < m; e++){
        file >> f;
        file >> t;
        // Edge *newEdge(f,t);
        Edge *a = new Edge(f,t);
        edgeVect.push_back(a);
    }
    file.close();
}

// used to update the fail times for a new test instance. Assumes the remaining network topology is the same
// format of input file is
// number of nodes (including sink) | number of targets | number of instances
void UpdateTimes(char* filename, vector<int> &nodeVect, int &numIterations){

    ifstream file;
    int numInstances;   //number of instances in the file
    int n;              //number of nodes
    int numTar;         //number of targets
    int bval;           //variable to store node's value

    file.open(filename);    //open network file

    file >> n;      //read number of nodes
    file >> numTar; //read number of targets
    file >> numInstances;

    numIterations = numInstances;

    //read node data, populate supply/demand vector
    for (int i = 0; i < (n + numTar)*numInstances; i++){
        file >> bval;
        nodeVect.push_back(bval);
    }

    file.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                         Node Initialization Section
/////////////////////////////////////////////////////////////////////////////////////////////////////////

Node::Node() {
    myValue = 0;
    myID = 0;
    myMaxCapacityPath = 0;
    prev = nullptr;
    next = nullptr;
}

Node::Node(int inputID){
    myValue = 0;
    myID = inputID;
    myMaxCapacityPath = 0;
    prev = nullptr;
    next = nullptr;
}

// Node adjacency list
void Node::AddReachableNode(int nodeR){
    nodeDegree = nodeDegree + 1;
    forwardStar.push_back(nodeR);
}

int Node::GetReachableNode(int i){
    return forwardStar[i];
}

// Identifier for node
void Node::SetID(int temp){
    myID = temp;
}

int Node::GetNodeID(){
    return  myID;
}

// value associated with node
int Node::GetValue() {
    return myValue;
}

void Node::SetValue(int data){
    myValue = data;
}

void Node::UpdateCapacityLabel(int newCapacity){
    myMaxCapacityPath = newCapacity;
}

int Node::GetMaxCapacityPath(){
    return myMaxCapacityPath;
}

void Node::SetPrev(Node *newPrev){
    prev = newPrev;
}

void Node::SetNext(Node *newNext){
    next = newNext;
}

void Node::SetFirstItem(Node *firstElement){
    firstInBucket = firstElement;
}

void Node::UpdateReached(){
    isReached = true;
}

bool Node::IsReached(){
    return isReached;
}

void Node::MarkPermanent(){
    isPermanent = true;
}

bool Node::IsPermanent(){
    return isPermanent;
}

void Node::ResetNodeReached(){
    isReached = false;
}
void Node::ResetNodePermanent(){
    isPermanent = false;
}

int Node::GetDegree(){
    return nodeDegree;
}

// reset the forward star of a node
void Node::ResetReachableNodes(){
    forwardStar.clear();
    nodeDegree = 0;
}

Node *Node::GetPrev(){
    return prev;
}

Node *Node::GetNext(){
    return next;
}

Node *Node::GetFirstElement(){
    return firstInBucket;
}

void Node::SetNewLoc(double xLoc, double yLoc) {
    myXLoc = xLoc;
    myYLoc = yLoc;
}

double Node::GetXLoc() {
    return myXLoc;
}

double Node::GetYLoc() {
    return myYLoc;
}

void Node::ResetNodesInBucket() {
    nodesInBucket.clear();
    myNodesInBucket = 0;
}

void Node::AddNodeToBucket(int nodeNew) {
    nodesInBucket.push_back(nodeNew);
    myNodesInBucket = myNodesInBucket + 1;
}

int Node::GetNodeInBucket(int i) {
    return nodesInBucket[i];
}

int Node::GetNumberInCell() {
    return myNodesInBucket;
}

bool Node::IsCovered() {
    return isCovered;
}

void Node::SetCovered() {
    isCovered = true;
}

void Node::ResetCovered() {
    isCovered = false;
}

bool Node::IsFailed() {
    return isFailed;
}

void Node::SetFailed() {
    isFailed = true;
}

void Node::ResetFailed(){
    isFailed = false;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                         MC section
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SimulationNodeMC::SimulationNodeMC() {
    this->SetValue(0);
    this->SetID(0);
    this->UpdateFailCapacity(0);
}

NodeLinkedListMC::NodeLinkedListMC(int numItems, int numSens, int numTarg){
    firstItem = nullptr;
    lastItem = nullptr;
    myNumElements = numItems;   // total number of nodes
    myNumSensors = numSens;     // number of sensors (including sink)
    myNumTargets = numTarg;     // number of targets
    InitializeBuckets(numItems);
}


NodeLinkedListMC::NodeLinkedListMC(int numItems, int numSens, int numTarg, vector<double> failTimes){
    firstItem = nullptr;
    lastItem = nullptr;
    myNumElements = numItems;   // total number of nodes
    myNumSensors = numSens;     // number of sensors (including sink)
    myNumTargets = numTarg;     // number of targets
    InitializeBuckets(numItems, failTimes);
}

void NodeLinkedListMC::InitializeBuckets(int numItems){

    myBuckets = new SimulationNodeMC[numItems];

    SimulationNodeMC *ptr = myBuckets;
    // update all values for the sink node
    ptr->SetID(0);
    ptr->UpdateReached();
    ptr->MarkPermanent();
    ptr->SetFirstItem((ptr + 1));

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    // values for sensor nodes
    for (int i = 1; i < myNumSensors ; ++i) {
        myBuckets[i].SetID(i);
        myBuckets[i].SetNext((ptr + 1));
        myBuckets[i].SetPrev((ptr - 1));
        myBuckets[i].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    // values for target nodes
    for (int j = myNumSensors ; j < numItems ; ++j) {

        myBuckets[j].SetID(j);
        myBuckets[j].SetNext((ptr + 1));
        myBuckets[j].SetFirstItem(nullptr);
        ptr = ptr + 1;

    }
    myBuckets[numItems-1].SetNext(nullptr);
}

void NodeLinkedListMC::InitializeBuckets(int numItems,  vector<double> failTimes){

    double maxFail = *max_element(failTimes.begin(), failTimes.end());

//    max_element(failTimes.begin(), failTimes.end());

    myBuckets = new SimulationNodeMC[numItems];

    SimulationNodeMC *ptr = myBuckets;
    // update all values for the sink node
    ptr->SetID(0);
    ptr->UpdateFailCapacity(maxFail);
    ptr->UpdateReached();
    ptr->MarkPermanent();
    ptr->SetFailTime(failTimes[0]);
    ptr->SetFirstItem((ptr + 1));

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    // values for sensor nodes
    for (int i = 1; i < myNumSensors ; ++i) {
        myBuckets[i].SetID(i);
        myBuckets[i].SetNext((ptr + 1));
        myBuckets[i].SetPrev((ptr - 1));
        myBuckets[i].SetFailTime(failTimes[i]);

        myBuckets[i].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    // values for target nodes
    for (int j = myNumSensors ; j < numItems ; ++j) {

        myBuckets[j].SetID(j);
        myBuckets[j].SetNext((ptr + 1));
        myBuckets[j].SetPrev((ptr - 1));
        myBuckets[j].SetFailTime(maxFail);
        myBuckets[j].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }
    myBuckets[numItems-1].SetNext(nullptr);
}

void SimulationNodeMC::ResetNode(){
    isFailed = false;
    isCovered = false;
}

int SimulationNodeMC::GetNodeID(){
    return  myID;
}

bool SimulationNodeMC::IsCovered() {
    return isCovered;
}

void SimulationNodeMC::SetCovered() {
    isCovered = true;
}

void SimulationNodeMC::ResetCovered() {
    isCovered = false;
}

bool SimulationNodeMC::IsFailed() {
    return isFailed;
}

void SimulationNodeMC::SetFailed() {
    isFailed = true;
}

void SimulationNodeMC::ResetFailed() {
    isFailed = false;
}

void SimulationNodeMC::AddNodeToBucket(int nodeNew) {
    nodesInBucket.push_back(nodeNew);
    myNodesInBucket = myNodesInBucket + 1;

}

int SimulationNodeMC::GetNodeInBucket(int i) {
    return nodesInBucket[i];
}

int SimulationNodeMC::GetNumberInCell() {
    return myNodesInBucket;
}

double SimulationNodeMC::GetXLoc() {
    return myXLoc;
}

double SimulationNodeMC::GetYLoc() {
    return myYLoc;
}

// reset the forward star of a node
void SimulationNodeMC::ResetReachableNodes(){
    forwardStar.clear();
    nodeDegree = 0;
}

void SimulationNodeMC::ResetNodesInBucket() {
    nodesInBucket.clear();
    myNodesInBucket = 0;
}

double SimulationNodeMC::GetFailTime() {
    return myFailTime;
}

void SimulationNodeMC::SetNewLoc(double xLoc, double yLoc) {
    myXLoc = xLoc;
    myYLoc = yLoc;
}

void SimulationNodeMC::SetFailTime(double t) {
    myFailTime = t;
}

void SimulationNodeMC::UpdateFailCapacity(double newTime) {
    myCapacityFailTime = newTime;
}

double SimulationNodeMC::GetFailCapacity() {
    return myCapacityFailTime;
}

void NodeLinkedListMC::DijkstrasAlgorithm(vector<Edge *> edgeList, int numNodes, int numTar, double &timeDuration) {

    // sink node is already permanent
    int numPermanent = 1;

    // vector to store ID's of temporary nodes
    vector<int> tempNodes;
    for (int i = 1; i <= (numNodes-1) ; ++i) {
        tempNodes.push_back(i);
    }

    // start time record
    auto start = chrono::high_resolution_clock::now();

    // update all nodes connected to sink
    for (int i = 0; i < myBuckets[0].GetDegree(); ++i) {
        int updateNode1 = myBuckets[0].GetReachableNode(i);
        double updateLabel = min(myBuckets[0].GetFailCapacity(), myBuckets[updateNode1].GetFailTime());
        myBuckets[updateNode1].UpdateFailCapacity(updateLabel);
    }

    while (numPermanent < numNodes) {

        int tempNode;
        int tempIndex = -1;
        double currMax = 0;
        int permNode = 0;

        // find the node with the largest distance label that is not permanent
        for (int j = 0; j < tempNodes.size(); ++j) {
            tempNode = tempNodes.at(j);
            if (myBuckets[tempNode].GetFailCapacity() > currMax) {
                // store node ID and max label

                currMax = myBuckets[tempNode].GetFailCapacity();
                permNode = myBuckets[tempNode].GetNodeID();
                tempIndex = j;
            }
        }

        if (tempIndex == -1) {
            // node not connected - should be able to break here
        } else {

            // remove from temporary nodes vector
            tempNodes.erase(tempNodes.begin() + tempIndex);

            // mark node permanent
            myBuckets[permNode].MarkPermanent();

            this->addFailOrder(permNode);
            this->addFailTime(currMax);

            // stores the failure order and time at which Targets are disconnected
            if (myBuckets[permNode].GetNodeID() >= numNodes - numTar) {
                this->addTargetFailOrder(myBuckets[permNode].GetNodeID());
                this->addTargetFailTime(myBuckets[permNode].GetFailCapacity());

                myBuckets[permNode].SetCovered();
            }
            else {
                this->addSensorFailOrder(myBuckets[permNode].GetNodeID());
                this->addSensorFailTime(myBuckets[permNode].GetFailCapacity());
            }

            // iterate through edges
            for (int k = 0; k < myBuckets[permNode].GetDegree(); ++k) {

                int updateNode1 = myBuckets[permNode].GetReachableNode(k);
                if (!myBuckets[updateNode1].IsPermanent()) {
                    double updateLabel = min(myBuckets[permNode].GetFailCapacity(), myBuckets[updateNode1].GetFailTime());
                    myBuckets[updateNode1].UpdateFailCapacity(updateLabel);
                }
            }
        }
        numPermanent = numPermanent + 1;
    }

    // If a target is not connected to the sink then it will never be marked permanent, so it will never be added to the
    // order of target failures. This for loop checks for such a case, and adds the targets to ths list with a fail time of zero
    if (orderedTargetFail.size() < numTar){
        SimulationNodeMC *ptr = myBuckets;
        ptr = ptr + (numNodes - numTar);
        for (int i = 0; i < numTar ; ++i) {
            if (!ptr->IsPermanent()){
                this->addTargetFailOrder(ptr->GetNodeID());
                this->addTargetFailTime(ptr->GetFailCapacity());
            }
            ptr = ptr + 1;
        }
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;
    timeDuration = elapsed.count();
}

// used to generate new fail times
void NodeLinkedListMC::GenerateNewFailTime(vector<double> &nodeTime){

    nodeTime.clear();
    nodeTime.reserve(myNumElements - 1);

    weibull_distribution<double> myWeibull(1.5,10.0);

    for (int i = 0; i < myNumSensors - 1 ; ++i) {
        double myTime = myWeibull(myDevTime);
        nodeTime.push_back(myTime);
    }

    // set fail time for sink
    double maxFail = *max_element(nodeTime.begin(), nodeTime.end());
    nodeTime.insert(nodeTime.begin(), maxFail);

    // initial time for all targets
    for (int j = 0; j < myNumTargets ; ++j) {
        nodeTime.push_back(maxFail);
    }
}

// used to generate new fail times
void NodeLinkedListMC::GenerateNewFailTimeMaintenance(vector<double> &nodeTime, double t, double curTime){

    weibull_distribution<double> myWeibull(1.5,10.0);

    for (int i = 1; i < myNumSensors - 1 ; ++i) {

        if (myBuckets[i].IsFailed()){
            double myTime = myWeibull(myDevTime);
            nodeTime[i] = myTime + curTime;
        }
    }


    nodeTime[0] = t;

    // initial time for all targets
    for (int j = myNumSensors; j < myNumElements; ++j) {
        nodeTime[j] = t;
    }
}

void NodeLinkedListMC::UpdateFailTimes(vector<double> newTimes){

    orderedCriticalFailTime.clear();
    orderedTargetCriticalFailTime.clear();
    orderedSensorCriticalFailTime.clear();

    this->ClearOrderTimes();

    double maxFail = *max_element(newTimes.begin(), newTimes.end());

    SimulationNodeMC *ptr = myBuckets;
    // update all values for the sink node
    ptr->SetID(0);
    ptr->UpdateFailCapacity(maxFail);
    ptr->UpdateReached();
    ptr->MarkPermanent();
    ptr->SetFailTime(maxFail);
    ptr->SetFirstItem((ptr + 1));

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    for (int i = 1; i < myNumSensors ; ++i) {

        // reset initial status
        myBuckets[i].ResetNodeReached();
        myBuckets[i].ResetNodePermanent();
        myBuckets[i].UpdateFailCapacity(0);

        myBuckets[i].ResetCovered();

        // re-initialize
        myBuckets[i].SetNext((ptr + 1));
        myBuckets[i].SetPrev((ptr - 1));
        myBuckets[i].SetFailTime(newTimes[i]);

        myBuckets[i].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    // Update value for target nodes
    for (int j = myNumSensors ; j < myNumElements ; ++j) {

        // reset initial status
        myBuckets[j].ResetNodeReached();
        myBuckets[j].ResetNodePermanent();
        myBuckets[j].UpdateFailCapacity(0);

        myBuckets[j].ResetCovered();

        myBuckets[j].SetID(j);
        myBuckets[j].SetNext((ptr + 1));
        myBuckets[j].SetPrev((ptr - 1));
        myBuckets[j].SetFailTime(newTimes[j]);

        myBuckets[j].SetFirstItem(nullptr);

        ptr = ptr + 1;

    }

    myBuckets[myNumElements-1].SetNext(nullptr);
}

void NodeLinkedListMC::addFailTime(double failTime){
    orderedCriticalFailTime.push_back(failTime);
}

// vector that stores the ordered failed times for targets
void NodeLinkedListMC::addTargetFailTime(double failTime){
    orderedTargetCriticalFailTime.push_back(failTime);
}

// vector that stores the ordered fail time for sensors
void NodeLinkedListMC::addSensorFailTime(double failTime){
    orderedSensorCriticalFailTime.push_back(failTime);
}

// NOTE: the list is in reverse order (largest to smallest)
double NodeLinkedListMC::GetTargetFailTime(int failNum){
    return orderedTargetCriticalFailTime.at(orderedTargetCriticalFailTime.size() - failNum);
}

void NodeLinkedListMC::addNetworkFailTime(double failTime) {
    networkFailTimes.push_back(failTime);
}

void NodeLinkedListMC::CalculateReliability(double t, int numReps) {

    double i = 0;
    while (i < t) {

        int numFailed = 0;

        for (int j = 0; j < networkFailTimes.size(); ++j) {

            if (networkFailTimes.at(j) <= i){
                numFailed = numFailed + 1;
            }

        }

        double rel = 1.0 - ((double) numFailed / (double) numReps);

        cout << "Reliability " << myNumSensors << " t = " << i << " -- " << rel << endl;
        i = i + 0.2;
    }
}

// Uniformly locate target nodes over network
void NodeLinkedListMC::UniformTargetLoc() {

    double n = sqrt(myNumTargets);

    double dist = 1 / (n - 1);

    SimulationNodeMC *ptr = myBuckets;
    ptr = ptr + myNumSensors;

    for (int i = 0; i < n ; ++i) {
        for (int j = 0; j < n ; ++j) {
            ptr->SetNewLoc(dist * (double) i, dist * (double) j);
            ptr = ptr + 1;
        }
    }
}

void NodeLinkedListMC::GenerateNewRGG_Bucket(){

//    random_device myDev;
    uniform_real_distribution<> myUnif(0,1);

    double myXNum;
    double myYNum;
    int myTempBin;

    int myBucketsPerRow = (int) ceil(1 / myCommunicationRange);

    SimulationNodeMC *ptr = myBuckets;
    //Sink Location
    ptr->SetNewLoc(0.5,0.5);
    // clear the forward start of the node, and nodes in bin zero
    ptr->ResetReachableNodes();
    ptr->ResetNodesInBucket();

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    for (int i = 1; i < myNumSensors ; ++i) {
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
        ptr->ResetFailed();
        ptr->ResetCovered();
        ptr = ptr + 1;

    }

    // Update value for target nodes
    for (int j = myNumSensors ; j < myNumElements ; ++j) {
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
        ptr->ResetFailed();
        ptr->ResetCovered();
        ptr = ptr + 1;
    }

    // start placing sensor in bins, starting with sink node
    myXNum = myBuckets[0].GetXLoc();
    myYNum = myBuckets[0].GetYLoc();
    // determine which grid number the new sensor location falls in
    myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

    myBuckets[myTempBin].AddNodeToBucket(0);

    ptr = myBuckets;
    ptr = ptr + 1;

    // generate new location for every sensor
    for (int i = 1; i < myNumSensors ; ++i) {
        myXNum = myUnif(myDevX);
        myYNum = myUnif(myDevY);

        ptr->SetNewLoc(myXNum, myYNum);

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow + 1);

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;
    }

    // place target nodes in new bin
    for (int k = myNumSensors; k < myNumElements ; ++k) {

        myXNum = ptr->GetXLoc() + 0.0001;
        myYNum = ptr->GetYLoc() + 0.0001;

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

        //cout << " Node " << ptr->GetNodeID() << " Added to bucket " << myTempBin << endl;

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;

    }

    int myXBucket;
    int myYBucket;
    double dist;
    vector<int> myNodesInCurBucket;

    int curBucket = myBucketsPerRow + 1;
    for (int l = 0; l < pow(myBucketsPerRow,2) ; ++l) {
        // determine which grid cell we are currently in
        myXBucket = curBucket % myBucketsPerRow;
        myYBucket = (curBucket / myBucketsPerRow) % myBucketsPerRow ;

        if (myXBucket == 0){
            myXBucket = myBucketsPerRow;
            if (curBucket / myBucketsPerRow == myBucketsPerRow){
                myYBucket = myBucketsPerRow - 1;
            }
            else {
                myYBucket = myYBucket - 1;
            }
        }

        if (myYBucket == 0){
            myYBucket = myBucketsPerRow;
        }

        // all within bin comparison
        for (int i = 0; i < (myBuckets[l].GetNumberInCell()) ; ++i) {

            int startNode = myBuckets[l].GetNodeInBucket(i);

            // if the starting node is the sink node, then we only need to worry about connecting to sensor nodes, which is based on communication range
            if (startNode == 0) {
                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);
                    // determine if end node is a sensor node or a target node
                    if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
                // if starting node is a sensor node
            else if (startNode > 0 && startNode < myNumSensors){

                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);

                    // start node is a sensor node, end node is sink node, check for reverse arc
                    if (endNode == 0){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sink comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[endNode].AddReachableNode(startNode);
                        }

                    }
                    else if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sensor comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                            myBuckets[endNode].AddReachableNode(startNode);
                        }
                    }
                    else {
                        // end node must be a target node
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= mySensorRange){
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
        }

        int myAdjacentCell;

        // look to the cell immediately above current cell
        int myTopBucket = myYBucket + 1;
        myAdjacentCell = myXBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow +1);
        if (l < myAdjacentCell && myTopBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately to the right
        int myRightBucket = myXBucket + 1;
        myAdjacentCell = myRightBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myRightBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately below
        int myBottomBucket = myYBucket - 1;
        myAdjacentCell = myXBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myBottomBucket < 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to left of current cell
        int myLeftBucket = myXBucket - 1;
        myAdjacentCell = myLeftBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myLeftBucket > 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top right
        myAdjacentCell = myRightBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom right
        myAdjacentCell = myRightBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top left
        myAdjacentCell = myLeftBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l< myAdjacentCell) && (myLeftBucket > 0 ) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom left
        myAdjacentCell = myLeftBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myLeftBucket > 0 ) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        curBucket = curBucket + 1;
    }
}

void NodeLinkedListMC::CompareAdjacentCells(int cell1, int cell2) {
    double dist;
    int startNode;
    int endNode;

    for (int i = 0; i < (myBuckets[cell1].GetNumberInCell()) ; ++i) {

        startNode = myBuckets[cell1].GetNodeInBucket(i);

        // if the starting node is the sink node, then we only need to worry about connecting to sensor nodes, which is based on communication range
        if (startNode == 0) {
            for (int j = 0; j < myBuckets[cell2].GetNumberInCell() ; ++j) {
                endNode = myBuckets[cell2].GetNodeInBucket(j);
                // determine if end node is a sensor node or a target node
                if (endNode < myNumSensors){
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    if (dist <= myCommunicationRange) {
                        myBuckets[startNode].AddReachableNode(endNode);
                    }
                }
            }
        }
            // if starting node is a sensor node
        else if (startNode > 0 && startNode < myNumSensors){

            for (int j = 0; j < myBuckets[cell2].GetNumberInCell() ; ++j) {
                endNode = myBuckets[cell2].GetNodeInBucket(j);

                if (endNode == 0){
                    // start node is a sensor node, end node is sink node, we check to add arc since nodes are in different grid cell
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    // sensor to sink comparison, we can arc from sink to sensor node
                    if (dist <= myCommunicationRange) {
                        myBuckets[endNode].AddReachableNode(startNode);
                    }
                }
                else if (endNode < myNumSensors){
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    // sensor to sensor comparison, we can add both arcs
                    if (dist <= myCommunicationRange) {
                        myBuckets[startNode].AddReachableNode(endNode);
                        myBuckets[endNode].AddReachableNode(startNode);
                    }
                }
                else {
                    // end node must be a target node
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    if (dist <= mySensorRange){
                        myBuckets[startNode].AddReachableNode(endNode);
                    }
                }
            }
        }
            // starting node is a target node, but we still need to check for sensor to target arcs
        else {
            for (int j = 0; j < myBuckets[cell2].GetNumberInCell() ; ++j) {
                endNode = myBuckets[cell2].GetNodeInBucket(j);

                if ((endNode == 0) || (endNode >= myNumSensors)){
                    // start node is target, end node is either a target or sink node, we skip
                }
                else{
                    // end node must be a sensor node
                    dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                    if (dist <= mySensorRange){
                        myBuckets[endNode].AddReachableNode(startNode);
                    }
                }
            }
        }
    }
}

void NodeLinkedListMC::UpdateRGG(double time) {
    for (int i = 0; i < myNumSensors ; ++i) {

        if (myBuckets[i].GetFailTime() < time){
            myBuckets[i].SetFailed();
//            while (myBuckets[i].GetDegree() > 0 ){
//                myBuckets[i].UpdateForwardStar(0);
//            }
        }
        else{

            int k = 0;
            while (k < myBuckets[i].GetDegree()){

                int tempNode = myBuckets[i].GetReachableNode(k);
                if (myBuckets[tempNode].GetFailTime() < time){
                    myBuckets[i].UpdateForwardStar(k);
                    myBuckets[tempNode].SetFailed();
                }

                k = k + 1;
            }
        }
    }
}

void SimulationNodeMC::UpdateForwardStar(int removeNum) {
    forwardStar.erase(forwardStar.begin() + removeNum);
    nodeDegree = nodeDegree - 1;
}

void NodeLinkedListMC::GenerateNewRGG_BucketMaint(double costFixed, double costVar, double &policyCost){

    uniform_real_distribution<> myUnif(0,1);

    double myXNum;
    double myYNum;
    int myTempBin;
    int numReplaced = 0;
    bool maintOccur = false;
    double maintCost = 0;

    int myBucketsPerRow = (int) ceil(1 / myCommunicationRange);

    SimulationNodeMC *ptr = myBuckets;
    //Sink Location
    ptr->SetNewLoc(0.5,0.5);
    // clear the forward start of the node, and nodes in bin zero
    ptr->ResetReachableNodes();
    ptr->ResetNodesInBucket();

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    for (int i = 1; i < myNumSensors ; ++i) {
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
//        ptr->ResetNode();
        ptr = ptr + 1;

    }

    // Update value for target nodes
    for (int j = myNumSensors ; j < myNumElements ; ++j) {
//        ptr->ResetNode();
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
        ptr = ptr + 1;
    }

    // start placing sensor in bins, starting with sink node
    myXNum = myBuckets[0].GetXLoc();
    myYNum = myBuckets[0].GetYLoc();
    // determine which grid number the new sensor location falls in
    myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

    myBuckets[myTempBin].AddNodeToBucket(0);

    ptr = myBuckets;
    ptr = ptr + 1;

    // generate new location for every sensor
    for (int i = 1; i < myNumSensors ; ++i) {

        myXNum = myBuckets[i].GetXLoc();
        myYNum = myBuckets[i].GetYLoc();

        if (myBuckets[i].IsFailed()){
 //           cout << " New location " << endl;
            myXNum = myUnif(myDevX);
            myYNum = myUnif(myDevY);
//            cout << " Node " << myBuckets[i].GetNodeID() <<  " New x : " << myXNum << " --- new y : " << myYNum << endl;

            ptr->ResetFailed();
            ptr->SetNewLoc(myXNum, myYNum);

            numReplaced = numReplaced + 1;
            maintOccur = true;
        }

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow + 1);

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;
    }

    // place target nodes in new bin
    for (int k = myNumSensors; k < myNumElements ; ++k) {

        myXNum = ptr->GetXLoc() + 0.0001;
        myYNum = ptr->GetYLoc() + 0.0001;

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

//        cout << " Node " << ptr->GetNodeID() << " -- XLoc : " << myXNum << " -- YLoc : " << myYNum << endl;

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;

    }


    // Calculate cost of replacing sensors
    if (numReplaced > 0){
        maintCost = ((double) numReplaced) * (double) costVar + costFixed;
        policyCost = maintCost;

    }
    else{
        policyCost = 0.0;
    }


    // reset failed nodes and covered nodes
    ptr = myBuckets;
    ptr = ptr + 1;
    for (int m = 1; m < myNumElements; ++m) {
        ptr->ResetFailed();
        ptr->ResetCovered();
        ptr = ptr + 1;
    }

    int myXBucket;
    int myYBucket;
    double dist;
    vector<int> myNodesInCurBucket;

    int curBucket = myBucketsPerRow + 1;
    for (int l = 0; l < pow(myBucketsPerRow,2) ; ++l) {
        // determine which grid cell we are currently in
        myXBucket = curBucket % myBucketsPerRow;
        myYBucket = (curBucket / myBucketsPerRow) % myBucketsPerRow ;

        if (myXBucket == 0){
            myXBucket = myBucketsPerRow;
            if (curBucket / myBucketsPerRow == myBucketsPerRow){
                myYBucket = myBucketsPerRow - 1;
            }
            else {
                myYBucket = myYBucket - 1;
            }
        }

        if (myYBucket == 0){
            myYBucket = myBucketsPerRow;
        }

        // all within bin comparison
        for (int i = 0; i < (myBuckets[l].GetNumberInCell()) ; ++i) {

            int startNode = myBuckets[l].GetNodeInBucket(i);

            // if the starting node is the sink node, then we only need to worry about connecting to sensor nodes, which is based on communication range
            if (startNode == 0) {
                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);
                    // determine if end node is a sensor node or a target node
                    if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
                // if starting node is a sensor node
            else if (startNode > 0 && startNode < myNumSensors){

                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);

                    // start node is a sensor node, end node is sink node, check for reverse arc
                    if (endNode == 0){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sink comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[endNode].AddReachableNode(startNode);
                        }

                    }
                    else if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sensor comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                            myBuckets[endNode].AddReachableNode(startNode);
                        }
                    }
                    else {
                        // end node must be a target node
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= mySensorRange){
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
        }

        int myAdjacentCell;

        // look to the cell immediately above current cell
        int myTopBucket = myYBucket + 1;
        myAdjacentCell = myXBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow +1);
        if (l < myAdjacentCell && myTopBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately to the right
        int myRightBucket = myXBucket + 1;
        myAdjacentCell = myRightBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myRightBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately below
        int myBottomBucket = myYBucket - 1;
        myAdjacentCell = myXBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myBottomBucket < 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to left of current cell
        int myLeftBucket = myXBucket - 1;
        myAdjacentCell = myLeftBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myLeftBucket > 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top right
        myAdjacentCell = myRightBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom right
        myAdjacentCell = myRightBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top left
        myAdjacentCell = myLeftBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l< myAdjacentCell) && (myLeftBucket > 0 ) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom left
        myAdjacentCell = myLeftBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myLeftBucket > 0 ) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        curBucket = curBucket + 1;
    }
}

void NodeLinkedListMC::GenerateNewRGG_BucketMaint(){

    uniform_real_distribution<> myUnif(0,1);

    double myXNum;
    double myYNum;
    int myTempBin;

    int myBucketsPerRow = (int) ceil(1 / myCommunicationRange);

    SimulationNodeMC *ptr = myBuckets;
    //Sink Location
    ptr->SetNewLoc(0.5,0.5);
    // clear the forward start of the node, and nodes in bin zero
    ptr->ResetReachableNodes();
    ptr->ResetNodesInBucket();

    ptr = ptr + 1;

    // sets the pointers for a node to be the node immediately prior and immediately after the node
    for (int i = 1; i < myNumSensors ; ++i) {
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
//        ptr->ResetNode();
        ptr = ptr + 1;

    }

    // Update value for target nodes
    for (int j = myNumSensors ; j < myNumElements ; ++j) {
 //       ptr->ResetNode();
        ptr->ResetReachableNodes();
        ptr->ResetNodesInBucket();
        ptr = ptr + 1;
    }

    // start placing sensor in bins, starting with sink node
    myXNum = myBuckets[0].GetXLoc();
    myYNum = myBuckets[0].GetYLoc();
    // determine which grid number the new sensor location falls in
    myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

    myBuckets[myTempBin].AddNodeToBucket(0);

    ptr = myBuckets;
    ptr = ptr + 1;

    // generate new location for every sensor
    for (int i = 1; i < myNumSensors ; ++i) {

        myXNum = myBuckets[i].GetXLoc();
        myYNum = myBuckets[i].GetYLoc();

        if (myBuckets[i].IsFailed()){
            myXNum = myUnif(myDevX);
            myYNum = myUnif(myDevY);
            ptr->ResetFailed();
            ptr->SetNewLoc(myXNum, myYNum);

        }

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow + 1);

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;
    }

    // place target nodes in new bin
    for (int k = myNumSensors; k < myNumElements ; ++k) {

        myXNum = ptr->GetXLoc() + 0.0001;
        myYNum = ptr->GetYLoc() + 0.0001;

        // determine which grid number the new sensor location falls in
        myTempBin = ((int) ceil(myXNum / myCommunicationRange) + (int) (ceil(myYNum / myCommunicationRange))*myBucketsPerRow) - (myBucketsPerRow +1);

        //cout << " Node " << ptr->GetNodeID() << " Added to bucket " << myTempBin << endl;

        // move sensor into a bin for this new location
        myBuckets[myTempBin].AddNodeToBucket(ptr->GetNodeID());

        ptr = ptr + 1;

    }

    // reset failed nodes and covered nodes
    ptr = myBuckets;
    ptr = ptr + 1;
    for (int m = 1; m < myNumElements; ++m) {
        ptr->ResetFailed();
        ptr->ResetCovered();
        ptr = ptr + 1;
    }

    int myXBucket;
    int myYBucket;
    double dist;
    vector<int> myNodesInCurBucket;

    int curBucket = myBucketsPerRow + 1;
    for (int l = 0; l < pow(myBucketsPerRow,2) ; ++l) {
        // determine which grid cell we are currently in
        myXBucket = curBucket % myBucketsPerRow;
        myYBucket = (curBucket / myBucketsPerRow) % myBucketsPerRow ;

        if (myXBucket == 0){
            myXBucket = myBucketsPerRow;
            if (curBucket / myBucketsPerRow == myBucketsPerRow){
                myYBucket = myBucketsPerRow - 1;
            }
            else {
                myYBucket = myYBucket - 1;
            }
        }

        if (myYBucket == 0){
            myYBucket = myBucketsPerRow;
        }

        // all within bin comparison
        for (int i = 0; i < (myBuckets[l].GetNumberInCell()) ; ++i) {

            int startNode = myBuckets[l].GetNodeInBucket(i);

            // if the starting node is the sink node, then we only need to worry about connecting to sensor nodes, which is based on communication range
            if (startNode == 0) {
                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);
                    // determine if end node is a sensor node or a target node
                    if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
                // if starting node is a sensor node
            else if (startNode > 0 && startNode < myNumSensors){

                for (int j = (i+1); j < myBuckets[l].GetNumberInCell() ; ++j) {
                    int endNode = myBuckets[l].GetNodeInBucket(j);

                    // start node is a sensor node, end node is sink node, check for reverse arc
                    if (endNode == 0){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sink comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[endNode].AddReachableNode(startNode);
                        }

                    }
                    else if (endNode < myNumSensors){
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        // sensor to sensor comparison, we can add both arcs
                        if (dist <= myCommunicationRange) {
                            myBuckets[startNode].AddReachableNode(endNode);
                            myBuckets[endNode].AddReachableNode(startNode);
                        }
                    }
                    else {
                        // end node must be a target node
                        dist = sqrt(pow(myBuckets[startNode].GetXLoc() - myBuckets[endNode].GetXLoc(), 2) +
                                    pow(myBuckets[startNode].GetYLoc() - myBuckets[endNode].GetYLoc(), 2));

                        if (dist <= mySensorRange){
                            myBuckets[startNode].AddReachableNode(endNode);
                        }
                    }
                }
            }
        }

        int myAdjacentCell;

        // look to the cell immediately above current cell
        int myTopBucket = myYBucket + 1;
        myAdjacentCell = myXBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow +1);
        if (l < myAdjacentCell && myTopBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately to the right
        int myRightBucket = myXBucket + 1;
        myAdjacentCell = myRightBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myRightBucket <= myBucketsPerRow){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell immediately below
        int myBottomBucket = myYBucket - 1;
        myAdjacentCell = myXBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myBottomBucket < 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to left of current cell
        int myLeftBucket = myXBucket - 1;
        myAdjacentCell = myLeftBucket + (myYBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if (l < myAdjacentCell && myLeftBucket > 0){
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top right
        myAdjacentCell = myRightBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom right
        myAdjacentCell = myRightBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myRightBucket <= myBucketsPerRow) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the top left
        myAdjacentCell = myLeftBucket + (myTopBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l< myAdjacentCell) && (myLeftBucket > 0 ) && (myTopBucket <= myBucketsPerRow)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        // look to the cell to the bottom left
        myAdjacentCell = myLeftBucket + (myBottomBucket * myBucketsPerRow) - (myBucketsPerRow + 1);
        if ((l < myAdjacentCell) && (myLeftBucket > 0 ) && (myBottomBucket > 0)) {
            this->CompareAdjacentCells(l, myAdjacentCell);
        }

        curBucket = curBucket + 1;
    }
}

double NodeLinkedListMC::TargetCoveragePercentage() {

    double numCovered = 0;

    for (int i = myNumSensors; i < myNumTargets ; ++i) {
        if (myBuckets[i].IsCovered()){
            numCovered = numCovered + 1;
        }
    }

    double percentCovered = (double) numCovered / (double) myNumTargets;

    cout << " numCovered " << numCovered << " -- myNumTargets " << myNumTargets << endl;

    return percentCovered;

}

double NodeLinkedListMC::BreadthFirstSearchReturn(int numNodes, int numTar, double &timeDuration) {

    // sink node is already permanent
    int numPermanent = 1;

    vector<int> myList;
    int numCovered = 0;

    // start time record
    auto start = chrono::high_resolution_clock::now();

    myBuckets[0].UpdateReached();

    // update all nodes connected to sink
    for (int i = 0; i < myBuckets[0].GetDegree(); ++i) {
        int updateNode1 = myBuckets[0].GetReachableNode(i);


        if (!myBuckets[updateNode1].IsFailed()){
            myList.push_back(updateNode1);

            if (myBuckets[updateNode1].IsFailed()){
                cout << " This Node Has Failed!!!! " << endl;
            }

        }

        double updateLabel = min(myBuckets[0].GetFailCapacity(), myBuckets[updateNode1].GetFailTime());
        myBuckets[updateNode1].UpdateFailCapacity(updateLabel);



        myBuckets[updateNode1].UpdateReached();
    }

    while (myList.size() > 0) {

        double currMax = 0;
        int permNode = 0;

        int tempIndex = myList.at(0);
        myList.erase(myList.begin());

        permNode = tempIndex;

        if (myBuckets[permNode].IsFailed()){
            cout << " This Node Has Failed!!!! " << endl;
        }

        if (tempIndex == -1) {
            // node not connected - should be able to break here
        } else {

            // mark node permanent
            myBuckets[permNode].MarkPermanent();
            myBuckets[permNode].UpdateReached();

            this->addFailOrder(permNode);
            this->addFailTime(currMax);

            // stores the failure order and time at which Targets are disconnected
            if (myBuckets[permNode].GetNodeID() >= numNodes - numTar) {
                this->addTargetFailOrder(myBuckets[permNode].GetNodeID());
                this->addTargetFailTime(myBuckets[permNode].GetFailCapacity());

                myBuckets[permNode].SetCovered();
                numCovered = numCovered + 1;

            }
            else {
                this->addSensorFailOrder(myBuckets[permNode].GetNodeID());
                this->addSensorFailTime(myBuckets[permNode].GetFailCapacity());
            }

            // iterate through edges
            for (int k = 0; k < myBuckets[permNode].GetDegree(); ++k) {

                int updateNode1 = myBuckets[permNode].GetReachableNode(k);

                if (!myBuckets[updateNode1].IsReached() && !myBuckets[updateNode1].IsFailed()){
                    myList.push_back(updateNode1);
                    myBuckets[updateNode1].UpdateReached();
                }

                if (!myBuckets[updateNode1].IsPermanent()) {
                    double updateLabel = min(myBuckets[permNode].GetFailCapacity(), myBuckets[updateNode1].GetFailTime());
                    myBuckets[updateNode1].UpdateFailCapacity(updateLabel);
                }
            }
        }
        numPermanent = numPermanent + 1;
    }

    // If a target is not connected to the sink then it will never be marked permanent, so it will never be added to the
    // order of target failures. This for loop checks for such a case, and adds the targets to ths list with a fail time of zero
    if (orderedTargetFail.size() < numTar){
        SimulationNodeMC *ptr = myBuckets;
        ptr = ptr + (numNodes - numTar);
        for (int i = 0; i < numTar ; ++i) {
            if (!ptr->IsPermanent()){
                this->addTargetFailOrder(ptr->GetNodeID());
                this->addTargetFailTime(ptr->GetFailCapacity());
            }
            ptr = ptr + 1;
        }
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;
    timeDuration = elapsed.count();

    double percentCovered = (double) numCovered / (double) myNumTargets;

    return percentCovered;

}

void NodeLinkedListMC::BreadthFirstSearch(int numNodes, int numTar, double &timeDuration) {

    // sink node is already permanent
    int numPermanent = 1;

    vector<int> myList;

    // start time record
    auto start = chrono::high_resolution_clock::now();

    myBuckets[0].UpdateReached();

    // update all nodes connected to sink
    for (int i = 0; i < myBuckets[0].GetDegree(); ++i) {
        int updateNode1 = myBuckets[0].GetReachableNode(i);

        myList.push_back(updateNode1);

        double updateLabel = min(myBuckets[0].GetFailCapacity(), myBuckets[updateNode1].GetFailTime());
        myBuckets[updateNode1].UpdateFailCapacity(updateLabel);

        myBuckets[updateNode1].UpdateReached();
    }

    while (myList.size() > 0) {

        double currMax = 0;
        int permNode = 0;

        int tempIndex = myList.at(0);
        myList.erase(myList.begin());

        permNode = tempIndex;

        if (tempIndex == -1) {
            // node not connected - should be able to break here
        } else {

            // mark node permanent
            myBuckets[permNode].MarkPermanent();
            myBuckets[permNode].UpdateReached();

            this->addFailOrder(permNode);
            this->addFailTime(currMax);

            // stores the failure order and time at which Targets are disconnected
            if (myBuckets[permNode].GetNodeID() >= numNodes - numTar) {
                this->addTargetFailOrder(myBuckets[permNode].GetNodeID());
                this->addTargetFailTime(myBuckets[permNode].GetFailCapacity());

                myBuckets[permNode].SetCovered();

            }
            else {
                this->addSensorFailOrder(myBuckets[permNode].GetNodeID());
                this->addSensorFailTime(myBuckets[permNode].GetFailCapacity());
            }

            // iterate through edges
            for (int k = 0; k < myBuckets[permNode].GetDegree(); ++k) {

                int updateNode1 = myBuckets[permNode].GetReachableNode(k);

                if (!myBuckets[updateNode1].IsReached()){
                    myList.push_back(updateNode1);
                    myBuckets[updateNode1].UpdateReached();
                }


                if (!myBuckets[updateNode1].IsPermanent()) {
                    double updateLabel = min(myBuckets[permNode].GetFailCapacity(), myBuckets[updateNode1].GetFailTime());
                    myBuckets[updateNode1].UpdateFailCapacity(updateLabel);
                }
            }
        }
        numPermanent = numPermanent + 1;
    }

    // If a target is not connected to the sink then it will never be marked permanent, so it will never be added to the
    // order of target failures. This for loop checks for such a case, and adds the targets to ths list with a fail time of zero
    if (orderedTargetFail.size() < numTar){
        SimulationNodeMC *ptr = myBuckets;
        ptr = ptr + (numNodes - numTar);
        for (int i = 0; i < numTar ; ++i) {
            if (!ptr->IsPermanent()){
                this->addTargetFailOrder(ptr->GetNodeID());
                this->addTargetFailTime(ptr->GetFailCapacity());
            }
            ptr = ptr + 1;
        }
    }

    // end time record
    auto finish = chrono::high_resolution_clock::now();
    // calculate execution time
    chrono::duration<double> elapsed = finish - start;
    timeDuration = elapsed.count();

}