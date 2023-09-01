#include "heap.h"

void Heap::pushUp(int idx) {
    while (idx > 0 && heap[idx].second < heap[(idx-1)/2].second) {
        nodeMap[heap[idx].first] = (idx - 1) / 2;
        nodeMap[heap[(idx-1)/2].first] = idx;
        std::swap(heap[idx], heap[(idx-1)/2]);
        idx = (idx - 1) / 2;
    }
}

void Heap::pushDown(int idx) {
    while (2 * idx + 1 < heap.size()) {
        int k = 2 * idx + 1;
        if (k + 1 < heap.size() && heap[k+1].second < heap[k].second) {
            k = k + 1;
        }
        if (heap[idx].second > heap[k].second) {
            nodeMap[heap[idx].first] = k;
            nodeMap[heap[k].first] = idx;
            std::swap(heap[idx], heap[k]);
            idx = k;
        } else {
            break;
        }
    }
}

void Heap::insert(int key, int value) {
    if (nodeMap.find(key) != nodeMap.end()) {
        int idx = nodeMap[key];
        int tmpValue = heap[idx].second;
        heap[idx].second = value;
        if (tmpValue > value) {
            pushUp(idx);
        } else {
            pushDown(idx);
        }
    } else {
        heap.push_back(std::pair<int, int>(key, value));
        nodeMap[key] = heap.size() - 1;
        pushUp(heap.size() - 1);
    }
}

void Heap::adjust(int key, int value) {
    insert(key, value);
}

int Heap::pop() {
    if (heap.size() == 0) {
        return -1;
    }
    int toReturn = heap[0].first;
    nodeMap.erase(heap[0].first);
    nodeMap[heap[heap.size()-1].first] = 0;
    std::swap(heap[0], heap[heap.size()-1]);
    heap.pop_back();
    pushDown(0);
    return toReturn;
}

bool Heap::in(int key) {
    return nodeMap.find(key) != nodeMap.end();
}

int Heap::size() {
    return nodeMap.size();
}

bool Heap::empty() {
    return size() <= 0;
}