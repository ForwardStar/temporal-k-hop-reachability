#ifndef COMMONFUNCTIONS
#define COMMONFUNCTIONS

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <map>
#include <iomanip>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <sstream>
#include <stack>
#include <string>

void putCCNumAndSize(int n, int *p);

void putProcess(double procedure, unsigned long long time_used);

std::stringstream timeFormatting(unsigned long long seconds);

unsigned long long currentTime();

#endif