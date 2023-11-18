#pragma once
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>

namespace vertexLib{

class Vertix{
public:
    Vertix();
    ~Vertix(){};
    Vertix(std::vector<double> vect);
    std::vector<double> vec;
    void printVertix();
};
}

vertexLib::Vertix CreateVertix2(std::vector<double> *a,std::vector<double> *b);


