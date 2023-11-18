#include "vertix.h"

using namespace vertexLib;

Vertix::Vertix(){};

Vertix::Vertix(std::vector<double> vect):vec(vect){
};

void Vertix::printVertix(){
    for (int j=0; j<vec.size(); j++) {
            std::cout<<' '<< vec[j];
        }
        std::cout<<std::endl;
};

Vertix CreateVertix2(std::vector<double> *a,std::vector<double> *b){
    Vertix avert;
    std::vector<double> c;
    for(int i = 0; i< a->size(); i++){
        c.push_back(a->at(i) + b->at(i));
    }
    avert.vec = c;
    return avert;
};