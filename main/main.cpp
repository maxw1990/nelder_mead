//
//  main.cpp
//  nelder_mead
//
//  Created by Max Werner on 02.02.21.
//
#include <algorithm>
#include <functional>
#include <array>
#include <string_view>
#include "3rdPartner/eigen-3.4.0/Eigen/Dense"
#include "lib/vertix.h"

std::vector<double> addvectors(std::vector<double> a,std::vector<double> b){
    std::vector<double> c;
    for(int i = 0; i< a.size(); i++){
        c.push_back(a[i] + b[i]);
    }
    return c;
}

std::vector<std::vector<double>> createEinheitMatrix(int size){
    int i = size;
    std::vector<std::vector<double>> Einheit;
    for (int j=0; j<i; j++) {//set up the identity matrix using the Kronecker Delta
            std::vector<double> temp2; // temporary vector to push onto identity
            for (int k=0; k<i; k++){
                if(j==k) {
                temp2.push_back(1);
                }
                else {
                temp2.push_back(0);
                }
            }
            Einheit.push_back(temp2);
    };
    return Einheit;
};

class Simplex{
public:
    Simplex(){
    };
    ~Simplex(){
        for(auto &i: eckpunkte){
            i.~Vertix();
        }
    };
    
    void createSimplexStartVertexes(vertexLib::Vertix *start, std::function<double(std::vector<double>)> f){
        std::vector<std::vector<double>> Einheit =   createEinheitMatrix((int)start->vec.size());
        eckpunkte.push_back(*start);
        for(int k = 0;k< start->vec.size();k++){
        eckpunkte.push_back(CreateVertix2(&start->vec, &Einheit[k]));
        }
        this->i = (int)start->vec.size();
    }
    
    std::vector<vertexLib::Vertix> *getVertices(){
        return &eckpunkte;
    }
    
    vertexLib::Vertix *getVertex(int i){
        return &eckpunkte[i];
    }
    
    int getDimension(){
        return i;
    };
    
    void setVertex2(int i, vertexLib::Vertix *vec){
        eckpunkte[i].vec = vec->vec;
    }
    
    void printfvertix(){
        std::cout<<"Dimension: "<< getDimension()<<std::endl;
        for (int j=0; j<eckpunkte.size(); j++) {
            for (int k=0; k<eckpunkte[0].vec.size(); k++) {
                std::cout<<' '<< eckpunkte[j].vec[k];
            }
            std::cout<<std::endl;
        }
    }

protected:
    std::vector<vertexLib::Vertix> eckpunkte;
    int i;
};

vertexLib::Vertix *createmean(Simplex *simp){
    vertexLib::Vertix *a;
    a = new vertexLib::Vertix;
    for(int i = 0;i<simp->getDimension(); ++i){
        a->vec.push_back(0);
    }
    for(int j = 0;j < simp->getDimension(); ++j){
    for(int i = 0; i< simp->getVertex(0)->vec.size(); i++){
        a->vec[i] += simp->getVertex(j)->vec[i];
    }
    }
    for(int i = 0; i< a->vec.size(); i++){
        a->vec[i] = 1/(double)simp->getDimension() * a->vec[i];
    }
    return a;
};

vertexLib::Vertix *reflect(Simplex *simp, vertexLib::Vertix *mean, double alpha){
    vertexLib::Vertix *reflect_m;
    reflect_m = new vertexLib::Vertix;
    for(int i = 0;i<simp->getDimension(); ++i){
        reflect_m->vec.push_back(0);
    }
    for(int i = 0; i< simp->getVertex(0)->vec.size(); ++i){
        reflect_m->vec[i] = mean->vec[i] + alpha*(mean->vec[i] - simp->getVertex(simp->getDimension())->vec[i]);
    }
    return reflect_m;
};

vertexLib::Vertix *expand(Simplex *simp, vertexLib::Vertix *reflect, vertexLib::Vertix *mean, double gamma){
    vertexLib::Vertix *expand_m;
    expand_m = new vertexLib::Vertix;
    for(int i = 0;i<simp->getDimension(); ++i){
        expand_m->vec.push_back(0);
    }
    for(int i = 0; i< simp->getVertex(0)->vec.size(); ++i){
        expand_m->vec[i] = mean->vec[i] + gamma*(reflect->vec[i] - mean->vec[i]);
    }
    return expand_m;
};

vertexLib::Vertix *contracted(Simplex *simp, vertexLib::Vertix *mean, double rho){
    vertexLib::Vertix *contracted_m;
    contracted_m = new vertexLib::Vertix;
    for(int i = 0;i<simp->getDimension(); ++i){
        contracted_m->vec.push_back(0);
    }
    for(int i = 0; i< simp->getVertex(0)->vec.size(); ++i){
        contracted_m->vec[i] = mean->vec[i] + rho*( simp->getVertex(simp->getDimension())->vec[i] - mean->vec[i] );
    }
    return contracted_m;
};

void sort_simplex(Simplex *simplex,std::function<double(std::vector<double>)> func){
    std::sort(simplex->getVertices()->begin(), simplex->getVertices()->end(), [&func](const vertexLib::Vertix &a, const vertexLib::Vertix &b) {
      return func(a.vec) < func(b.vec);
    });
}

double std_deviation(Simplex * simplex, vertexLib::Vertix *mean, std::function<double(std::vector<double>)> func){
    double std_dev = 0.0;
    double variance = 0.0;
    for(int i = 0; i < simplex->getVertices()->size(); ++i){
        variance += std::pow(func(simplex->getVertex(i)->vec) - func(mean->vec),2);
    };
    variance = (1/(double)simplex->getDimension()) * variance;
    std_dev = std::sqrt(variance);
    return std_dev;
};

void shrink(Simplex *simp, double theta){
    vertexLib::Vertix a;
    std::vector<double> vector;
    vector.resize(simp->getDimension());
    for(int i = 0; i< simp->getDimension(); ++i){
        for(int j = 0; j < simp->getVertex(0)->vec.size(); ++j){
            vector[i] = simp->getVertex(0)->vec[j] + theta*( simp->getVertex(i)->vec[j] - simp->getVertex(0)->vec[j]);
        }
        a.vec = vector;
        simp->setVertex2(i, &a);
    }
};


double func(std::vector<double> vec){
return(vec[0]*vec[0]+vec[1]-11)*(vec[0]*vec[0]+vec[1]-11)+(vec[0]+vec[1]*vec[1]-7)*(vec[0]+vec[1]*vec[1]-7);
}

double funcSim(std::vector<double> vec){
return(vec[0]*vec[0])+(vec[1]*vec[1]);
}

double booth(std::vector<double> vec){
return(vec[0]+2.0*vec[1]-7)*(vec[0]+2.0*vec[1]-7)+(vec[1]+2.0*vec[0]-5)*(vec[1]+2.0*vec[0]-5);
}

double func4(std::vector<double> vec){
    return (vec[0]+10*vec[1])*(vec[0]+10*vec[1])+5*(vec[2]-vec[3])*(vec[2]-vec[3])+std::pow((vec[1]-2*vec[2]),4)+10*std::pow((vec[0]-vec[3]),4);
}

double func2(std::vector<double> vec){
    return(vec[0]*vec[0]+vec[1]-11)*(vec[0]*vec[0]+vec[1]-11)+(vec[0]+vec[1]*vec[1]-7)*(vec[0]+vec[1]*vec[1]-7);
}

class Minimization{
public:
    ~Minimization(){simp->~Simplex();};
    virtual void print() = 0;
    virtual void Start_minimizing() = 0;
protected:
    Simplex * simp;
};

class NelderMead : public Minimization{
public:
    NelderMead(std::function<double(std::vector<double>)> minimize,std::vector<double> start, double alpha, double beta, double gamma, double std_dev_min):
     alpha(alpha), beta(beta), gamma(gamma),std_dev_min(std_dev_min)
    {
        m_minimize = minimize;
        vertexLib::Vertix *vert = new vertexLib::Vertix(start);
        Simplex *simp1 = new Simplex();
        simp1->createSimplexStartVertexes(vert, minimize);
        this->simp = simp1;
    };
    
    void Start_minimizing() override{
        vertexLib::Vertix *meanstart = createmean(simp);
        double std_dev = std_deviation(simp, meanstart,m_minimize);
        std::cout << std_deviation(simp, meanstart,m_minimize) << std::endl;
        int count = 0;
        while(std_dev > std_dev_min){
        //first step
            std::cout << count << std::endl;
            count++;
        sort_simplex(simp, m_minimize);
        //second step
        vertexLib::Vertix *mean = createmean(simp);
        std_dev = std_deviation(simp, mean,m_minimize);
        std::cout << std_deviation(simp, mean,m_minimize) << std::endl;
        //third step reflection
        vertexLib::Vertix *reflect_x = reflect(simp,mean, alpha);
        if(m_minimize(simp->getVertex(0)->vec) <= m_minimize(reflect_x->vec) && m_minimize(reflect_x->vec) < m_minimize(simp->getVertex(1)->vec)){
            simp->setVertex2(simp->getDimension(), reflect_x);
            delete reflect_x;
            continue;
        }
        //firth step expansion
        if(m_minimize(reflect_x->vec) < m_minimize(simp->getVertex(0)->vec)){
            vertexLib::Vertix *expand_m = expand(simp, reflect_x, mean, beta);
            if(m_minimize(expand_m->vec) < m_minimize(reflect_x->vec)){
                //simp1.setVertex(simp1.getDimension(), expand_m->vec);
                simp->setVertex2(simp->getDimension(), expand_m);
                delete expand_m;
                delete reflect_x;
                continue;
            }
            else{
                simp->setVertex2(simp->getDimension(), reflect_x);
                delete expand_m;
                delete reflect_x;
                continue;
            }
        }
        //fith step
        vertexLib::Vertix *contracted_m = contracted(simp, mean, gamma);
        if(m_minimize(contracted_m->vec) < m_minimize(simp->getVertex(simp->getDimension())->vec)){
            simp->setVertex2(simp->getDimension(), contracted_m);
            delete contracted_m;
            continue;
        }
        
        //sixth step
            shrink(simp, gamma);
        }
        this->print();
    }
    
    void print() override {
        std::cout<<"Vertix: "<<std::endl;
        simp->printfvertix();
         std::cout<<"Value: "<<std::endl;
        std::cout << m_minimize(simp->getVertex(0)->vec) << std::endl;
        std::cout << m_minimize(simp->getVertex(1)->vec) << std::endl;
        std::cout << m_minimize(simp->getVertex(2)->vec) << std::endl;
    };
protected:
    double alpha, beta, gamma, std_dev_min;
    std::function<double(std::vector<double>)> m_minimize;
};

int main(int argc, const char * argv[]) {

    //f_display = func2;
    std::vector<double> start = {-4,-5};
    //Minimization *minimizing = new NelderMead(start,1,2,0.5,1e-08);
    NelderMead *minimizing = new NelderMead(booth,start,1,2,0.5,1e-08);
    minimizing->Start_minimizing();
    //Minimization* pc = dynamic_cast<Minimization*>(minimizing);
    //pc->Start_minimizing();
    //minimizing->print();
    //minimizing->Start_minimizing();
    return 0;
    
}


