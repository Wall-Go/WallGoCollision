#include "kinematics.h"


//*****************Basic routines that work for generic four-vectors***********************


// Scalar product of two four-vector
double SP4(FourVector a, FourVector b){
    double* aC=a.components();
    double* bC=b.components();
return aC[0]*bC[0]-aC[1]*bC[1]-aC[2]*bC[2]-aC[3]*bC[3];
}


// Vector product of the spatial components
double SP3(FourVector a, FourVector b){
    double* aC=a.components();
    double* bC=b.components();
return aC[1]*bC[1]+aC[2]*bC[2]+aC[3]*bC[3];
}



// Adds two four-vectors
FourVector addFV(FourVector a, FourVector b){
    FourVector vec;

    //We extract the components and then add them to create the new four-vector
    double* aC=a.components();
    double* bC=b.components();
    double comp[4];

    //Help vector for sum of the two vectors
    comp[0]=*(aC)+*(bC);
    comp[1]=*(aC+1)+*(bC+1);
    comp[2]=*(aC+2)+*(bC+2);
    comp[3]=*(aC+3)+*(bC+3);

    //Sets the new vector to be the sum of a and b
    vec.setComponents(comp);

return vec;
}

// Subtracts the four-vector b from a. So vec=a-b
FourVector subtractFV(FourVector a, FourVector b){
    FourVector vec;

    //We extract the components and then add them to create the new four-vector
    double* aC=a.components();
    double* bC=b.components();
    double comp[4];

    //Help vector for sum of the two vectors
    comp[0]=*(aC)-*(bC);
    comp[1]=*(aC+1)-*(bC+1);
    comp[2]=*(aC+2)-*(bC+2);
    comp[3]=*(aC+3)-*(bC+3);

    //Sets the new vector to be the sum of a and b
    vec.setComponents(comp);

return vec;
}




//*****************Specialized routines for our collision integrals*********************



// Calculates the norm for the p2 vector
double P2Norm(FourVector k,FourVector p, FourVector p2){
return SP4(p,k)/(SP4(p2,p)+SP4(p2,k));
}



// Creates s,t,y Mandelstam invariants
void CreateInvariants(double mandelstam[3],FourVector p,FourVector k, FourVector p2, FourVector k2){
//help four-vector
FourVector helpFV;

//The s invariant
helpFV=addFV(p,k);
mandelstam[0]=helpFV.norm();

//The t invariant
helpFV=subtractFV(p,p2);
mandelstam[1]=helpFV.norm();

//The u invariant
helpFV=subtractFV(p,k2);
mandelstam[2]=helpFV.norm();
}




