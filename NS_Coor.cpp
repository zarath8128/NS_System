#include "NS.h"

using namespace NS;

NS_Coor::NS_Coor(double delta, double origin)
	:delta(delta), origin(origin)
{}

NS_Coor::operator double(){return origin + delta*i;}
double NS_Coor::operator()(int i){return origin + delta*i;}
