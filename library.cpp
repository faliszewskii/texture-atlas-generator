#include "library.h"

#include <iostream>
#include <eigen3/Eigen/Core>

int main()
{
   Eigen::Vector3f f{1,1,1};
   f.cross3(f);
}