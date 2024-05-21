#include <iostream>
#include "../eigen/Eigen/Dense"

int main()
{
    Eigen::Matrix3d A {{2.0, -1.0, 3.0},
                       {1.0, -2.0, -5.0},
                       {6.0, 3.0, -10.0}};
    Eigen::Vector3d b {5.0, 3.0, 2.0};

    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);

    std::cout << "A*x is: \n" << A * x << std::endl;
    std::cout << "b is: \n" << b << std::endl;

    return 0;
}
