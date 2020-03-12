#include "include/Bezier.h"
#include <armadillo>

int main()
{

// Create object for a cubic Bezier curve
Bezier cubic_curve(3);

arma::mat control_points;

control_points = {  {0, 1}, 
                    {2, 3}, 
                    {4, 6}, 
                    {7, 4}};

cubic_curve.Set_ControlPoints(control_points);
std::cout << "Cubic curve defined by the following control points:': " << std::endl << cubic_curve.Get_ControlPoints() << std::endl;

// Compute curve 2D point coordinates at 't = 0.5'
std::cout << "Coordinates at 't = 0.5': " << cubic_curve.Curva(0.5) << std::endl;


// Create object for a quintic Bezier curve
Bezier quintic_curve(5);
}