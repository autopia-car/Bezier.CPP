#ifndef BEZIER_H
#define BEZIER_H

#include <armadillo>
#include <stdexcept>

using namespace std;
using namespace arma;

class Bezier
{
    public:
        Bezier(int degree);
        Bezier(mat P);

        void Set_t(vec t);
        void Set_t(int t);
        void Set_t(double t);
        void Set_ControlPoints(mat P);

        mat Get_ControlPoints();
        mat Get_ControlPoints(uvec pc);
        mat Curva();
        mat Curva(vec t);
        mat Curva(int t);
        mat Curva(double t);
        mat Tangente();
        mat Tangente(vec t);
        mat Tangente(int t);
        mat Tangente(double t);
        mat Normal();
        mat Normal(vec t);
        mat Normal(int t);
        mat Normal(double t);
        mat Curvatura();
        mat Curvatura(vec t);
        mat Curvatura(int t);
        mat Curvatura(double t);

    protected:

    private:
        int _degree;
        vec _Cx, _Cy; // For matrix calculations
        vec _cum_d; // Cumulative distances among discretized points
        int _t_length;

        mat _factores;
        mat _t;
        mat _T;
        mat _T_dot;
        mat _T_dotdot;
        mat _t_equi;

        mat _pcontrol;
        mat _curva;
        mat _tangente;
        mat _normal;
        mat _curvatura;
        double _curve_length;

        bool _flag_control_points = false;
        bool _flag_curva = false;
        bool _flag_curvatura = false;
        bool _flag_tangente = false;
        bool _flag_normal = false;
        bool _flag_equi = false;
        bool _flag_curvaLength = false;

        void Set_Tmatrix();
        void Set_Factores();
};

#endif // BEZIER_H
