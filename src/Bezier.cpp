#include "Bezier.h"

using namespace arma;

Bezier::Bezier(int degree)
{
    _degree = degree;
    Set_Factores();
    Set_t(100);
}

Bezier::Bezier(mat P)
{
    _degree = P.n_rows - 1;
    Set_Factores();
    Set_ControlPoints(P);
    Set_t(100);
}

void Bezier::Set_Factores()
{
    switch (_degree)
    {
        case 3:
            _factores = {{-1, 3, -3, 1},
                        {3, -6, 3, 0},
                        {-3, 3, 0, 0},
                        {1, 0, 0, 0}};
            break;
        case 5:
            _factores = {{-1, 5, -10, 10, -5, 1},
                        {5, -20, 30, -20, 5, 0},
                        {-10, 30, -30, 10, 0, 0},
                        {10, -20, 10, 0, 0, 0},
                        {-5, 5, 0, 0, 0, 0},
                        {1, 0, 0, 0, 0, 0}};
            break;
        default: throw std::string("Bezier degree must be 3 or 5");
    }
}

void Bezier::Set_ControlPoints(mat P)
{
    _pcontrol = P;
    _Cx = _factores*_pcontrol.col(0);
    _Cy = _factores*_pcontrol.col(1);

    _flag_control_points = true;
    _flag_curva = false;
    _flag_curvatura = false;
    _flag_tangente = false;
    _flag_normal = false;
    _flag_equi = false;
    _flag_curvaLength = false;
}

mat Bezier::Get_ControlPoints()
{
    return _pcontrol;
}

mat Bezier::Get_ControlPoints(uvec pc)
{
    return _pcontrol.rows(pc);
}

void Bezier::Set_t(vec t)
{
    _t = t;
    _t_length = t.n_elem;

    Set_Tmatrix();
}

void Bezier::Set_t(int t)
{
    _t = linspace(0, 1, t);
    _t_length = t;

    Set_Tmatrix();
}

void Bezier::Set_t(double t)
{
    _t = t;
    _t_length = 1;

    Set_Tmatrix();
}

void Bezier::Set_Tmatrix()
{
    _T.set_size(_t_length, _degree+1);
    _T_dot.set_size(_t_length, _degree+1);
    _T_dotdot.set_size(_t_length, _degree+1);

    vec t2, t3;
    t3 = pow(_t, 3);
    t2 = pow(_t, 2);

    _T.col(_degree-3) = t3;
    _T.col(_degree-2) = t2;
    _T.col(_degree-1) = _t;
    _T.col(_degree).ones();

    _T_dot.col(_degree-3) = 3 * t2;
    _T_dot.col(_degree-2) = 2 * _t;
    _T_dot.col(_degree-1).ones();
    _T_dot.col(_degree).zeros();

    _T_dotdot.col(_degree-3) =  6 * _t;
    _T_dotdot.col(_degree-2).fill(2);
    _T_dotdot.cols(_degree-1,_degree).zeros();

    if (_degree == 5)
    {
        _T.col(0) = pow(_t, 5);
        _T.col(1) = pow(_t, 4);

        _T_dot.col(0) = 5 * _T.col(1);
        _T_dot.col(1) = 4 * t3;

        _T_dotdot.col(0) = 20 * t3;
        _T_dotdot.col(1) = 12 * t2;
    }
    _flag_curva = false;
    _flag_curvatura = false;
    _flag_tangente = false;
    _flag_normal = false;
    _flag_equi = false;
    _flag_curvaLength = false;
}

mat Bezier::Curva()
{
    if (!_flag_curva)
    {
        vec Bx(_t_length), By(_t_length);

        Bx = _T*_Cx;
        By = _T*_Cy;

        _curva = join_rows(Bx, By);
        _flag_curva = true;
    }
    return _curva;
}

mat Bezier::Curva(vec t)
{
    Set_t(t);
    return Curva();
}

mat Bezier::Curva(int t)
{
    if (t != _t_length) {Set_t(t);}
    return Curva();
}

mat Bezier::Curva(double t)
{
    Set_t(t);
    return Curva();
}

mat Bezier::Tangente()
{
    if (!_flag_tangente)
    {
        vec Bx, By;

        Bx = _T_dot*_Cx;
        By = _T_dot*_Cy;

        _tangente = join_rows(Bx, By);
        _flag_tangente = true;
    }
    return _tangente;
}

mat Bezier::Tangente(vec t)
{
    Set_t(t);
    return Tangente();
}

mat Bezier::Tangente(int t)
{
    if (t != _t_length) {Set_t(t);}
    return Tangente();
}

mat Bezier::Tangente(double t)
{
    Set_t(t);
    return Tangente();
}

mat Bezier::Normal()
{
    if (!_flag_tangente)
    {
        vec Bx, By;

        Bx = _T_dot*_Cx;
        By = _T_dot*_Cy;

        _tangente = join_rows(Bx, By);
        _flag_tangente = true;
        _normal = join_rows(-By, Bx);
        _flag_normal = true;
    }
    return _normal;
}

mat Bezier::Normal(vec t)
{
    Set_t(t);
    return Normal();
}

mat Bezier::Normal(int t)
{
    if (t != _t_length) {Set_t(t);}
    return Normal();
}

mat Bezier::Normal(double t)
{
    Set_t(t);
    return Normal();
}

mat Bezier::Curvatura()
{
    if (!_flag_curvatura)
    {
        vec dx, dy, ddx, ddy;

        dx = _T_dot*_Cx;
        dy = _T_dot*_Cy;

        ddx = _T_dotdot*_Cx;
        ddy = _T_dotdot*_Cy;
        _curvatura = (dx%ddy - dy%ddx)/pow(square(dx) + square(dy),1.5);
        _flag_curvatura = true;
    }
    return _curvatura;
}

mat Bezier::Curvatura(vec t)
{
    Set_t(t);
    return Curvatura();
}

mat Bezier::Curvatura(int t)
{
    if (t != _t_length) {Set_t(t);}
    return Curvatura();
}

mat Bezier::Curvatura(double t)
{
    Set_t(t);
    return Curvatura();
}

mat Bezier::CurvaEquid(double dbp)
{
    if (!_flag_equi)
    {
        if (!_flag_curvaLength) {CurvaLength();}

        // divide 't' para tener puntos "~equidistantes" (en t)
        unsigned int np = round(_cum_d(_cum_d.n_elem-1)/dbp);
        if (np > 10000)
        {
            _curva.clear();
            _flag_curva = false;
            return _curva;
        }

        vec dp = linspace(0,_cum_d(_cum_d.n_elem-1),np);
        interp1(_cum_d, _t, dp, _t_equi);
        //Set_t(_t_equi);
        _curva = Curva(_t_equi);
        _tangente = Tangente(_t_equi);
        _curvatura = Curvatura(_t_equi);

        _flag_equi = true;
    }
    return _curva;
}

