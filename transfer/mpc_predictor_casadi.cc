#include "mpc_predictor_casadi.h"

#include <casadi/casadi.hpp>
#include <fstream>
#include <math.h>
#include <string>

using namespace casadi;

// Declare model variables
//MX x;
//MX xs;
//MX us;

// Optimization Problem
//Opti* opti;

int counter = 0;

// dx/dt = f(x,u)
MX f(const MX& x, const MX& u) {
    MX xdot = vertcat(u(0), u(1));
    xdot = vertcat(xdot, u(2));
    xdot = vertcat(xdot, u(3));
    xdot = vertcat(xdot, u(4));
    xdot = vertcat(xdot, u(5));
    return vertcat(xdot, u(6));
}

MX endeff_pos(const MX& xs) {
    MX ys = - (2.0*cos(xs(0))*sin(xs(1)))/5.0 - 39.0*sin(xs(3))*(sin(xs(0))*sin(xs(2)) 
            - cos(xs(0))*cos(xs(1))*cos(xs(2)))/100.0 - (39.0*cos(xs(0))*cos(xs(3))*sin(xs(1)))/100.0;
    ys = vertcat(ys, 
            (39.0*sin(xs(3))*(cos(xs(0))*sin(xs(2)) + cos(xs(1))*cos(xs(2))*sin(xs(0))))/100.0 
            - (2.0*sin(xs(0))*sin(xs(1)))/5.0 - (39.0*cos(xs(3))*sin(xs(0))*sin(xs(1)))/100.0);
    return vertcat(ys,
            (2.0*cos(xs(1)))/5.0 + (39.0*cos(xs(1))*cos(xs(3)))/100.0 
            + (39.0*cos(xs(2))*sin(xs(1))*sin(xs(3)))/100.0 + 31.0/100.0);
}

MpcPredictor::MpcPredictor(double dt_) { // dt: sampling time, T: horizon length
    dt = dt_;
    T = N * dt;
}

void MpcPredictor::predict(std::vector<float>& x0_, std::vector<float>& y_ref, std::vector<float>& return_vector) {
    
    ++counter;

    auto opti = Opti();
    Slice all;

    std::vector<double> yt = {y_ref[0], y_ref[1], y_ref[2]};

    //if (counter <= 20)
    //    yt = {0, 0, 1.1};
    //else 
    //    yt = {0, -1, -1};

    auto x = opti.variable(7, N+1);
    auto xs = opti.variable(7, 1);
    MX J;
    auto u = opti.variable(7,N); // control trajectory 
    auto ys = opti.variable(3, 1);
    //auto us = opti.variable(7, 1); // trivially zero

    // Initial condition
    MX x0 = x0_[0];
    x0 = vertcat(x0,x0_[1]);
    x0 = vertcat(x0,x0_[2]);
    x0 = vertcat(x0,x0_[3]);
    x0 = vertcat(x0,x0_[4]);
    x0 = vertcat(x0,x0_[5]);
    x0 = vertcat(x0,x0_[6]);
    opti.subject_to(x(all,0)==x0);
    x(all,0) = x0;

    J = 0;

    // dynamic constraints
    for (int k=0;k<N;++k) {
        auto k1 = f(x(all,k),         u(all,k));
        auto k2 = f(x(all,k)+dt/2*k1, u(all,k));
        auto k3 = f(x(all,k)+dt/2*k2, u(all,k));
        auto k4 = f(x(all,k)+dt*k3,   u(all,k));
        auto x_next = x(all,k) + dt/6*(k1+2*k2+2*k3+k4);
        //opti.subject_to(x(all,k+1)==x_next); // close the gaps
        x(all,k+1) = x_next;
        //J += pow((u(0,k)-0.2),2);
        J += dot((x(all,k)-xs).T(), (x(all,k)-xs).T()) + dot((u(all, k)).T(), (u(all, k)).T());
    }

    // Artificial state constraints
    opti.subject_to(xs<=M_PI);
    opti.subject_to(xs>=-M_PI);

    for (int i = 1; i < N; ++i) {
        // Path constraints
        opti.subject_to(x(0,i)<=2.96);
        opti.subject_to(x(0,i)>=-2.96);
        opti.subject_to(x(1,i)<=-0.1 + 1.57);
        opti.subject_to(x(1,i)>=-3.1 + 1.57);
        opti.subject_to(x(2,i)<=4.0 - 1.05);
        opti.subject_to(x(2,i)>=-1.9 - 1.05);
        opti.subject_to(x(3,i)<=2.09);
        opti.subject_to(x(3,i)>=-2.09);
        opti.subject_to(x(4,i)<=1.35 + 1.57);
        opti.subject_to(x(4,i)>=-3.1 + 1.57);
        opti.subject_to(x(5,i)<=2.09);
        opti.subject_to(x(5,i)>=-2.09);
        opti.subject_to(x(6,i)<=2.96);
        opti.subject_to(x(6,i)>=-2.96);
    }
    for (int i = 0; i < N; ++i) {
        // Input constraints
        opti.subject_to(u(all,i)<=0.5);
        opti.subject_to(u(all,i)>=-0.5);
    }

    // Terminal cost
    //J += dot(x(all,N)-xs, x(all,N)-xs); \\ Zero anyway for equality constraint

    // Terminal constraint
    opti.subject_to(x(all,N)-xs==0);

    // Artificial endeff_pos
    opti.subject_to(ys==endeff_pos(xs));
    // Constraint on end effector
    for (int i = 1; i <= N; ++i)
        opti.subject_to(endeff_pos(x(all,i))(0)<=0.1);

    // Add deviation from artificial endeff_pos to objective
    J += 20*dot((ys-yt).T(), (ys-yt).T()); 

    opti.minimize(J);

    //DMDict arg1;
    //DMDict arg2 = {{"accept_after_max_steps",10}};
    opti.solver("ipopt", {{"iteration_callback_ignore_errors",true}}, {{"max_iter",30}}); // numerical backend
    
    
    auto sol = opti.solve(); //solve 

    std::cout << "\n" << sol.get_str() << "\n";

    std::vector<double> u_star = std::vector<double>(sol.value(u(all,0)));
    for (int i = 0; i < 7; ++i)
        return_vector[i] = u_star[i];
    //std::cout << std::vector<double>(sol.value(u(all,0))) << std::endl;
    //std::cout << std::vector<double>(sol.value(xs)) << std::endl;



}
