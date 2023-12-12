#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from dynamics_simulation_example2 import Dynamics_sim

"""
Created on Sat Jan 19 18:52:50 2019
@author: Guang Yang
"""
class STL_CBF_Planning:
    '''
    The Problem class contains problem parameters
    '''
    def __init__(self):
        self.model = Model("Double Integrators")
        self.M = 99999 #Big M
        self.N = 10 # Number of horizon steps
        self.t_f = 1 # End of STL time interval
        self.tau_list = [self.t_f/(self.N-1) for i in range(self.N-1)] # List that stores problem parameters
        self.k_cbf= 10 # CBF coefficient
        self.num_state = 2 # two states for double intergrator x1 and x2
        self.num_jordan_block = 2
        self.size_jordan_block = 2 #s(i)
        self.X_0 = (1,-1) #initial state
        self.num_predicates = 2
        self.num_CBF_constraints = 1
        self.delta = 0.2

    def initialize_model(self):
        '''
        OUTPUT: model with state variables, system Dynamics constraints and cost function
        The function takes in Gurobi model and generate state variables (x,u) according to the horizon length

        Gurobi Variable Format:
        {(0, 0): <gurobi.Var C0>,   x1(t_0)
         (0, 1): <gurobi.Var C1>,   x2(t_0)
         (1, 0): <gurobi.Var C2>,   x1(t_1)
         (1, 1): <gurobi.Var C3>,   x2(t_1)
         (2, 0): <gurobi.Var C4>,   x1(t_2)
         (2, 1): <gurobi.Var C5>}   x2(t_2)
         ....
         (N, 0): <gurobi.Var C 2N>,   x1(t_N)
         (N, 1): <gurobi.Var C 2N+1>} x2(t_N)
        '''
        # Initialize State Variables, Total State Variables = 2N+(N-1) = 3N-1
        print("Initialize state variables, initial condition, integer variable and cost function")
        self.U = self.model.addVars(self.N-1,name="u", vtype=GRB.CONTINUOUS,lb=-200,ub=200)
        self.X = self.model.addVars(self.N,self.num_state,name="x",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.z_STL = self.model.addVar(vtype=GRB.BINARY, name="z_STL") #Final STL indicator variable, z_STL = 1 <-> STL is satisfied
        self.Z = self.model.addVars(self.N,self.num_predicates, name="Z_integer",vtype=GRB.BINARY)
        self.mu = self.model.addVars(self.N,self.num_predicates, name="mu",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY) #Predicate
        #self.rho= self.model.addVars(name="STL_Robustness", obj = -1.0)
        self.Z_cbf_x = self.model.addVars(self.N,self.num_CBF_constraints, name="Z_cbf_x",vtype=GRB.BINARY)
        self.Z_cbf_u = self.model.addVars(self.N,self.num_CBF_constraints, name="Z_cbf_u",vtype=GRB.BINARY)

        # Cost Function = \intergral_{0,t_f} u^Tu
        self.model.setObjective(quicksum((self.U[k]*self.U[k])*self.tau_list[k]*0.5  for k in range(self.N-1)),GRB.MINIMIZE)

        # Initilaize Dynamics Constraint
        self.model.addConstr(self.X[0,0] == self.X_0[0],"initial_state_constr_x1")
        self.model.addConstr(self.X[0,1] == self.X_0[1],"initial_state_constr_x2")

        for i in range(self.N-1):
            self.model.addConstr(self.X[i,0] + self.X[i,1]*self.tau_list[i] + (self.tau_list[i]**2/2)*self.U[i] == self.X[i+1,0],"dynamics_constr_x1_step"+str(i))
            self.model.addConstr(self.X[i,1] + self.U[i]*self.tau_list[i] == self.X[i+1,1],"dynamics_constr_x2_step"+str(i))


    def save_model(self):
        self.model.write("test_model.lp")

    def always_predicate(self):
        '''
        G[3,4] x2>3+delta
        G[5,6] x2<-1-delta

        where delta increases the robustness of the signal

        predicate 0: x2-3 >= 0
        predicate 1: -x2 -1 >=0
        '''
        for t_step in range(0,self.N):
            self.model.addConstr(self.mu[t_step,0] == (-0.5*self.X[t_step,0]-2+self.X[t_step,1]),"Always_Predicate_0")

        for t_step in range(0,self.N):
            self.model.addConstr(self.mu[t_step,1] == (-self.X[t_step,1]-8-self.delta),"Always_Predicate_1")
        '''
        Integer to Predicates Encoding,
        Z[t,0] = 1 <==> mu[t,0] = x_1(t) > 0
        '''
        for t_step in range(0,self.N):
            for i in range(self.num_predicates):
                self.model.addConstr(-self.mu[t_step,i] <= self.M * (1-self.Z[t_step,i]))
                self.model.addConstr(self.mu[t_step,i] <= self.M * self.Z[t_step,i])

        '''
        Always Encoding
        z_0 <= Z[t,0], t=0,...,N
        '''
        self.z_0 = self.model.addVar(name="Z_integer",vtype=GRB.BINARY)
        self.z_1 = self.model.addVar(name="Z_integer",vtype=GRB.BINARY)

        for t_step in range(3,5): # enable predicate 0 at t=3 and t=4
            self.model.addConstr(self.z_0 <= self.Z[t_step,0],"Always_Integer_Constraint_predicate_0")

        for t_step in range(7,9): # enable predicate 0 at t=3 and t=4
            self.model.addConstr(self.z_1 <= self.Z[t_step,1],"Always_Integer_Constraint_predicate_0")

        self.model.addConstr(self.z_0 == self.z_STL, "Enforce constraint")
        self.model.addConstr(self.z_1 == self.z_STL, "Enforce constraint")



    def CBF_constraint(self):
        '''
        CBF Constraint:
        h(x) = mu = x-3
        Active between t_3, t_4
        '''
        self.beta_X = self.model.addVars(self.N,name="beta_x_list",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U = self.model.addVars(self.N,name="beta_u_list",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)

        for t_step in range(0,self.N-1):
            self.model.addConstr((self.beta_X[t_step] + self.beta_U[t_step] == 0.0),"beta_constr_for_max")

        # Integer Varibale for removing postive contributing term
        # Z_cbf_x[k,i,j], where k is the time step, i is the eigen value indicator, j is the jordan block indicator
        self.Z_cbf_x = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_x_list",vtype=GRB.BINARY)
        self.Z_cbf_u = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_u_list",vtype=GRB.BINARY)

        '''
         Constraint on x, zeta_x
         if: -k*x(t_k) < 0, then enforce -k*x(t_k) + alpha_x*(k*x2_max) >= 0

         if: -ut < 0, then enforce -u - k*u*tau - alpha_u*(k*x2_min) >= 0

         if: k*x(t_k) < 0, then enforce k*x(t_k) - alpha_x_2*(k*x2_min) >= 0

         if: ut < 0, then enforce u + k*u*tau - alpha_u_2*(k*x2_min) >= 0
        '''
        for t_step in range(3,4):
            # x2(t) >= 3
            self.model.addConstr(((self.k_cbf*self.X[t_step,1] - self.k_cbf*3.0 + self.beta_X[t_step]) >= 0) ,"CBF_x_constr_step"+str(t_step))
            self.model.addConstr(((self.U[t_step]+(self.k_cbf*self.U[t_step]*2*self.tau_list[t_step])+ self.beta_U[t_step])>= 0),"CBF_u_constr"+str(t_step))

        self.model.addConstr(self.U[3] == self.U[4])

    def solve_MILP(self):
        x_traj = []
        u_traj = []
        x1_traj = []
        x2_traj = []

        self.model.addConstr(self.z_STL == 1)
        self.model.optimize()

        status = self.model.status
        if status == GRB.Status.OPTIMAL:
            for i in range(self.N):
                for j in range(self.num_state):
                    x_traj.append(self.X[i,j].x)


            for i in range(self.N-1):
                u_traj.append(self.U[i].x)

            x1_traj = x_traj[0::2]
            x2_traj = x_traj[1::2]

            # Print result
            for variable in self.model.getVars():
                print('The decision variable %s has value of %g' % (variable.varName, variable.x))
        elif status == GRB.Status.INFEASIBLE:
            print('Optimization was stopped with status %d' % status)
            self.model.computeIIS()
            self.model.write("infeasible_model.ilp")
            for c in self.model.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)

        return status,x1_traj, x2_traj, u_traj


    def plot_signals(self,x1_traj,x2_traj,u_traj):
        t = np.linspace(0,self.t_f,self.N)
        plt.plot(x1_traj,x2_traj,'ro')
        plt.title("x_1 vs x_2")
        plt.grid()


if __name__ =="__main__":
    planning_obj = STL_CBF_Planning()
    planning_obj.initialize_model()


    planning_obj.always_predicate()
    #planning_obj.CBF_constraint()

    # Solve and save the MILP model
    planning_obj.save_model()
    status,x1_traj_discrete, x2_traj_discrete, u_traj = planning_obj.solve_MILP()

    if status == GRB.Status.OPTIMAL:
        planning_obj.plot_signals(x1_traj_discrete,x2_traj_discrete, u_traj)

        # Plot continous trajectory
        sim_obj = Dynamics_sim(planning_obj.t_f, planning_obj.N,planning_obj.X_0,u_traj,x1_traj_discrete, x2_traj_discrete, planning_obj.tau_list)
        sim_obj.generate_traj_plots()

        # Plot continous CBF margin
