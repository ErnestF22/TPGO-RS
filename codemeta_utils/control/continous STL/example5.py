#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from dynamics_simulation_example5 import Dynamics_sim
import time

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
        self.N = 20 # Number of horizon steps
        self.t_f = 10 # End of STL time interval
        self.tau_list = [self.t_f/(self.N-1) for i in range(self.N-1)] # List that stores problem parameters
        self.k_cbf= [10,10 ]# CBF coefficient
        self.x2_max = 10 # CBF constraint
        self.x2_min = -10
        self.tau_A= 1.0
        self.u_min = -1
        self.u_max = 1
        self.num_state = 2 # two states for double intergrator x1 and x2
        self.num_jordan_block = 2
        self.size_jordan_block = 2 #s(i)
        self.X_0 = (0.8,0) #initial state
        self.num_predicates = 2
        self.num_CBF_constraints = 2
        self.delta = 0.3

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
        self.U = self.model.addVars(self.N-1,name="u", vtype=GRB.CONTINUOUS,lb=self.u_min,ub=self.u_max)
        self.X = self.model.addVars(self.N,self.num_state,name="x",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.z_STL = self.model.addVar(vtype=GRB.BINARY, name="z_STL") #Final STL indicator variable, z_STL = 1 <-> STL is satisfied
        self.Z = self.model.addVars(self.N,self.num_predicates, name="Z_integer",vtype=GRB.BINARY)
        self.mu = self.model.addVars(self.N,self.num_predicates, name="mu",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY) #Predicate
        #self.rho= self.model.addVars(name="STL_Robustness", obj = -1.0)
        #self.Z_cbf_x = self.model.addVars(self.N,self.num_CBF_constraints, name="Z_cbf_x",vtype=GRB.BINARY)
        #self.Z_cbf_u = self.model.addVars(self.N,self.num_CBF_constraints, name="Z_cbf_u",vtype=GRB.BINARY)

        # Cost Function
        self.model.setObjective(quicksum((self.U[k]*self.U[k])*self.tau_list[k]*0.5  for k in range(self.N-1)),GRB.MINIMIZE)

        # Initilaize Dynamics Constraint
        self.model.addConstr(self.X[0,0] == self.X_0[0],"initial_state_constr_x1")
        self.model.addConstr(self.X[0,1] == self.X_0[1],"initial_state_constr_x2")

        for i in range(self.N-1):
            self.model.addConstr(self.X[i,0] + self.X[i,1]*self.tau_list[i] + (self.tau_list[i]**2/2)*self.U[i] == self.X[i+1,0],"dynamics_constr_x1_step"+str(i))
            self.model.addConstr(self.X[i,1] + self.U[i]*self.tau_list[i] == self.X[i+1,1],"dynamics_constr_x2_step"+str(i))


    def save_model(self):
        self.model.write("example1_model.lp")


    def add_STL_constr(self):
        '''Input: STL predicate, Time interval in steps that the spec must be valid
           Output: MILP Constraint
        '''
        for t_step in range(0,self.N):
            for i in range(self.num_predicates):
                self.model.addConstr(-self.mu[t_step,i] <= self.M * (1-self.Z[t_step,i]))
                self.model.addConstr(self.mu[t_step,i] <= self.M * self.Z[t_step,i])

        for t_step in range(0,self.N):
            #x1>3
            self.model.addConstr(self.mu[t_step,0] == (self.X[t_step,0] - 3.0-self.delta),"Predicate_1")
            #x1<-2
            self.model.addConstr(self.mu[t_step,1] == (-self.X[t_step,0] - 2.5-self.delta),"Predicate_2")

        for predicate_indicator in range(0,self.num_predicates):
            self.model.addConstr(self.z_STL <= quicksum(self.Z[t_step,predicate_indicator] for t_step in range(0,self.N)),"Eventually_Integer_Constraint")

        for t_step in range(0,self.N):
            self.model.addConstr(self.z_STL >= self.Z[t_step,predicate_indicator],"Eventually_Integer_Constraint")


        # Always [0,2] mu_0 and Always [0,2] mu_1



    def always(self):
        # Relaxtion term
        # alpha_x and alpha_u forms a convex combination for each step, i.e., alpha_x +  alpha_u = 1.0
        #self.alpha_X = self.model.addVars(self.N,name="alpha_x_list",vtype=GRB.CONTINUOUS,lb=0,ub=1)
        #self.alpha_U = self.model.addVars(self.N,name="alpha_u_list",vtype=GRB.CONTINUOUS,lb=0,ub=1)
        self.beta_X = self.model.addVars(self.N,name="beta_x_list",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U = self.model.addVars(self.N,name="beta_u_list",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        #self.beta_X_2 = self.model.addVars(self.N,name="beta_x_list_2",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        #self.beta_U_2 = self.model.addVars(self.N,name="beta_u_list_2",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)

        for t_step in range(0,self.N):
            self.model.addConstr((self.beta_X[t_step] + self.beta_U[t_step] == 0.0),"beta_x_u_Constr")
            #self.model.addConstr((self.beta_X_2[t_step] + self.beta_U_2[t_step] == 0.0),"beta_constr_for_min")



        #self.z_h1 = self.model.addVars(self.N,name="Z_integer_h1",vtype=GRB.BINARY)
        #self.z_h2 = self.model.addVars(self.N,name="Z_integer_h2",vtype=GRB.BINARY)

        # CBF Constraint for h1(x)=x1-3>0 and h2(x)=-x1-2
        for t_step in range(0,self.N-1):
            # x1(t)-3.0 >= 0
            self.model.addConstr((self.Z[t_step,0] == 1) >> (self.k_cbf[0]*self.X[t_step,1]+self.k_cbf[1]*self.X[t_step,1]*self.tau_A + self.k_cbf[1]*self.X[t_step,0]-3.0*self.k_cbf[1]+self.beta_X[t_step] >= 0),"CBF_x_constr_step")
            self.model.addConstr((self.Z[t_step,0] == 1) >> (self.U[t_step]+self.k_cbf[0]*self.U[t_step]*self.tau_A+0.5*self.k_cbf[0]*self.U[t_step]*self.tau_A*self.tau_A+self.beta_U[t_step] >= 0),"CBF_u_constr")


    def add_CBF_constraint(self):
        '''
        CBF Constraint:
        x2 <= x2_max for all time
        x2 >= x2_min for all time
        '''
        # Relaxtion term
        # alpha_x and alpha_u forms a convex combination for each step, i.e., alpha_x +  alpha_u = 1.0
        #self.alpha_X = self.model.addVars(self.N,name="alpha_x_list",vtype=GRB.CONTINUOUS,lb=0,ub=1)
        #self.alpha_U = self.model.addVars(self.N,name="alpha_u_list",vtype=GRB.CONTINUOUS,lb=0,ub=1)
        self.beta_X = self.model.addVars(self.N,name="beta_x_list",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U = self.model.addVars(self.N,name="beta_u_list",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_X_2 = self.model.addVars(self.N,name="beta_x_list_2",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U_2 = self.model.addVars(self.N,name="beta_u_list_2",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)

        #for t_step in range(1,self.N):
        #    self.model.addConstr((self.alpha_X[t_step] + self.alpha_U[t_step] == 1.0),"alpha_constr")

        for t_step in range(0,self.N):
            self.model.addConstr((self.beta_X[t_step] + self.beta_U[t_step] == 0.0),"beta_constr_for_max")
            self.model.addConstr((self.beta_X_2[t_step] + self.beta_U_2[t_step] == 0.0),"beta_constr_for_min")

        #self.alpha_X_2 = self.model.addVars(self.N,name="alpha_x_list",vtype=GRB.CONTINUOUS,lb=0,ub=1)
        #self.alpha_U_2 = self.model.addVars(self.N,name="alpha_u_list",vtype=GRB.CONTINUOUS,lb=0,ub=1)
        #for t_step in range(1,self.N):
        #    self.model.addConstr((self.alpha_X_2[t_step] + self.alpha_U_2[t_step] == 1.0),"convex_combination_constr")


        # Integer Varibale for removing postive contributing term
        # Z_cbf_x[k,i,j], where k is the time step, i is the eigen value indicator, j is the jordan block indicator
        self.Z_cbf_x = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_x_list",vtype=GRB.BINARY)
        self.Z_cbf_u = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_u_list",vtype=GRB.BINARY)
        self.Z_cbf_x_2 = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_x2_list",vtype=GRB.BINARY)
        self.Z_cbf_u_2 = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_u2_list",vtype=GRB.BINARY)

        for t_step in range(0,self.N-1):
        #    for i in range(self.num_CBF_constraints):
                #Integer Encoding for x
                self.model.addConstr(-(self.k_cbf*self.X[t_step,1]) <= self.M * (1-self.Z_cbf_x[t_step,0,0]))
                self.model.addConstr(self.k_cbf*self.X[t_step,1] <= self.M * self.Z_cbf_x[t_step,0,0])

                self.model.addConstr(-(self.k_cbf*self.X[t_step,1]) <= self.M * (1-self.Z_cbf_x_2[t_step,0,0]))
                self.model.addConstr(self.k_cbf*self.X[t_step,1] <= self.M * self.Z_cbf_x_2[t_step,0,0])

                #Integer Encoding for u
                self.model.addConstr(-self.k_cbf*self.U[t_step] <= self.M * (1-self.Z_cbf_u[t_step,0,1]))
                self.model.addConstr(self.k_cbf*self.U[t_step]<= self.M * self.Z_cbf_u[t_step,0,1])

                self.model.addConstr(self.k_cbf*self.U[t_step] <= self.M * (1-self.Z_cbf_u_2[t_step,0,1]))
                self.model.addConstr(-self.k_cbf*self.U[t_step] <= self.M * self.Z_cbf_u_2[t_step,0,1])

        '''
         Constraint on x, zeta_x
         if: -k*x(t_k) < 0, then enforce -k*x(t_k) + alpha_x*(k*x2_max) >= 0

         if: -ut < 0, then enforce -u - k*u*tau - alpha_u*(k*x2_min) >= 0

         if: k*x(t_k) < 0, then enforce k*x(t_k) - alpha_x_2*(k*x2_min) >= 0

         if: ut < 0, then enforce u + k*u*tau - alpha_u_2*(k*x2_min) >= 0
        '''
        for t_step in range(0,self.N-1):
            # X2_max Constraints
            #self.model.addConstr((self.Z_cbf_x[t_step,0,0]==1) >> ((-self.k_cbf*self.X[t_step,1] + self.alpha_X[t_step]*self.k_cbf*self.x2_max+self.beta_X[t_step]) >= 0) ,"CBF_x_constr_step"+str(t_step))
            #self.model.addConstr((self.Z_cbf_u[t_step,0,1]==1) >> ((-self.U[t_step]+(-self.k_cbf*self.U[t_step]*self.tau_list[t_step])+self.alpha_U[t_step]*self.k_cbf*self.x2_max+self.beta_U[t_step]) >= 0),"CBF_u_constr"+str(t_step))
            #self.model.addConstr((-self.U[t_step]-self.k_cbf*self.U[t_step]*self.tau_list[t_step]+self.k_cbf*(-self.X[t_step,1]+self.x2_max) >= 0),"CBF_u_constr"+str(t_step))

            # X2_min Constraints
            self.model.addConstr((self.Z_cbf_u_2[t_step,0,1]==1) >> ((self.k_cbf*self.X[t_step,1] - self.k_cbf*self.x2_min + self.beta_X_2[t_step]) >= 0) ,"CBF_x_constr_step"+str(t_step))
            self.model.addConstr((self.Z_cbf_u_2[t_step,0,1]==1) >> ((self.U[t_step]+(self.k_cbf*self.U[t_step]*self.tau_list[t_step])+ self.beta_U_2[t_step])>= 0),"CBF_u_constr"+str(t_step))

            #self.model.addConstr((self.Z_cbf_x[t_step,0,0]==1) >> ((-self.k_cbf*self.X[t_step,1] + self.k_cbf*self.x2_max +self.beta_X[t_step]) >= 0) ,"CBF_x_constr_step"+str(t_step))
            self.model.addConstr((self.Z_cbf_u[t_step,0,1]==1) >> ((-self.U[t_step]-(self.k_cbf*self.U[t_step]*self.tau_list[t_step])+self.beta_U[t_step]) >= 0),"CBF_u_constr"+str(t_step))
            self.model.addConstr((self.Z_cbf_u[t_step,0,1]==1) >> ((-self.k_cbf*self.X[t_step,1] + self.k_cbf*self.x2_max +self.beta_X[t_step]) >= 0) ,"CBF_x_constr_step"+str(t_step))
            #self.model.addConstr(((-self.U[t_step]-(self.k_cbf*self.U[t_step]*self.tau_list[t_step])+self.beta_U[t_step]) >= 0),"CBF_u_constr"+str(t_step))

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

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharey=False)
        ax1.plot(t, x1_traj,'ro')
        ax1.set_title('$x_1$ vs t')
        ax1.set_xlim([0, self.t_f])
        ax1.set_ylim([-10, 10])

        ax2.plot(t,x2_traj,'ro')
        ax2.set_title('$x_2$ vs t')
        ax2.set_xlim([0, self.t_f])
        ax2.set_ylim([-100, 100])


        #ax3.step(t[0:-1],u_traj)
        ax3.step(t[:-1],u_traj,where='post')
        ax3.set_title('$u$ over time')
        ax3.set_xlim([0, self.t_f])
        ax3.set_ylim([-50, 50])
        fig.tight_layout()

        plt.show()




if __name__ =="__main__":
    planning_obj = STL_CBF_Planning()

    planning_obj.initialize_model()
    # Add STL constraint
    planning_obj.add_STL_constr()
    #planning_obj.always()

    # Add CBF constraint
    #planning_obj.add_CBF_constraint()
    # Solve and save the MIQP model
    planning_obj.save_model()

    start_time = time.time()
    status, x1_traj_discrete, x2_traj_discrete, u_traj_discrete = planning_obj.solve_MILP()
    elapsed_time = time.time() - start_time

    if status == GRB.Status.OPTIMAL:
        #planning_obj.plot_signals(x1_traj_discrete, x2_traj_discrete, u_traj_discrete)

        # Plot continous trajectory
        sim_obj = Dynamics_sim(planning_obj.t_f, planning_obj.N,planning_obj.X_0,u_traj_discrete,x1_traj_discrete,x2_traj_discrete,planning_obj.tau_list,planning_obj.x2_max,planning_obj.x2_min)
        sim_obj.generate_traj_plots()

        print("The MIQP is solved in ",elapsed_time)

        # Plot continous CBF margin
