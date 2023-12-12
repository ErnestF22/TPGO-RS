#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from dynamics_simulation_example4 import Dynamics_sim
import time

"""
Created on Sept 14 2:47:50 2019
@author: Guang Yang

Eventually Always Example for ACC 2019
"""

class STL_CBF_Planning:
    '''
    The Problem class contains problem parameters
    '''
    def __init__(self):
        self.model = Model("2d Double Integrators")
        self.M = 99999 #Big M
        self.N = 10 # Number of horizon steps
        self.t_f = 5 # End of STL time interval
        self.tau_list = [self.t_f/(self.N-1) for i in range(self.N-1)] # List that stores problem parameters
        self.k_cbf= [30,30] # CBF coefficients for x1 and x3 Positions
        self.u_min = -5
        self.u_max = 5
        self.num_state = 4 # 2-D Double Intergrator with positions on x1 and x3
        self.num_control = 2 # two states for control inputs
        self.num_jordan_block = 2
        self.size_jordan_block = 2 #s(i)
        self.X_0 = [1, 0, 1, 0] #initial state
        self.num_predicates = 4
        self.total_num_STL_predicates = self.num_predicates*self.N #The STL specification is valid for all time
        self.num_CBF_constraints = 1

    def initialize_model(self):
        '''
        OUTPUT: model with state variables, system Dynamics constraints and cost function
        The function takes in Gurobi model and generate state variables (x,u) according to the horizon length

        Gurobi Variable Format:
        {(0, 0): <gurobi.Var C0>,   x1(t_0) position x1
         (0, 1): <gurobi.Var C1>,   x2(t_0) velocity x2
         (0, 2): <gurobi.Var C0>,   x3(t_0) position x3
         (0, 3): <gurobi.Var C1>,   x4(t_0) velocity x4
         (1, 0): <gurobi.Var C2>,   x1(t_1)
         (1, 1): <gurobi.Var C3>,   x2(t_1)
         (1, 2): <gurobi.Var C0>,   x3(t_1)
         (1, 3): <gurobi.Var C1>,   x4(t_1)

         ....
         (0, 0):   u1(t_0) acceleration in x1
         (0, 1):   u2(t_0) accleration in x3
        '''
        # Initialize State Variables, Total State Variables = 2N+(N-1) = 3N-1
        print("Initialize state variables, initial condition, integer variable and cost function")
        self.U = self.model.addVars(self.N-1,self.num_control,name="u", vtype=GRB.CONTINUOUS,lb=self.u_min,ub=self.u_max)
        self.X = self.model.addVars(self.N,self.num_state,name="x",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.z_STL = self.model.addVar(vtype=GRB.BINARY, name="z_STL") #Final STL indicator variable, z_STL = 1 <-> STL is satisfied
        self.Z = self.model.addVars(self.N,self.num_predicates, name="Z_integer",vtype=GRB.BINARY)
        self.mu = self.model.addVars(self.N,self.num_predicates, name="mu",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY) #Predicate

        self.Z_cbf_x = self.model.addVars(self.N,self.num_CBF_constraints, name="Z_cbf_x",vtype=GRB.BINARY)
        self.Z_cbf_u = self.model.addVars(self.N,self.num_CBF_constraints, name="Z_cbf_u",vtype=GRB.BINARY)

        # Cost Function = \intergral_{0,t_f} u^Tu
        self.model.setObjective(quicksum((self.U[k,0]*self.U[k,0]+self.U[k,1]*self.U[k,1])*self.tau_list[k]*0.5  for k in range(self.N-1)),GRB.MINIMIZE)

        # Initilaize Dynamics Constraint
        self.model.addConstr(self.X[0,0] == self.X_0[0],"initial_state_constr_x1")
        self.model.addConstr(self.X[0,1] == self.X_0[1],"initial_state_constr_x2")
        self.model.addConstr(self.X[0,2] == self.X_0[2],"initial_state_constr_x3")
        self.model.addConstr(self.X[0,3] == self.X_0[3],"initial_state_constr_x4")

        for i in range(self.N-1):
            self.model.addConstr(self.X[i,0] + self.X[i,1]*self.tau_list[i] + (self.tau_list[i]**2/2)*self.U[i,0] == self.X[i+1,0],"dynamics_constr_x1_step"+str(i))
            self.model.addConstr(self.X[i,1] + self.U[i,0]*self.tau_list[i] == self.X[i+1,1],"dynamics_constr_x2_step"+str(i))
            self.model.addConstr(self.X[i,2] + self.X[i,3]*self.tau_list[i] + (self.tau_list[i]**2/2)*self.U[i,1] == self.X[i+1,2],"dynamics_constr_x3_step"+str(i))
            self.model.addConstr(self.X[i,3] + self.U[i,1]*self.tau_list[i] == self.X[i+1,3],"dynamics_constr_x4_step"+str(i))

    def STL_predicates(self):
        # Predicate Encoding
        for t_step in range(1,self.N):
            for i in range(self.num_predicates):
                self.model.addConstr(-self.mu[t_step,i] <= self.M * (1-self.Z[t_step,i]))
                self.model.addConstr(self.mu[t_step,i] <= self.M * self.Z[t_step,i])

        # Robustness Encoding
        for t_step in range(0,self.N):
            self.model.addConstr(self.mu[t_step,0] == (self.X[t_step,0]-1.7),"Predicate_x0_2")
            self.model.addConstr(self.mu[t_step,1] == (self.X[t_step,2]-2.7),"Predicate_x3")
            self.model.addConstr(self.mu[t_step,2] == (-self.X[t_step,2]+0.3),"Predicate_x3")
            self.model.addConstr(self.mu[t_step,3] == (-self.X[t_step,0]-0.2),"Predicate_x0")

        for t_step in range(0,self.N):
            self.model.addConstr((self.Z[t_step,0] == 1) >> (self.Z[t_step,1] == 1))
            self.model.addConstr((self.Z[t_step,1] == 1) >> (self.Z[t_step,0] == 1))
            self.model.addConstr((self.Z[t_step,2] == 1) >> (self.Z[t_step,3] == 1))
            self.model.addConstr((self.Z[t_step,3] == 1) >> (self.Z[t_step,2] == 1))



        # Integer Encoding
        self.model.addConstr(self.z_STL <= quicksum(self.Z[t_step,0] for t_step in [1,2,3,4,5,6,7,8,9]),"Eventually_Integer_Constraint")
        self.model.addConstr(self.z_STL <= quicksum(self.Z[t_step,1] for t_step in [1,2,3,4,5,6,7,8,9]),"Eventually_Integer_Constraint")
        self.model.addConstr(self.z_STL <= quicksum(self.Z[t_step,2] for t_step in [1,2,3,4,5,6,7,8,9]),"Eventually_Integer_Constraint")
        self.model.addConstr(self.z_STL <= quicksum(self.Z[t_step,3] for t_step in [1,2,3,4,5,6,7,8,9]),"Eventually_Integer_Constraint")


        for t_step in [1,2,3,4,5,6,7,8,9]:
            self.model.addConstr(self.z_STL >= self.Z[t_step,0],"Eventually_Integer_Constraint")
            self.model.addConstr(self.z_STL >= self.Z[t_step,1],"Eventually_Integer_Constraint")
            self.model.addConstr(self.z_STL >= self.Z[t_step,2],"Eventually_Integer_Constraint")
            self.model.addConstr(self.z_STL >= self.Z[t_step,3],"Eventually_Integer_Constraint")


    def CBF_Constraints(self):

        self.beta_X_h1 = self.model.addVars(self.N,name="beta_x_list_h1",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U_h1 = self.model.addVars(self.N,name="beta_u_list_h1",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)

        self.beta_X_h2 = self.model.addVars(self.N,name="beta_x_list_h2",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U_h2 = self.model.addVars(self.N,name="beta_u_list_h2",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)

        self.beta_X_h3 = self.model.addVars(self.N,name="beta_x_list_h3",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U_h3 = self.model.addVars(self.N,name="beta_u_list_h3",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)

        self.beta_X_h4 = self.model.addVars(self.N,name="beta_x_list_h3",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)
        self.beta_U_h4 = self.model.addVars(self.N,name="beta_u_list_h3",vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=GRB.INFINITY)

        for t_step in range(0,self.N-1):
            self.model.addConstr((self.beta_X_h1[t_step] + self.beta_U_h1[t_step] == 0.0),"beta_constr_h1")
            self.model.addConstr((self.beta_X_h2[t_step] + self.beta_U_h2[t_step] == 0.0),"beta_constr_h2")
            self.model.addConstr((self.beta_X_h3[t_step] + self.beta_U_h3[t_step] == 0.0),"beta_constr_h3")
            self.model.addConstr((self.beta_X_h3[t_step] + self.beta_U_h4[t_step] == 0.0),"beta_constr_h4")


        # Integer Varibale for removing postive contributing term
        # Z_cbf_x[k,i,j], where k is the time step, i is the eigen value indicator, j is the jordan block indicator
        #self.Z_cbf_x = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_x_list",vtype=GRB.BINARY)
        #self.Z_cbf_u = self.model.addVars(self.N,self.num_jordan_block,self.size_jordan_block, name="Z_CBF_u_list",vtype=GRB.BINARY)
        self.z_h1 = self.model.addVars(self.N,name="Z_integer_h1",vtype=GRB.BINARY)
        self.z_h2 = self.model.addVars(self.N,name="Z_integer_h2",vtype=GRB.BINARY)
        self.z_h3 = self.model.addVars(self.N,name="Z_integer_h3",vtype=GRB.BINARY)
        self.z_h4 = self.model.addVars(self.N,name="Z_integer_h4",vtype=GRB.BINARY)

        for i in range(0,self.N):
            self.model.addConstr((self.z_h1[i] == 1) >> (self.z_h2[i] == 1))
            self.model.addConstr((self.z_h2[i] == 1) >> (self.z_h1[i] == 1))
        '''
         Constraint on x, zeta_x
         if: -k*x(t_k) < 0, then enforce -k*x(t_k) + alpha_x*(k*x2_max) >= 0

         if: -ut < 0, then enforce -u - k*u*tau - alpha_u*(k*x2_min) >= 0

         if: k*x(t_k) < 0, then enforce k*x(t_k) - alpha_x_2*(k*x2_min) >= 0

         if: ut < 0, then enforce u + k*u*tau - alpha_u_2*(k*x2_min) >= 0
        '''

        for t_step in range(0,self.N-1):
            # x0(t)-1.5 >= 0
            self.model.addConstr((self.z_h1[t_step] == 1) >> (self.k_cbf[0]*self.X[t_step,1]+self.k_cbf[1]*self.X[t_step,1]*self.tau_list[t_step] + self.k_cbf[1]*self.X[t_step,0]-1.5*self.k_cbf[1]+self.beta_X_h1[t_step] >= 0) ,"CBF_x_constr_step"+str(t_step))
            self.model.addConstr((self.z_h1[t_step] == 1) >> (self.U[t_step,0]+self.k_cbf[0]*self.U[t_step,0]*self.tau_list[t_step]+0.5*self.k_cbf[0]*self.U[t_step,0]*self.tau_list[t_step]*self.tau_list[t_step] +self.beta_U_h1[t_step] >= 0) ,"CBF_u_constr"+str(t_step))


            # x3(t)-2 >= 0
            self.model.addConstr((self.z_h2[t_step] == 1) >> (self.k_cbf[0]*self.X[t_step,2]+self.k_cbf[1]*self.X[t_step,2]*self.tau_list[t_step]+ self.k_cbf[1]*self.X[t_step,2]-2.5*self.k_cbf[1]+self.beta_X_h2[t_step] >= 0) ,"CBF_x_constr_step"+str(t_step))
            self.model.addConstr((self.z_h2[t_step] == 1) >> (self.U[t_step,1]+self.k_cbf[0]*self.U[t_step,1]*self.tau_list[t_step]+0.5*self.k_cbf[0]*self.U[t_step,1]*self.tau_list[t_step]*self.tau_list[t_step] +self.beta_U_h2[t_step] >= 0),"CBF_u_constr"+str(t_step))

            # x3 <= 1.0
            #self.model.addConstr((self.z_h3[t_step] == 1) >> (-self.k_cbf[0]*self.X[t_step,3]-self.k_cbf[1]*self.X[t_step,3]*self.tau_list[t_step]- self.k_cbf[1]*self.X[t_step,2]+1.0*self.k_cbf[1]+self.beta_X_h3[t_step] >= 0) ,"CBF_x_constr_step"+str(t_step))
            #self.model.addConstr((self.z_h3[t_step] == 1) >> (-self.U[t_step,1]-self.k_cbf[0]*self.U[t_step,1]*self.tau_list[t_step]-0.5*self.k_cbf[0]*self.U[t_step,1]*self.tau_list[t_step]*self.tau_list[t_step] +self.beta_U_h3[t_step] >= 0),"CBF_u_constr"+str(t_step))





        # Safety Encoding
        for t_step in range(0,self.N-1):
            # Integer encoding for u
            self.model.addConstr(-(self.k_cbf[0]*self.U[t_step,0]) <= self.M * (1-self.z_h1[t_step]))
            self.model.addConstr((self.k_cbf[0]*self.U[t_step,0]) <= self.M * self.z_h1[t_step])

            # Integer encoding for x
            self.model.addConstr((self.k_cbf[0]*self.X[t_step,1]) <= self.M * (1-self.z_h1[t_step]))
            self.model.addConstr(-self.k_cbf[0]*self.X[t_step,1] <= self.M * self.z_h1[t_step])

            # Integer encoding for u
            self.model.addConstr(-(self.k_cbf[0]*self.U[t_step,1]) <= self.M * (1-self.z_h2[t_step]))
            self.model.addConstr((self.k_cbf[0]*self.U[t_step,1]) <= self.M * self.z_h2[t_step])

            # Integer encoding for x
            self.model.addConstr((self.k_cbf[0]*self.X[t_step,1]) <= self.M * (1-self.z_h2[t_step]))
            self.model.addConstr(-self.k_cbf[0]*self.X[t_step,1] <= self.M * self.z_h2[t_step])



        #for t_step in range(0,self.N):
            #self.model.addConstr(self.z_STL <= (self.z_h1[t_step] + self.z_h3[t_step]))
            #self.model.addConstr(self.z_STL >= self.z_h1[t_step])
            #self.model.addConstr(self.z_STL >= self.z_h3[t_step])
            #self.model.addConstr(self.z_STL <= self.z_h1[t_step])
            #self.model.addConstr(self.z_STL <= self.z_h3[t_step])



    def save_model(self):
        self.model.write("example4_model.lp")

    def solve_MILP(self):
        x_traj = []
        u_traj = []
        x1_traj = []
        x2_traj = []
        x3_traj = []
        x4_traj = []

        self.model.addConstr(self.z_STL == 1)
        self.model.optimize()

        status = self.model.status
        if status == GRB.Status.OPTIMAL:
            for i in range(self.N):
                for j in range(self.num_state):
                    x_traj.append(self.X[i,j].x)


            for i in range(self.N-1):
                for j in range(self.num_control):
                    u_traj.append(self.U[i,j].x)

            u1_traj = u_traj[0::2]
            u2_traj = u_traj[1::2]

            x1_traj = x_traj[0::4]
            x2_traj = x_traj[1::4]
            x3_traj = x_traj[2::4]
            x4_traj = x_traj[3::4]

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

        return status,x1_traj, x2_traj, x3_traj,x4_traj,u1_traj,u2_traj


if __name__ =="__main__":
    planning_obj = STL_CBF_Planning()
    planning_obj.initialize_model()
    planning_obj.STL_predicates()
    planning_obj.CBF_Constraints()
    # Solve and save the MILP model
    planning_obj.save_model()

    start_time = time.time()
    status, x1_traj_discrete, x2_traj_discrete, x3_traj_discrete,x4_traj_discrete,u1_traj, u2_traj = planning_obj.solve_MILP()
    elapsed_time = time.time() - start_time
    print("The MIQP is solved in ",elapsed_time)

    sim_obj = Dynamics_sim(planning_obj.t_f, planning_obj.N,planning_obj.X_0,u1_traj,u2_traj,x1_traj_discrete,x2_traj_discrete,x3_traj_discrete,x4_traj_discrete,planning_obj.tau_list)
    sim_obj.generate_traj_plots()


    #fig, ax = plt.subplots()
    #ax.plot(x1_traj_discrete,x3_traj_discrete,'--')
    #ax.grid(True, which='both')

    #ax.axhline(y=0, color='k')
    #ax.axvline(x=0, color='k')
    #plt.show()
