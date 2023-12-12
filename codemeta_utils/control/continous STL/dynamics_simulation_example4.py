from casadi import *
import matplotlib.pyplot as plt
import numpy as np
# -*- coding: utf-8 -*-

class Dynamics_sim:
    '''
    This Class takes in a set of discrete input (ZOH) and generate
    state trajectories based on continous system dynamics
    '''
    def __init__(self,time_horizon,control_interval,X_0,u1,u2,x1_traj_discrete,x2_traj_discrete,x3_traj_discrete,x4_traj_discrete,tau_list):
        self.X_tk = X_0
        self.N = control_interval
        self.x1_traj_discrete = x1_traj_discrete
        self.x2_traj_discrete = x2_traj_discrete
        self.x3_traj_discrete = x3_traj_discrete
        self.x4_traj_discrete = x4_traj_discrete
        self.u1_traj = u1
        self.u2_traj = u2
        self.x1_traj = [X_0[0]]
        self.x2_traj = [X_0[1]]
        self.x3_traj = [X_0[2]]
        self.x4_traj = [X_0[3]]
        self.simulation_interval = 30
        self.t_k = 0
        self.t_f = time_horizon
        self.t_k_sim = 0
        self.t_f_sim = tau_list[0]
        self.x1 = MX.sym('x1')
        self.x2 = MX.sym('x2')
        self.x3 = MX.sym('x3')
        self.x4 = MX.sym('x4')
        self.x = vertcat(self.x1, self.x2,self.x3,self.x4)
        self.u1 = MX.sym('u1')
        self.u2 = MX.sym('u2')
        self.u = vertcat(self.u1,self.u2)
        self.xdot = vertcat(self.x2, self.u1,self.x4,self.u2)
        self.k_cbf = 1.0

        for t_step in range(0,control_interval-1):
            tau_sim = [(self.t_f_sim-self.t_k_sim)/(self.simulation_interval-1) for i in range(self.simulation_interval-1)]
            for sim_step in range(0,self.simulation_interval-1):
                self.dae = {'x':self.x, 'p':self.u, 'ode':self.xdot}
                self.opts = {'t0': 0,'tf':tau_sim[sim_step]}
                self.F = integrator('F', 'cvodes', self.dae, self.opts)
                self.Fk = self.F(x0=self.X_tk,p=vertcat(self.u1_traj[t_step],self.u2_traj[t_step]))
                self.t_f_sim += tau_sim[sim_step]
                self.x1_traj.append(self.Fk['xf'][0])
                self.x2_traj.append(self.Fk['xf'][1])
                self.x3_traj.append(self.Fk['xf'][2])
                self.x4_traj.append(self.Fk['xf'][3])
                self.X_tk = self.Fk['xf']

                #self.t_k_sim += tau_sim[sim_step] #incrementally updating intial time

            self.t_k_sim += tau_list[t_step]
            self.t_f_sim = self.t_k_sim + tau_list[t_step]



    def generate_traj_plots(self):
        #t = np.linspace(0,self.t_f,(self.N-1)*(self.simulation_interval-1)+1)
        #t_discrete = np.linspace(0,self.t_f,self.N)

        # Eventually -3 and 3
        #spec_y_min = np.ones(len(t))*-2
        #spec_y_max = np.ones(len(t))*2

        #Fill CBF unsafe region
        #x_2_max_t = [0,0, self.t_f, self.t_f]
        #x_2_max_y = [self.x2_max, 100, 100, self.x2_max]

        #x_2_min_t = [0,0, self.t_f, self.t_f]
        #x_2_min_y = [self.x2_min, -100, -100, self.x2_min]

        #fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharey=False)
        #ax1.plot(t, self.x1_traj)
        #ax1.plot(t_discrete,self.x1_traj_discrete, 'bo')
        #ax1.plot(t,spec_y_min, 'k--')
        #ax1.plot(t,spec_y_max, 'k--')
        #ax1.set_title('State Trajectory for $x_1$')
        #ax1.set_xlim([0, self.t_f])
        #ax1.set_ylim([-5, 5])
        #ax1.set(xlabel="Time (s)",ylabel="$x_1(t)$")

        #ax2.plot(t,self.x3_traj)
        #ax2.plot(t_discrete,self.x3_traj_discrete, 'bo')
        #ax2.set_title('State Trajectory for $x_2$')
        #ax2.set_xlim([0, self.t_f])
        #ax2.set_ylim([-20, 20])
        #ax2.fill(x_2_max_t,x_2_max_y,'r')
        #ax2.fill(x_2_min_t,x_2_min_y,'r')
        #ax2.set(xlabel="Time (s)",ylabel="$x_2(t)$")

        #ax3.step(t_discrete[0:-1],self.u_traj)
        #ax3.set_title('Control Inputs $u$')
        #ax3.set_xlim([0, self.t_f])
        #ax3.set_ylim([-50, 50])
        #ax3.set(xlabel="Time (s)",ylabel="$u(t)$")


        #fig.subplots_adjust(hspace=0.5)
        #Unsafe Region
        x1_coord = [-10,-10, 0, 0]
        x3_coord = [-10,0, 0, -10]

        fig, ax = plt.subplots()
        ax.plot(self.x1_traj,self.x3_traj,'b')
        ax.plot(self.x1_traj_discrete,self.x3_traj_discrete,'bo')
        ax.grid(True, which='both')
        ax.set_xlabel('$x_1$ Position',fontsize=15)
        ax.set_ylabel('$x_3$ Position',fontsize=15)
        ax.fill(x1_coord,x3_coord,'r')
        ax.set_xlim([-0.6,3.0])
        ax.set_ylim([-0.6,3.0])
        ax.set_title('2-D Double Integrator System with the CBF Constraints',fontsize=18)
        ax.set_aspect(aspect=1.0)


        ax.axhline(y=0, color='k')
        ax.axvline(x=0, color='k')
        plt.show()
