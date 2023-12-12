from casadi import *
import matplotlib.pyplot as plt
import numpy as np
# -*- coding: utf-8 -*-

class Dynamics_sim:
    '''
    This Class takes in a set of discrete input (ZOH) and generate
    state trajectories based on continous system dynamics. Note the
    Update Interval in the simulation is based on real system dynamics that might
    not be the same in the planned interval
    '''
    def __init__(self,time_horizon,control_interval,X_0,control,tau_list):
        self.X_tk = X_0
        self.N = 12
        self.u_traj = control
        self.x1_traj = [X_0[0]]
        self.x2_traj = [X_0[1]]
        self.simulation_interval = 100
        self.t_k = 0
        self.t_f = time_horizon
        self.t_k_sim = 0
        self.t_f_sim = tau_list[0]
        self.x1 = MX.sym('x1')
        self.x2 = MX.sym('x2')
        self.x = vertcat(self.x1, self.x2)
        self.u = MX.sym('u')
        self.xdot = vertcat(self.x2, self.u)
        self.k_cbf = 5.0

        for t_step in range(0,control_interval-1):
            tau_sim = [(self.t_f_sim-self.t_k_sim)/(self.simulation_interval-1) for i in range(self.simulation_interval-1)]
            for sim_step in range(0,self.simulation_interval-1):
                self.dae = {'x':self.x, 'p':self.u, 'ode':self.xdot}
                #self.opts = {'t0': self.t_k_sim,'tf':self.t_f_sim}
                self.opts = {'t0': 0,'tf':tau_sim[sim_step]}
                self.F = integrator('F', 'cvodes', self.dae, self.opts)
                self.Fk = self.F(x0=self.X_tk,p=self.u_traj[t_step])
                self.t_f_sim += tau_sim[sim_step]
                self.x1_traj.append(self.Fk['xf'][0])
                self.x2_traj.append(self.Fk['xf'][1])
                self.X_tk = self.Fk['xf']

            self.t_k_sim += tau_list[t_step]
            self.t_f_sim = self.t_k_sim + tau_list[t_step]

    def cbf(self,max_min_indicator,u,x2):
        if max_min_indicator == 'min':
            return u + self.k_cbf*(x2-self.x2_min)

        if max_min_indicator == 'max':
            return -u + self.k_cbf*(-x2+self.x2_max)

    #def generate_cbf_plots(self):



    def generate_traj_plots(self):
        t = np.linspace(0,self.t_f,(self.N-1)*(self.simulation_interval-1)+1)
        t_discrete = np.linspace(0,self.t_f,self.N)

        spec_y_min = np.ones(len(t))*3

        x_2_min_t = [0.63, 0.63, 0.8, 0.8]
        x_2_min_y = [-100.0, 3, 3, -100.0]

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharey=False)
        ax1.plot(t, self.x1_traj)
        #ax1.plot(t, self.x1_traj_discrete, 'bo')
        ax1.set_title('State Trajectory for $x_1$')
        ax1.set_xlim([0, self.t_f])
        ax1.set_ylim([-5, 5])
        ax1.set(xlabel="Time (s)",ylabel="$x_1(t)$")

        ax2.plot(t,self.x2_traj)
        #ax2.plot(t,self.x2_traj_discrete, 'bo')
        ax2.plot(t,spec_y_min, 'k--')
        ax2.set_title('State Trajectory for $x_2$')
        ax2.set_xlim([0, self.t_f])
        ax2.fill(x_2_min_t,x_2_min_y,'r')
        ax2.set_ylim([-6, 6])
        ax2.set(xlabel="Time (s)",ylabel="$x_2(t)$")

        ax3.step(t_discrete[0:-1],self.u_traj)
        ax3.set_title('Control Input $u$')
        ax3.set_xlim([0, self.t_f])
        ax3.set_ylim([-50, 50])
        ax3.set(xlabel="Time (s)",ylabel="$u(t)$")


        fig.subplots_adjust(hspace=0.5)

        fig, ax = plt.subplots(1, 1, sharey=False)
        ax.plot(self.x1_traj, self.x2_traj)
        plt.show()


if __name__ == "__main__":
    test_t_f = 1
    test_N = 5
    tau_list = [test_t_f/(test_N-1) for i in range(test_N-1)]
    X_0 = [1,0]
    u = [10,-10,10,-10]
    #u=2
    test_obj = Dynamics_sim(test_t_f,test_N,X_0,u,tau_list)
    test_obj.generate_traj_plots()
