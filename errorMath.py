import numpy as np
from scipy import integrate
from matplotlib.pylab import *
 
def diffEq(t, y):
    global f, g
    # Output from ODE function must be a COLUMN vector, with n rows
    low = y[0]
    high = y[1]

    n = len(y)      # 1: implies its a single ODE
    dydt = np.zeros((n,1))
    dydt[0] = -(1-f)*low + f*high + g*low
    dydt[1] = +(1-f)*low - f*high
    return dydt

if __name__ == '__main__':
    global f,g
    gI, fI = -1,-1
    numFact = 10
    arrayErrors = np.zeros((numFact,numFact))
    gFactors = np.linspace(0,0.1,numFact)
    fFactors = 1. - np.logspace(-5,0,numFact)
    for g in gFactors:
        print gI
        gI += 1
        fI = -1
        for f in fFactors:
            fI += 1
            # Start by specifying the integrator:
            # use ``vode`` with "backward differentiation formula"
            r = integrate.ode(diffEq).set_integrator('vode', method='bdf')

            # Set the time range
            t_start = 0.0
            t_final = 10.0
            delta_t = 0.01
            # Number of time steps: 1 extra for initial condition
            num_steps = np.floor((t_final - t_start)/delta_t) + 1

            # Set initial condition(s): for integrating variable and time!
            prot = [1,1]
            r.set_initial_value(prot, t_start)

            # Additional Python step: create vectors to store trajectories
            t = np.zeros((num_steps, 1))
            CA = np.zeros((num_steps, 2))
            t[0] = t_start
            CA[0,:] = prot

            # Integrate the ODE(s) across each delta_t timestep
            k = 1
            while r.successful() and k < num_steps:
                r.integrate(r.t + delta_t)

                # Store the results to plot later
                t[k] = r.t
                CA[k,:] = r.y[0:2]
                k += 1

            # All done!  Plot the trajectories:
            # plot(t, CA[:,0]/(CA[:,0]+CA[:,1]))
            error = (CA[-1,0]/(CA[-1,0]+CA[-1,1]) - f)
            arrayErrors[fI,gI] = error


    #grid('on')
    #xlabel('Time [minutes]')
    #ylabel('Fraction')
    #ylim((0,1))
    plt.imshow(arrayErrors,interpolation='none',cmap='Blues')
    plt.xticks(np.arange(numFact),gFactors,rotation=90)
    plt.yticks(np.arange(numFact),fFactors)
    colorbar()
    ylabel("Fraction of population")
    xlabel("Difference in growth")
    show()