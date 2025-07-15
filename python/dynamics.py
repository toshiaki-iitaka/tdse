import numpy as np
from scipy.fftpack import dst 
import math
import matplotlib.pyplot as plt
import sys

# --- Parameters ---
PI = math.pi
NMAX = 100  # Maximum number of rooms/points
#MTIME = 1500 # Total number of time steps for evolution loop
MTIME = 10000 # Total number of time steps for evolution loop
iu = 0.0 + 1.0j

# --- Subroutine EIGEN: Sets initial wave function and calculates initial energy ---
def eigen_subroutine(m_val, dx_val, nmax_val):
    """
    Calculates the initial wave function (PH0) and its corresponding energy (EM).

    Args:
        m_val (int): Mode number for the initial wave function.
        dx_val (float): Spatial step size.
        nmax_val (int): Maximum number of points (NMAX).

    Returns:
        tuple: A tuple containing:
            - ph0 (np.ndarray): Initial wave function (complex array).
            - em (float): Eigenenergy.
    """
    ph0 = np.zeros(nmax_val, dtype=np.complex128)

    ph0[m_val-1]=1
    ph0=dst(ph0,type=1,norm='ortho')
    em = 2 * (-0.5 / (dx_val * dx_val)) * np.cos((PI * m_val) / (nmax_val + 1))

    return ph0, em

# --- Parameters --- Operate Hamiltonian on quantum state p
def H(p):
    q=np.empty_like(p)
    q[0] = ha[0] * p[0] + hb[0] * p[1] 
    for n in range(1, NMAX - 1):
        q[n] = hb[n-1] * p[n-1] + ha[n] * p[n] + hb[n] * p[n+1] 
    q[NMAX-1] = hb[NMAX-2] * p[NMAX-2] + ha[NMAX-1] * p[NMAX-1]
    return q

# --- Time evolution with Euler method ----
def MSD1(p,dt):
    global iu
    return p -iu * dt * H(p)

# --- Main Program: NROOM ---
def nroom_program():
    """
    Main program to simulate N-room quantum mechanics.
    """

    global ha,hb,dt,iu

    # --- Initial Data Calculations ---
    dx = 1.0
    dx2 = dx * dx
    emax = 1.0 / dx2

    m = 1
    alpha = 0.1
    
    dt = alpha / emax

    # --- Calculation of Hamiltonian ---
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX - 1, dtype=np.float64)
    for n in range(NMAX - 1):
        hb[n] = -0.5 / dx2

    # plot emvironment
    plt.gca().clear()
    plt.title('Dynamics MSD')
    plt.xlabel(r'n')
    plt.ylabel(r'prob')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlim([0,NMAX])
    plt.ylim([0,20/NMAX])

    # dgree of the method
    for NDEG in [1,2,4,6]:
        NDEGH = int(np.floor(NDEG/2))
        NDEG1 = NDEG+1

        # --- Setting Initial Wave Function ---
        phi = np.zeros((NDEG1,NMAX), dtype=np.complex128)
        
        x0 = 0.3 * dx * NMAX
        sg = 0.05 * dx * NMAX
        pa_momentum = 0.1 * 2 * PI / dx
        for n_idx in range(NMAX):
            x = (n_idx + 1) * dx 
            phi[np.mod(0,NDEG1),n_idx] = np.exp(-0.5 * (x - x0)**2 / sg**2) * np.exp(iu * pa_momentum * (x - x0))
        anorm = np.linalg.norm(phi[np.mod(0,NDEG1),:])
        phi[np.mod(0,NDEG1),:] = phi[np.mod(0,NDEG1),:] / anorm

        # -- prepare wave function for t<0
        for itime in range(1,NDEG+1):
            phi[np.mod(-itime,NDEG1),:] = MSD1( phi[np.mod(-itime+1,NDEG1),:] , -dt)
            anorm = np.linalg.norm(phi[np.mod(-itime,NDEG1),:])
            phi[np.mod(-itime,NDEG1),:] = phi[np.mod(-itime,NDEG1),:]

        #ph0, em = eigen_subroutine(m, dx, NMAX)
            
        # --- Time Evolution Loop ---
        try:
            for itime in range(MTIME):
                t = dt * itime
                anorm = np.linalg.norm(phi[np.mod(itime,NDEG1)])**2
                #ovr = np.vdot(phi[np.mod(itime,NDEG1)], ph0)  * np.exp(-iu * em * t)
                
                if(NDEG == 1):
                    #epsilon_theory = abs(em*dt)**2 / 2.0 
                    q = phi[np.mod(itime,NDEG1)] /2
                elif(NDEG == 2):
                    #epsilon_theory = abs(em*dt)**3 / 6.0 
                    q = phi[np.mod(itime,NDEG1)]
                elif(NDEG == 4):
                    #epsilon_theory = abs(em*dt)**5 * 7.0 / 90.0
                    q = 2 * (-1.0/3.0 * phi[np.mod(itime-(NDEGH-1),NDEG1)]
                             + 2.0/3.0 * ( phi[np.mod(itime-(NDEGH-1)+1,NDEG1)] + phi[np.mod(itime-(NDEGH-1)-1,NDEG1)]) )
                elif(NDEG == 6):
                    #epsilon_theory = abs(em*dt)**7 * 41.0 / 840.0
                    q = 3 * ( 13.0/10.0  *   phi[np.mod(itime-(NDEGH-1),NDEG1)]
                              - 7.0/10.0  * ( phi[np.mod(itime-(NDEGH-1)+1,NDEG1)] + phi[np.mod(itime-(NDEGH-1)-1,NDEG1)]) 
                              + 11.0/20.0 * ( phi[np.mod(itime-(NDEGH-1)+2,NDEG1)] + phi[np.mod(itime-(NDEGH-1)-2,NDEG1)])  )
                else:
                    print(f'wrong NDEG = {NDEG}')
                    sys.exit()

                phi[np.mod(itime+1,NDEG1)] = -2 * iu * dt * H(q) + phi[np.mod(itime-(NDEG-1),NDEG1)]
                    
        except Exception as e:
            print(f"An error occurred during file operations: {e}")
                    
        plot_y = []
        for n_idx in range(NMAX):
            plot_y.append(np.abs(phi[np.mod(MTIME,NDEG1), n_idx])**2)
        plt.plot(plot_y, linestyle='solid', label='MSD'+str(NDEG))
    #
    # plot
    #
    plt.legend()
    plt.show()
    plt.savefig("dynamics.png")
        
    print("Program finished.")

# Call the main program function
if __name__ == '__main__':
    nroom_program()
