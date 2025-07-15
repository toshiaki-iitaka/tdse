import numpy as np
from scipy.fftpack import dst 
import math
import matplotlib.pyplot as plt
import sys

# --- Parameters ---
PI = math.pi
NMAX = 100  # Maximum number of rooms/points

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
    # Define complex imaginary unit
    iu = 0.0 + 1.0j
    # Ensure M is within a valid range for the cosine argument to avoid division by zero or large values

    ## Get input from user
    #print('M ALPHA')
    #try:
    #    m_input, alpha_input = map(float, input().split())
    #    m = int(m_input)
    #    alpha = float(alpha_input)
    #except ValueError:
    #    print("Invalid input. Please enter two numbers separated by a space (e.g., '1 0.1').")
    #    return
    #
    #if m <= 0 or m > NMAX:
    #    print(f"Warning: M ({m}) should be between 1 and NMAX ({NMAX}) for stable calculation of DT.")

    # fix to M=1 and ALPHA=0.1
    m = 1
    alpha = 0.1
    
    # Calculate DT (time step)
    # Fortran: DT=ALPHA/(EMAX*COS(M*PI/(NMAX+1)))
    cos_term = np.cos(m * PI / (NMAX + 1))
    if cos_term == 0:
        print("Error: Division by zero in DT calculation (COS term is zero). Adjust M or NMAX.")
        return
    dt = alpha / (emax * cos_term)

    # Calculate MTIME (total time steps)
    mtime = (2 * PI / alpha) * 10
    mtime = min(int(mtime), 10000) # Ensure MTIME is an integer

    # --- Calculation of Hamiltonian ---
    # HA and HB arrays
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX - 1, dtype=np.float64)

    # dgree of the method
    #NDEG  = 1
    for NDEG in [1,2,4,6]:
        NDEGH = int(np.floor(NDEG/2))
        NDEG1 = NDEG+1
        for n in range(NMAX - 1):
            hb[n] = -0.5 / dx2

        # --- Setting Initial Wave Function ---
        ph0, em = eigen_subroutine(m, dx, NMAX)
        phi = np.zeros((NDEG1,NMAX), dtype=np.complex128)

        for itime in range(-NDEG+1,1):
            phi[np.mod(itime,NDEG1)] = ph0 * np.exp(-iu * em * dt *itime)
            
        # --- Time Evolution Loop ---
        # Open the output file
        try:
            with open('tst'+str(NDEG)+'.dat', 'w') as f_out:
                plot_y = [[] for _ in range(6) ]
                for itime in range(mtime):
                    t = dt * itime
                    anorm = np.linalg.norm(phi[np.mod(itime,NDEG1)])**2
                    ovr = np.vdot(phi[np.mod(itime,NDEG1)], ph0)  * np.exp(-iu * em * t)
                    
                    
                    if(NDEG == 1):
                        epsilon_theory = abs(em*dt)**2 / 2.0 
                        q = phi[np.mod(itime,NDEG1)] /2
                    elif(NDEG == 2):
                        epsilon_theory = abs(em*dt)**3 / 6.0 
                        q = phi[np.mod(itime,NDEG1)]
                    elif(NDEG == 4):
                        epsilon_theory = abs(em*dt)**5 * 7.0 / 90.0
                        q = 2 * (-1.0/3.0 * phi[np.mod(itime-(NDEGH-1),NDEG1)]
                                 + 2.0/3.0 * ( phi[np.mod(itime-(NDEGH-1)+1,NDEG1)] + phi[np.mod(itime-(NDEGH-1)-1,NDEG1)]) )
                    elif(NDEG == 6):
                        epsilon_theory = abs(em*dt)**7 * 41.0 / 840.0
                        q = 3 * ( 13.0/10.0  *   phi[np.mod(itime-(NDEGH-1),NDEG1)]
                                  - 7.0/10.0  * ( phi[np.mod(itime-(NDEGH-1)+1,NDEG1)] + phi[np.mod(itime-(NDEGH-1)-1,NDEG1)]) 
                                  + 11.0/20.0 * ( phi[np.mod(itime-(NDEGH-1)+2,NDEG1)] + phi[np.mod(itime-(NDEGH-1)-2,NDEG1)])  )
                    else:
                        print(f'wrong NDEG = {NDEG}')
                        sys.exit()
                    #phi[np.mod(itime+1,NDEG1)] = -2 * iu * dt * H(q) + phi[np.mod(itime-(NDEGH+NDEGH-1),NDEG1)]
                    phi[np.mod(itime+1,NDEG1)] = -2 * iu * dt * H(q) + phi[np.mod(itime-(NDEG-1),NDEG1)]
                    
                    
                    file_output = (
                        f"{itime:4d} {emax*t:15.7E} {abs(anorm-1):15.7E} "
                        f"{abs(1.0 - abs(ovr)):15.7E} {abs(np.arctan2(ovr.imag, ovr.real)):15.7E} "
                        f"{epsilon_theory * itime:15.7E}"
                    )
                    print(file_output)
                    f_out.write(file_output + "\n")
                    
                    plot_y[0].append(itime)
                    plot_y[1].append(emax*t)
                    plot_y[2].append(abs(anorm-1))
                    plot_y[3].append(abs(1.0 - abs(ovr)))
                    plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))
                    plot_y[5].append(epsilon_theory * itime)

                #
                # plot
                #
                plt.gca().clear()
                plt.title('Error MSD'+str(NDEG))
                plt.xlabel(r'Emax*t')
                plt.ylabel(r'error')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlim([1,100])
                plt.ylim([1e-13,1e2])
                plt.plot(plot_y[1],plot_y[2], marker='o',    linestyle='None',  label='ε_norm')
                plt.plot(plot_y[1],plot_y[4], marker='o',    linestyle='None',  label='ε_phase')
                plt.plot(plot_y[1],plot_y[5], marker='None', linestyle='dashed',label='ε_theory')
                plt.legend()
                plt.show()
                plt.savefig("error"+str(NDEG)+".png")
                
        except Exception as e:
            print(f"An error occurred during file operations: {e}")
        
        print("Program finished. Results written to error"+str(NDEG)+".dat")

# Call the main program function
if __name__ == '__main__':
    nroom_program()
