import numpy as np
import math
import matplotlib.pyplot as plt

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

    # Calculate initial wave function values based on sine function
    # Fortran N is 1-indexed, Python n is 0-indexed.
    # So, Fortran N becomes (n_idx + 1) in Python for calculations involving N.
    for n_idx in range(nmax_val):
        ph0[n_idx] = np.sin((PI * m_val * (n_idx + 1)) / (nmax_val + 1))

    # Normalization of wave function
    pa = 0.0
    for n_idx in range(nmax_val):
        pa += abs(ph0[n_idx])**2
    pa = np.sqrt(pa)

    for n_idx in range(nmax_val):
        ph0[n_idx] = ph0[n_idx] / pa

    # Calculate eigenenergy (EM)
    # The active line in Fortran: EM=2*(-0.5/DX/DX) * COS((PI*M)/(NMAX+1))
    em = 2 * (-0.5 / (dx_val * dx_val)) * np.cos((PI * m_val) / (nmax_val + 1))

    return ph0, em

# --- Main Program: TSTCN (Crank-Nicolson Scheme) ---
def tstcn_program():
    """
    Main program to simulate 1D quantum mechanics using the Crank-Nicolson scheme.
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # Get input from user
    print('M, ALPHA')
    try:
        m_input, alpha_input = map(float, input().split())
        m = int(m_input)
        alpha = float(alpha_input)
    except ValueError:
        print("Invalid input. Please enter two numbers separated by a space (e.g., '1 0.1').")
        return

    # --- Initial Data Calculations ---
    dx = 1.0
    dx2 = dx * dx
    emax = 1.0 / dx2 # EMAX from Fortran's initial calculation

    # Calculate DT (time step)
    # Fortran: DT=ALPHA/(EMAX*COS(M*PI/(NMAX+1)))
    cos_term = np.cos(m * PI / (NMAX + 1))
    if abs(cos_term) < 1e-9: # Check for near-zero to prevent division by zero
        print("Error: COS term in DT calculation is near zero. Adjust M or NMAX to avoid singularity.")
        return
    dt = alpha / (emax * cos_term)

    # Calculate MTIME (total time steps)
    # Fortran: MTIME=(2*PI/ALPHA) * 10
    # Fortran: MTIME=MIN(MTIME,10000)
    mtime = (2 * PI / alpha) * 10
    mtime = min(int(mtime), 10000) # Ensure MTIME is an integer

    # --- Calculation of Hamiltonian and Crank-Nicolson coefficients ---
    # Note: Fortran HB is DIMENSION HB(NMAX) here, unlike previous programs.
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX, dtype=np.float64) # HB has NMAX elements
    a1 = np.zeros(NMAX, dtype=np.complex128)
    b1 = np.zeros(NMAX, dtype=np.complex128)
    a2 = np.zeros(NMAX, dtype=np.complex128)
    b2 = np.zeros(NMAX, dtype=np.complex128)

    # DO 50 N=1,NMAX
    # C        HA(N) = V(X)+1/DX2 - E0
    # C        HA(N) = 1/DX2-E0
    # HA(N)=0
    # HB(N) = -0.5/DX2
    # A1(N) = 1.-0.5*IU*DT*HA(N)
    # A2(N) = 1.+0.5*IU*DT*HA(N)
    # B1(N) = -0.5*IU*DT*HB(N)
    # B2(N) = +0.5*IU*DT*HB(N)
    # 50 CONTINUE
    for n_idx in range(NMAX):
        ha[n_idx] = 0.0 # As per Fortran code
        hb[n_idx] = -0.5 / dx2

        a1[n_idx] = 1.0 - 0.5 * iu * dt * ha[n_idx]
        a2[n_idx] = 1.0 + 0.5 * iu * dt * ha[n_idx]
        b1[n_idx] = -0.5 * iu * dt * hb[n_idx]
        b2[n_idx] = +0.5 * iu * dt * hb[n_idx]

    # --- Setting Initial Wave Function ---
    ph0, em = eigen_subroutine(m, dx, NMAX)

    # Initialize P array with the initial wave function
    p = np.zeros(NMAX, dtype=np.complex128)
    # DO 60 N=1,NMAX
    # P(N)=PH0(N)
    # 60 CONTINUE
    for n_idx in range(NMAX):
        p[n_idx] = ph0[n_idx]

    # --- Time Evolution Loop ---
    # Open the output file
    try:
        with open('tstcn.dat', 'w') as f_out:
            plot_x = range(mtime + 1)
            plot_y = [[] for _ in range(6) ]
            # DO 1000 ITIME=0,MTIME
            for itime in range(mtime + 1): # Loop from 0 to MTIME inclusive
                t = itime * dt

                # --- Output Wave Function (every 10th step, starting from 1) ---
                # IF(MOD(ITIME,10).EQ. 1) THEN
                if (itime % 10) == 1:
                    anorm = 0.0
                    ovr = 0.0 + 0.0j
                    # DO 150 N=1,NMAX
                    # ANORM=ANORM+ABS(P(N))**2
                    # OVR=OVR+DCONJG(P(N))*PH0(N)*EXP(-IU*EM*T)
                    # 150   CONTINUE
                    for n_idx in range(NMAX):
                        anorm += abs(p[n_idx])**2
                        ovr += np.conjugate(p[n_idx]) * ph0[n_idx] * np.exp(-iu * em * t)

                    # Fortran console output:
                    # WRITE(*,100) ITIME,EMAX*T,ANORM-1,
                    # &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
                    # &,ABS((EM*DT)**3/12.) * ITIME
                    # Format: I4, 30E15.7 (for 6 values, 30E15.7 is generous)
                    error_term = abs((em * dt)**3 / 12.0) * itime
                    console_output = (
                        f"{itime:4d} {emax*t:15.7E} {anorm-1:15.7E} "
                        f"{1.0 - abs(ovr):15.7E} {np.arctan2(ovr.imag, ovr.real):15.7E} "
                        f"{error_term:15.7E}"
                    )
                    print(console_output)

                    # Fortran file output:
                    # WRITE(20,100) ITIME,EMAX*T,ABS(ANORM-1),
                    # &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
                    # &,ABS((EM*DT)**3/12.) * ITIME
                    file_output = (
                        f"{itime:4d} {emax*t:15.7E} {abs(anorm-1):15.7E} "
                        f"{abs(1.0 - abs(ovr)):15.7E} {abs(np.arctan2(ovr.imag, ovr.real)):15.7E} "
                        f"{error_term:15.7E}"
                    )
                    f_out.write(file_output + "\n")
                    
                    plot_y[0].append(itime)
                    plot_y[1].append(emax*t)
                    plot_y[2].append(abs(anorm-1))
                    plot_y[3].append(abs(1.0 - abs(ovr)))
                    plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))
                    plot_y[5].append(abs(em*dt)**3 * 1.0/12.0 * itime)

                # --- Multiplication of Symmetric Tridiagonal Matrix: O = (AB)*P ---
                # This calculates the right-hand side vector for the Crank-Nicolson equation.
                # The matrix (AB) here uses A1 (diagonal) and B1 (off-diagonal).
                o = np.zeros(NMAX, dtype=np.complex128)

                # O(1) = A1(1)*P(1)+B1(1)*P(2)
                o[0] = a1[0] * p[0] + b1[0] * p[1]
                # DO 70 N=2,NMAX-1
                # O(N) = B1(N-1)*P(N-1)+A1(N)*P(N)+B1(N)*P(N+1)
                # 70 CONTINUE
                for n_idx in range(1, NMAX - 1): # Python index 1 to NMAX-2
                    o[n_idx] = b1[n_idx-1] * p[n_idx-1] + a1[n_idx] * p[n_idx] + b1[n_idx] * p[n_idx+1]
                # O(NMAX) = B1(NMAX-1)*P(NMAX-1)+A1(NMAX)*P(NMAX)
                o[NMAX-1] = b1[NMAX-2] * p[NMAX-2] + a1[NMAX-1] * p[NMAX-1]

                # --- Inversion of Symmetric Tridiagonal Matrix: P = O / (AB) ---
                # This solves the linear system using a variant of the Thomas algorithm.
                # The matrix (AB) here uses A2 (diagonal) and B2 (off-diagonal).
                e_solver = np.zeros(NMAX, dtype=np.complex128)
                f_solver = np.zeros(NMAX, dtype=np.complex128)

                # E(1) = -A2(1)/B2(1)
                # F(1) = O(1)/B2(1)
                # Handle potential division by zero for B2[0]
                if abs(b2[0]) < 1e-18:
                    print(f"Error: B2[0] is near zero at ITIME={itime}. Cannot proceed with tridiagonal solver.")
                    return
                e_solver[0] = -a2[0] / b2[0]
                f_solver[0] = o[0] / b2[0]

                # DO 80 N=2,NMAX
                # E(N)= -(B2(N-1)/E(N-1)+A2(N)) / B2(N)
                # F(N)= O(N)/B2(N) + (B2(N-1)/B2(N))*(F(N-1)/E(N-1))
                # 80 CONTINUE
                for n_idx in range(1, NMAX): # Python index 1 to NMAX-1
                    # Handle potential division by zero for E[n_idx-1] and B2[n_idx]
                    if abs(e_solver[n_idx-1]) < 1e-18:
                        print(f"Error: E_solver[{n_idx-1}] is near zero at ITIME={itime}. Cannot proceed with tridiagonal solver.")
                        return
                    if abs(b2[n_idx]) < 1e-18:
                        print(f"Error: B2[{n_idx}] is near zero at ITIME={itime}. Cannot proceed with tridiagonal solver.")
                        return

                    e_solver[n_idx] = -(b2[n_idx-1] / e_solver[n_idx-1] + a2[n_idx]) / b2[n_idx]
                    f_solver[n_idx] = o[n_idx] / b2[n_idx] + (b2[n_idx-1] / b2[n_idx]) * (f_solver[n_idx-1] / e_solver[n_idx-1])

                # P(NMAX)=-F(NMAX)/E(NMAX)
                # Handle potential division by zero for E[NMAX-1]
                if abs(e_solver[NMAX-1]) < 1e-18:
                    print(f"Error: E_solver[{NMAX-1}] is near zero at ITIME={itime}. Cannot proceed with tridiagonal solver.")
                    return
                p[NMAX-1] = -f_solver[NMAX-1] / e_solver[NMAX-1]

                # DO 90 N=NMAX-1,1,-1
                # P(N)=(P(N+1)-F(N))/E(N)
                # 90 CONTINUE
                for n_idx in range(NMAX - 2, -1, -1): # Python index NMAX-2 down to 0
                    # Handle potential division by zero for E[n_idx]
                    if abs(e_solver[n_idx]) < 1e-18:
                        print(f"Error: E_solver[{n_idx}] is near zero at ITIME={itime}. Cannot proceed with tridiagonal solver.")
                        return
                    p[n_idx] = (p[n_idx+1] - f_solver[n_idx]) / e_solver[n_idx]

            #
            # plot
            #
            plt.gca().clear()
            plt.title('Crank-Nicolson')
            plt.xlabel(r'Emax*t')
            plt.ylabel(r'error')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlim([1,100])
            plt.ylim([1e-20,1])
            plt.plot(plot_y[1],plot_y[2], marker='o',    linestyle='None',  label='ε_norm')
            plt.plot(plot_y[1],plot_y[4], marker='o',    linestyle='None',  label='ε_phase')
            plt.plot(plot_y[1],plot_y[5], marker='None', linestyle='dashed',label='ε_theory')
            plt.legend()
            plt.show()
            plt.savefig("btstcn.png")
                
    except Exception as e:
        print(f"An error occurred during file operations: {e}")

    print("Program finished. Results written to tstcn.dat")

# Call the main program function
if __name__ == '__main__':
    tstcn_program()
