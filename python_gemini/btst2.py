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
    for n in range(nmax_val):
        # Fortran N is 1-indexed, Python n is 0-indexed.
        # So, Fortran N becomes (n + 1) in Python.
        ph0[n] = np.sin((PI * m_val * (n + 1)) / (nmax_val + 1))

    # Normalization of wave function
    pa = 0.0
    for n in range(nmax_val):
        pa += abs(ph0[n])**2
    pa = np.sqrt(pa)

    for n in range(nmax_val):
        ph0[n] = ph0[n] / pa

    # Calculate eigenenergy (EM)
    # The commented out line in Fortran: EM=-4*(-0.5/DX/DX) * SIN((PI*M)/(NMAX+1)/2.0)**2
    # The active line in Fortran: EM=2*(-0.5/DX/DX) * COS((PI*M)/(NMAX+1))
    em = 2 * (-0.5 / (dx_val * dx_val)) * np.cos((PI * m_val) / (nmax_val + 1))

    return ph0, em

# --- Main Program: NROOM ---
def nroom_program():
    """
    Main program to simulate N-room quantum mechanics.
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
    emax = 1.0 / dx2

    # Ensure M is within a valid range for the cosine argument to avoid division by zero or large values
    if m <= 0 or m > NMAX:
        print(f"Warning: M ({m}) should be between 1 and NMAX ({NMAX}) for stable calculation of DT.")
        # Attempt to proceed with a default or adjusted M if necessary, or exit.
        # For now, let's just warn and continue, as Fortran might not have explicit checks.

    # Calculate DT (time step)
    # Fortran: DT=ALPHA/(EMAX*COS(M*PI/(NMAX+1)))
    cos_term = np.cos(m * PI / (NMAX + 1))
    if cos_term == 0:
        print("Error: Division by zero in DT calculation (COS term is zero). Adjust M or NMAX.")
        return
    dt = alpha / (emax * cos_term)

    # Calculate MTIME (total time steps)
    # Fortran: MTIME=(2*PI/ALPHA) * 10
    # Fortran: MTIME=MIN(MTIME,10000)
    mtime = (2 * PI / alpha) * 10
    mtime = min(int(mtime), 10000) # Ensure MTIME is an integer

    # --- Calculation of Hamiltonian ---
    # HA and HB arrays
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX - 1, dtype=np.float64)

    # Fortran code had V(N) and part of HA(N) commented out.
    # HA(N) is effectively 0 based on the provided Fortran.
    # DO 50 N=1,NMAX
    # C        V(N)=0
    # C        HA(N) = V(N) + 1/DX2
    # HA(N)=0
    # 50 CONTINUE
    # -> ha remains all zeros as initialized

    # DO 51 N=1,NMAX-1
    #     HB(N) = -0.5/DX2
    # 51 CONTINUE
    for n in range(NMAX - 1):
        hb[n] = -0.5 / dx2

    # --- Setting Initial Wave Function ---
    ph0, em = eigen_subroutine(m, dx, NMAX)

    # Initialize P and Q arrays
    p = np.zeros(NMAX, dtype=np.complex128)
    q = np.zeros(NMAX, dtype=np.complex128)

    # DO 80 N=1,NMAX
    # P(N)=PH0(N)
    # Q(N)=PH0(N)*EXP(-IU*EM*DT)
    # 80 CONTINUE
    for n in range(NMAX):
        p[n] = ph0[n]
        q[n] = ph0[n] * np.exp(-iu * em * dt)

    # --- Time Evolution Loop ---
    # Open the output file
    try:
        with open('tst2.dat', 'w') as f_out:
            plot_x = range(1, mtime + 1, 3)
            plot_y = [[] for _ in range(6) ]
            # DO 1000 ITIME=1,MTIME,3
            for itime in range(1, mtime + 1, 3):
                t = dt * itime

                # --- Output Wave Function (every 3rd step, starting from 1) ---
                # IF(MOD(ITIME,3).EQ. 1) THEN
                if (itime % 3) == 1:
                    anorm = 0.0
                    ovr = 0.0 + 0.0j
                    # DO 150 N=1,NMAX
                    # ANORM=ANORM+ABS(Q(N))**2
                    # OVR=OVR+DCONJG(Q(N))*PH0(N)*EXP(-IU*EM*T)
                    # 150   CONTINUE
                    for n in range(NMAX):
                        anorm += abs(q[n])**2
                        ovr += np.conjugate(q[n]) * ph0[n] * np.exp(-iu * em * t)

                    # Fortran console output:
                    # WRITE(*,100) ITIME,EMAX*T,ANORM-1,
                    # &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
                    # &,ABS(EM*DT)**3/6.* ITIME
                    # Format: I4, 30E15.7 (for 6 values, the 30E15.7 is generous)
                    console_output = (
                        f"{itime:4d} {emax*t:15.7E} {anorm-1:15.7E} "
                        f"{1.0 - abs(ovr):15.7E} {np.arctan2(ovr.imag, ovr.real):15.7E} "
                        f"{abs(em*dt)**3 / 6.0 * itime:15.7E}"
                    )
                    print(console_output)

                    # Fortran file output:
                    # WRITE(20,100) ITIME,EMAX*T,ABS(ANORM-1),
                    # &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
                    # &,ABS(EM*DT)**3/6. * ITIME
                    file_output = (
                        f"{itime:4d} {emax*t:15.7E} {abs(anorm-1):15.7E} "
                        f"{abs(1.0 - abs(ovr)):15.7E} {abs(np.arctan2(ovr.imag, ovr.real)):15.7E} "
                        f"{abs(em*dt)**3 / 6.0 * itime:15.7E}"
                    )
                    f_out.write(file_output + "\n")

                    plot_y[0].append(itime)
                    plot_y[1].append(emax*t)
                    plot_y[2].append(abs(anorm-1))
                    plot_y[3].append(abs(1.0 - abs(ovr)))
                    plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))
                    plot_y[5].append(abs(em*dt)**3 / 6.0 * itime)

                # --- Time Evolution Equations ---
                # These are multiple steps of the D2 method (second order differencing scheme)
                # Fortran variables P, Q, R are used cyclically.
                # R(N) depends on Q(N-1), Q(N), Q(N+1) and P(N)
                # P(N) depends on R(N-1), R(N), R(N+1) and Q(N)
                # Q(N) depends on P(N-1), P(N), P(N+1) and R(N)

                r = np.zeros(NMAX, dtype=np.complex128)

                # R(1) = -2*IU*DT*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
                r[0] = -2 * iu * dt * (ha[0] * q[0] + hb[0] * q[1]) + p[0]

                # DO 70 N=2,NMAX-1
                # R(N)=-2*IU*DT*(HB(N-1)*Q(N-1)+HA(N)*Q(N)+HB(N)*Q(N+1))+P(N)
                # 70 CONTINUE
                for n in range(1, NMAX - 1): # Python index n corresponds to Fortran N+1
                    r[n] = -2 * iu * dt * (hb[n-1] * q[n-1] + ha[n] * q[n] + hb[n] * q[n+1]) + p[n]

                # R(NMAX)=-2*IU*DT*(HB(NMAX-1)*Q(NMAX-1)+HA(NMAX)*Q(NMAX))+P(NMAX)
                r[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * q[NMAX-2] + ha[NMAX-1] * q[NMAX-1]) + p[NMAX-1]

                # Update P using R and Q
                # P(1) = -2*IU*DT*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
                p[0] = -2 * iu * dt * (ha[0] * r[0] + hb[0] * r[1]) + q[0]

                # DO 71 N=2,NMAX-1
                # P(N)=-2*IU*DT*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))+Q(N)
                # 71 CONTINUE
                for n in range(1, NMAX - 1):
                    p[n] = -2 * iu * dt * (hb[n-1] * r[n-1] + ha[n] * r[n] + hb[n] * r[n+1]) + q[n]

                # P(NMAX)=-2*IU*DT*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))+Q(NMAX)
                p[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * r[NMAX-2] + ha[NMAX-1] * r[NMAX-1]) + q[NMAX-1]

                # Update Q using P and R
                # Q(1) = -2*IU*DT*(HA(1)*P(1)+HB(1)*P(2))+R(1)
                q[0] = -2 * iu * dt * (ha[0] * p[0] + hb[0] * p[1]) + r[0]

                # DO 72 N=2,NMAX-1
                # Q(N)=-2*IU*DT*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))+R(N)
                # 72 CONTINUE
                for n in range(1, NMAX - 1):
                    q[n] = -2 * iu * dt * (hb[n-1] * p[n-1] + ha[n] * p[n] + hb[n] * p[n+1]) + r[n]

                # Q(NMAX)=-2*IU*DT*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))+R(NMAX)
                q[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * p[NMAX-2] + ha[NMAX-1] * p[NMAX-1]) + r[NMAX-1]

            #
            # plot
            #
            plt.gca().clear()
            plt.title('ST2')
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
            plt.savefig("btst2.png")
                
    except Exception as e:
        print(f"An error occurred during file operations: {e}")

    print("Program finished. Results written to tst2.dat")

# Call the main program function
if __name__ == '__main__':
    nroom_program()
