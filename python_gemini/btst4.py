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
    # So, Fortran N becomes (n + 1) in Python for calculations involving N.
    for n in range(nmax_val):
        ph0[n] = np.sin((PI * m_val * (n + 1)) / (nmax_val + 1))

    # Normalization of wave function
    pa = 0.0
    for n in range(nmax_val):
        pa += abs(ph0[n])**2
    pa = np.sqrt(pa)

    for n in range(nmax_val):
        ph0[n] = ph0[n] / pa

    # Calculate eigenenergy (EM)
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

    # Initialize P, Q, R, S arrays
    # T and U are temporary arrays used within the time evolution loop
    p = np.zeros(NMAX, dtype=np.complex128)
    q = np.zeros(NMAX, dtype=np.complex128)
    r = np.zeros(NMAX, dtype=np.complex128)
    s = np.zeros(NMAX, dtype=np.complex128)

    # DO 90 N=1,NMAX
    # P(N)=PH0(N)
    # Q(N)=PH0(N)*EXP(-IU*EM*DT)
    # R(N)=PH0(N)*EXP(-IU*EM*2*DT)
    # S(N)=PH0(N)*EXP(-IU*EM*3*DT)
    # 90 CONTINUE
    for n in range(NMAX):
        p[n] = ph0[n]
        q[n] = ph0[n] * np.exp(-iu * em * dt)
        r[n] = ph0[n] * np.exp(-iu * em * 2 * dt)
        s[n] = ph0[n] * np.exp(-iu * em * 3 * dt)

    # --- Time Evolution Loop ---
    # Open the output file
    try:
        with open('tst4.dat', 'w') as f_out:
            plot_x = range(3, mtime + 1, 5)
            plot_y = [[] for _ in range(6) ]
            # DO 1000 ITIME=3,MTIME,5
            # Python range end is exclusive, so mtime + 1
            for itime in range(3, mtime + 1, 5):
                tt = dt * itime

                # --- Output Wave Function (every 10th step, starting from 3) ---
                # IF(MOD(ITIME,10).EQ. 3) THEN
                if (itime % 10) == 3:
                    anorm = 0.0
                    ovr = 0.0 + 0.0j
                    # DO 150 N=1,NMAX
                    # ANORM=ANORM+ABS(S(N))**2
                    # OVR=OVR+DCONJG(S(N))*PH0(N)*EXP(-IU*EM*TT)
                    # 150   CONTINUE
                    for n in range(NMAX):
                        anorm += abs(s[n])**2
                        ovr += np.conjugate(s[n]) * ph0[n] * np.exp(-iu * em * tt)

                    # Fortran console output:
                    # WRITE(*,100) ITIME,EMAX*TT,ANORM-1,
                    # &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
                    # &,ABS((EM*DT)**5*(7./90.)) * ITIME
                    # Format: I4, 10E15.7 (for 6 values, 10E15.7 is generous)
                    error_term = abs((em * dt)**5 * (7.0 / 90.0)) * itime
                    console_output = (
                        f"{itime:4d} {emax*tt:15.7E} {anorm-1:15.7E} "
                        f"{1.0 - abs(ovr):15.7E} {np.arctan2(ovr.imag, ovr.real):15.7E} "
                        f"{error_term:15.7E}"
                    )
                    print(console_output)

                    # Fortran file output:
                    # WRITE(20,100) ITIME,EMAX*TT,ABS(ANORM-1),
                    # &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
                    # &,ABS((EM*DT)**5*(7./90.)) * ITIME
                    file_output = (
                        f"{itime:4d} {emax*tt:15.7E} {abs(anorm-1):15.7E} "
                        f"{abs(1.0 - abs(ovr)):15.7E} {abs(np.arctan2(ovr.imag, ovr.real)):15.7E} "
                        f"{error_term:15.7E}"
                    )
                    f_out.write(file_output + "\n")

                    plot_y[0].append(itime)
                    plot_y[1].append(emax*tt)
                    plot_y[2].append(abs(anorm-1))
                    plot_y[3].append(abs(1.0 - abs(ovr)))
                    plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))
                    plot_y[5].append(abs(em*dt)**5 * 7/90.0 * itime)

                # --- Time Evolution Steps (Fifth-order method) ---
                # These involve cyclic updates of P, Q, R, S, T, U.
                # U is a temporary array for intermediate calculations.
                u = np.zeros(NMAX, dtype=np.complex128)
                t_new = np.zeros(NMAX, dtype=np.complex128) # Store T values for this iteration
                p_new = np.zeros(NMAX, dtype=np.complex128) # Store P values for this iteration
                q_new = np.zeros(NMAX, dtype=np.complex128) # Store Q values for this iteration
                r_new = np.zeros(NMAX, dtype=np.complex128) # Store R values for this iteration
                s_new = np.zeros(NMAX, dtype=np.complex128) # Store S values for this iteration


                # Step 1: Update T using P, Q, S and intermediate U
                # DO 80 N=1,NMAX
                # U(N) = (2*Q(N)-R(N)+2*S(N))/3.0
                # 80 CONTINUE
                for n in range(NMAX):
                    u[n] = (2 * q[n] - r[n] + 2 * s[n]) / 3.0

                # T(1)= -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2)) + P(1)
                t_new[0] = -4 * iu * dt * (ha[0] * u[0] + hb[0] * u[1]) + p[0]
                # DO 70 N=2,NMAX-1
                # T(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1)) + P(N)
                # 70 CONTINUE
                for n in range(1, NMAX - 1):
                    t_new[n] = -4 * iu * dt * (hb[n-1] * u[n-1] + ha[n] * u[n] + hb[n] * u[n+1]) + p[n]
                # T(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX)) + P(NMAX)
                t_new[NMAX-1] = -4 * iu * dt * (hb[NMAX-2] * u[NMAX-2] + ha[NMAX-1] * u[NMAX-1]) + p[NMAX-1]


                # Step 2: Update P using Q, R, T and intermediate U
                # DO 81 N=1,NMAX
                # U(N) = (2*R(N)-S(N)+2*T(N))/3.0
                # 81 CONTINUE
                for n in range(NMAX):
                    u[n] = (2 * r[n] - s[n] + 2 * t_new[n]) / 3.0 # Use t_new here

                # P(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2)) +Q(1)
                p_new[0] = -4 * iu * dt * (ha[0] * u[0] + hb[0] * u[1]) + q[0]
                # DO 71 N=2,NMAX-1
                # P(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+Q(N)
                # 71 CONTINUE
                for n in range(1, NMAX - 1):
                    p_new[n] = -4 * iu * dt * (hb[n-1] * u[n-1] + ha[n] * u[n] + hb[n] * u[n+1]) + q[n]
                # P(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX)) + Q(NMAX)
                p_new[NMAX-1] = -4 * iu * dt * (hb[NMAX-2] * u[NMAX-2] + ha[NMAX-1] * u[NMAX-1]) + q[NMAX-1]


                # Step 3: Update Q using R, S, P and intermediate U
                # DO 82 N=1,NMAX
                # U(N) = (2*S(N)-T(N)+2*P(N))/3.0
                # 82 CONTINUE
                for n in range(NMAX):
                    u[n] = (2 * s[n] - t_new[n] + 2 * p_new[n]) / 3.0 # Use t_new and p_new

                # Q(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2)) + R(1)
                q_new[0] = -4 * iu * dt * (ha[0] * u[0] + hb[0] * u[1]) + r[0]
                # DO 72 N=2,NMAX-1
                # Q(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+R(N)
                # 72 CONTINUE
                for n in range(1, NMAX - 1):
                    q_new[n] = -4 * iu * dt * (hb[n-1] * u[n-1] + ha[n] * u[n] + hb[n] * u[n+1]) + r[n]
                # Q(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX))+R(NMAX)
                q_new[NMAX-1] = -4 * iu * dt * (hb[NMAX-2] * u[NMAX-2] + ha[NMAX-1] * u[NMAX-1]) + r[NMAX-1]


                # Step 4: Update R using S, T, Q and intermediate U
                # DO 83 N=1,NMAX
                # U(N) = (2*T(N)-P(N)+2*Q(N))/3.0
                # 83 CONTINUE
                for n in range(NMAX):
                    u[n] = (2 * t_new[n] - p_new[n] + 2 * q_new[n]) / 3.0 # Use t_new, p_new, q_new

                # R(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2))+S(1)
                r_new[0] = -4 * iu * dt * (ha[0] * u[0] + hb[0] * u[1]) + s[0]
                # DO 73 N=2,NMAX-1
                # R(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+S(N)
                # 73 CONTINUE
                for n in range(1, NMAX - 1):
                    r_new[n] = -4 * iu * dt * (hb[n-1] * u[n-1] + ha[n] * u[n] + hb[n] * u[n+1]) + s[n]
                # R(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX))+S(NMAX)
                r_new[NMAX-1] = -4 * iu * dt * (hb[NMAX-2] * u[NMAX-2] + ha[NMAX-1] * u[NMAX-1]) + s[NMAX-1]


                # Step 5: Update S using P, Q, R and intermediate U
                # DO 84 N=1,NMAX
                # U(N) = (2*P(N)-Q(N)+2*R(N))/3.0
                # 84 CONTINUE
                for n in range(NMAX):
                    u[n] = (2 * p_new[n] - q_new[n] + 2 * r_new[n]) / 3.0 # Use p_new, q_new, r_new

                # S(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2))+T(1)
                s_new[0] = -4 * iu * dt * (ha[0] * u[0] + hb[0] * u[1]) + t_new[0] # Use t_new
                # DO 74 N=2,NMAX-1
                # S(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+T(N)
                # 74 CONTINUE
                for n in range(1, NMAX - 1):
                    s_new[n] = -4 * iu * dt * (hb[n-1] * u[n-1] + ha[n] * u[n] + hb[n] * u[n+1]) + t_new[n] # Use t_new
                # S(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX))+T(NMAX)
                s_new[NMAX-1] = -4 * iu * dt * (hb[NMAX-2] * u[NMAX-2] + ha[NMAX-1] * u[NMAX-1]) + t_new[NMAX-1] # Use t_new

                # Update the main wave function arrays for the next iteration
                p[:] = p_new[:]
                q[:] = q_new[:]
                r[:] = r_new[:]
                s[:] = s_new[:]

            #
            # plot
            #
            plt.gca().clear()
            plt.title('ST4')
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
            plt.savefig("btst4.png")

    except Exception as e:
        print(f"An error occurred during file operations: {e}")

    print("Program finished. Results written to tst4.dat")

# Call the main program function
if __name__ == '__main__':
    nroom_program()
