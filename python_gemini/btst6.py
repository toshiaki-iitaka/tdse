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

    # Fortran code had V0(N) and part of HA(N) commented out.
    # HA(N) is effectively 0 based on the provided Fortran.
    # DO 50 N=1,NMAX
    # C        V0(N)=0
    # C        HA(N) = V0(N) + 1/DX2
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

    # Initialize P, Q, R, S, T, U arrays
    # V and W are temporary arrays used within the time evolution loop
    p = np.zeros(NMAX, dtype=np.complex128)
    q = np.zeros(NMAX, dtype=np.complex128)
    r = np.zeros(NMAX, dtype=np.complex128)
    s = np.zeros(NMAX, dtype=np.complex128)
    t = np.zeros(NMAX, dtype=np.complex128)
    u = np.zeros(NMAX, dtype=np.complex128)

    # DO 90 N=1,NMAX
    # P(N)=PH0(N)
    # Q(N)=PH0(N)*EXP(-IU*EM*DT)
    # R(N)=PH0(N)*EXP(-IU*EM*2*DT)
    # S(N)=PH0(N)*EXP(-IU*EM*3*DT)
    # T(N)=PH0(N)*EXP(-IU*EM*4*DT)
    # U(N)=PH0(N)*EXP(-IU*EM*5*DT)
    # 90 CONTINUE
    for n in range(NMAX):
        p[n] = ph0[n]
        q[n] = ph0[n] * np.exp(-iu * em * dt)
        r[n] = ph0[n] * np.exp(-iu * em * 2 * dt)
        s[n] = ph0[n] * np.exp(-iu * em * 3 * dt)
        t[n] = ph0[n] * np.exp(-iu * em * 4 * dt)
        u[n] = ph0[n] * np.exp(-iu * em * 5 * dt)

    # --- Time Evolution Loop ---
    # Open the output file
    try:
        with open('tst6.dat', 'w') as f_out:
            plot_x = range(5, mtime + 1, 7)
            plot_y = [[] for _ in range(6) ]
            # DO 1000 ITIME=5,MTIME,7
            # Python range end is exclusive, so mtime + 1
            for itime in range(5, mtime + 1, 7):
                tt = dt * itime

                # --- Output Wave Function (every 14th step, starting from 5) ---
                # IF(MOD(ITIME,14).EQ. 5) THEN
                if (itime % 14) == 5:
                    anorm = 0.0
                    ovr = 0.0 + 0.0j
                    # DO 150 N=1,NMAX
                    # ANORM=ANORM+ABS(U(N))**2
                    # OVR=OVR+DCONJG(U(N))*PH0(N)*EXP(-IU*EM*TT)
                    # 150   CONTINUE
                    for n in range(NMAX):
                        anorm += abs(u[n])**2
                        ovr += np.conjugate(u[n]) * ph0[n] * np.exp(-iu * em * tt)

                    # Fortran console output:
                    # WRITE(*,100) ITIME,EMAX*TT,ANORM-1,
                    # &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
                    # &,ABS((EM*DT)**7*(41./840.)) * ITIME
                    # Format: I4, 10E15.7 (for 6 values, 10E15.7 is generous)
                    error_term = abs((em * dt)**7 * (41.0 / 840.0)) * itime
                    console_output = (
                        f"{itime:4d} {emax*tt:15.7E} {anorm-1:15.7E} "
                        f"{1.0 - abs(ovr):15.7E} {np.arctan2(ovr.imag, ovr.real):15.7E} "
                        f"{error_term:15.7E}"
                    )
                    print(console_output)

                    # Fortran file output:
                    # WRITE(20,100) ITIME,EMAX*TT,ABS(ANORM-1),
                    # &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
                    # &,ABS((EM*DT)**7*(41./840.)) * ITIME
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
                    plot_y[5].append(abs(em*dt)**7 * 41.0/840.0 * itime)

                # --- Time Evolution Steps (Seventh-order method) ---
                # These involve cyclic updates of P, Q, R, S, T, U, V.
                # W is a temporary array for intermediate calculations.
                w = np.zeros(NMAX, dtype=np.complex128)
                v_new = np.zeros(NMAX, dtype=np.complex128)
                p_new = np.zeros(NMAX, dtype=np.complex128)
                q_new = np.zeros(NMAX, dtype=np.complex128)
                r_new = np.zeros(NMAX, dtype=np.complex128)
                s_new = np.zeros(NMAX, dtype=np.complex128)
                t_new = np.zeros(NMAX, dtype=np.complex128)
                u_new = np.zeros(NMAX, dtype=np.complex128)


                # Step 1: Update V using P, Q, R, S, T, U and intermediate W
                # DO 80 N=1,NMAX
                # W(N) = (11*Q(N)-14*R(N)+26*S(N)-14*T(N)+11*U(N))/20.0
                # 80 CONTINUE
                for n in range(NMAX):
                    w[n] = (11 * q[n] - 14 * r[n] + 26 * s[n] - 14 * t[n] + 11 * u[n]) / 20.0

                # V(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + P(1)
                v_new[0] = -6 * iu * dt * (ha[0] * w[0] + hb[0] * w[1]) + p[0]
                # DO 70 N=2,NMAX-1
                # V(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + P(N)
                # 70 CONTINUE
                for n in range(1, NMAX - 1):
                    v_new[n] = -6 * iu * dt * (hb[n-1] * w[n-1] + ha[n] * w[n] + hb[n] * w[n+1]) + p[n]
                # V(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + P(NMAX)
                v_new[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w[NMAX-2] + ha[NMAX-1] * w[NMAX-1]) + p[NMAX-1]


                # Step 2: Update P using Q, R, S, T, U, V and intermediate W
                # DO 81 N=1,NMAX
                # W(N) = (11*R(N)-14*S(N)+26*T(N)-14*U(N)+11*V(N))/20.0
                # 81 CONTINUE
                for n in range(NMAX):
                    w[n] = (11 * r[n] - 14 * s[n] + 26 * t[n] - 14 * u[n] + 11 * v_new[n]) / 20.0 # Use v_new

                # P(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + Q(1)
                p_new[0] = -6 * iu * dt * (ha[0] * w[0] + hb[0] * w[1]) + q[0]
                # DO 71 N=2,NMAX-1
                # P(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + Q(N)
                # 71 CONTINUE
                for n in range(1, NMAX - 1):
                    p_new[n] = -6 * iu * dt * (hb[n-1] * w[n-1] + ha[n] * w[n] + hb[n] * w[n+1]) + q[n]
                # P(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + Q(NMAX)
                p_new[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w[NMAX-2] + ha[NMAX-1] * w[NMAX-1]) + q[NMAX-1]


                # Step 3: Update Q using R, S, T, U, V, P and intermediate W
                # DO 82 N=1,NMAX
                # W(N) = (11*S(N)-14*T(N)+26*U(N)-14*V(N)+11*P(N))/20.0
                # 82 CONTINUE
                for n in range(NMAX):
                    w[n] = (11 * s[n] - 14 * t[n] + 26 * u[n] - 14 * v_new[n] + 11 * p_new[n]) / 20.0 # Use v_new, p_new

                # Q(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + R(1)
                q_new[0] = -6 * iu * dt * (ha[0] * w[0] + hb[0] * w[1]) + r[0]
                # DO 72 N=2,NMAX-1
                # Q(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + R(N)
                # 72 CONTINUE
                for n in range(1, NMAX - 1):
                    q_new[n] = -6 * iu * dt * (hb[n-1] * w[n-1] + ha[n] * w[n] + hb[n] * w[n+1]) + r[n]
                # Q(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + R(NMAX)
                q_new[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w[NMAX-2] + ha[NMAX-1] * w[NMAX-1]) + r[NMAX-1]


                # Step 4: Update R using S, T, U, V, P, Q and intermediate W
                # DO 83 N=1,NMAX
                # W(N) = (11*T(N)-14*U(N)+26*V(N)-14*P(N)+11*Q(N))/20.0
                # 83 CONTINUE
                for n in range(NMAX):
                    w[n] = (11 * t[n] - 14 * u[n] + 26 * v_new[n] - 14 * p_new[n] + 11 * q_new[n]) / 20.0 # Use v_new, p_new, q_new

                # R(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + S(1)
                r_new[0] = -6 * iu * dt * (ha[0] * w[0] + hb[0] * w[1]) + s[0]
                # DO 73 N=2,NMAX-1
                # R(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + S(N)
                # 73 CONTINUE
                for n in range(1, NMAX - 1):
                    r_new[n] = -6 * iu * dt * (hb[n-1] * w[n-1] + ha[n] * w[n] + hb[n] * w[n+1]) + s[n]
                # R(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + S(NMAX)
                r_new[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w[NMAX-2] + ha[NMAX-1] * w[NMAX-1]) + s[NMAX-1]


                # Step 5: Update S using T, U, V, P, Q, R and intermediate W
                # DO 84 N=1,NMAX
                # W(N) = (11*U(N)-14*V(N)+26*P(N)-14*Q(N)+11*R(N))/20.0
                # 84 CONTINUE
                for n in range(NMAX):
                    w[n] = (11 * u[n] - 14 * v_new[n] + 26 * p_new[n] - 14 * q_new[n] + 11 * r_new[n]) / 20.0 # Use v_new, p_new, q_new, r_new

                # S(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + T(1)
                s_new[0] = -6 * iu * dt * (ha[0] * w[0] + hb[0] * w[1]) + t[0]
                # DO 74 N=2,NMAX-1
                # S(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + T(N)
                # 74 CONTINUE
                for n in range(1, NMAX - 1):
                    s_new[n] = -6 * iu * dt * (hb[n-1] * w[n-1] + ha[n] * w[n] + hb[n] * w[n+1]) + t[n]
                # S(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + T(NMAX)
                s_new[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w[NMAX-2] + ha[NMAX-1] * w[NMAX-1]) + t[NMAX-1]


                # Step 6: Update T using U, V, P, Q, R, S and intermediate W
                # DO 85 N=1,NMAX
                # W(N) = (11*V(N)-14*P(N)+26*Q(N)-14*R(N)+11*S(N))/20.0
                # 85 CONTINUE
                for n in range(NMAX):
                    w[n] = (11 * v_new[n] - 14 * p_new[n] + 26 * q_new[n] - 14 * r_new[n] + 11 * s_new[n]) / 20.0 # Use v_new, p_new, q_new, r_new, s_new

                # T(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + U(1)
                t_new[0] = -6 * iu * dt * (ha[0] * w[0] + hb[0] * w[1]) + u[0]
                # DO 75 N=2,NMAX-1
                # T(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + U(N)
                # 75 CONTINUE
                for n in range(1, NMAX - 1):
                    t_new[n] = -6 * iu * dt * (hb[n-1] * w[n-1] + ha[n] * w[n] + hb[n] * w[n+1]) + u[n]
                # T(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + U(NMAX)
                t_new[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w[NMAX-2] + ha[NMAX-1] * w[NMAX-1]) + u[NMAX-1]


                # Step 7: Update U using V, P, Q, R, S, T and intermediate W
                # DO 86 N=1,NMAX
                # W(N) = (11*P(N)-14*Q(N)+26*R(N)-14*S(N)+11*T(N))/20.0
                # 86 CONTINUE
                for n in range(NMAX):
                    w[n] = (11 * p_new[n] - 14 * q_new[n] + 26 * r_new[n] - 14 * s_new[n] + 11 * t_new[n]) / 20.0 # Use p_new, q_new, r_new, s_new, t_new

                # U(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + V(1)
                u_new[0] = -6 * iu * dt * (ha[0] * w[0] + hb[0] * w[1]) + v_new[0] # Use v_new
                # DO 76 N=2,NMAX-1
                # U(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + V(N)
                # 76 CONTINUE
                for n in range(1, NMAX - 1):
                    u_new[n] = -6 * iu * dt * (hb[n-1] * w[n-1] + ha[n] * w[n] + hb[n] * w[n+1]) + v_new[n] # Use v_new
                # U(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + V(NMAX)
                u_new[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w[NMAX-2] + ha[NMAX-1] * w[NMAX-1]) + v_new[NMAX-1] # Use v_new

                # Update the main wave function arrays for the next iteration
                p[:] = p_new[:]
                q[:] = q_new[:]
                r[:] = r_new[:]
                s[:] = s_new[:]
                t[:] = t_new[:]
                u[:] = u_new[:]
                # V is the newest, so it becomes the 'oldest' for the next cycle
                # The Fortran code does not explicitly update V (it's updated into V_new and then used)
                # The cycle is P,Q,R,S,T,U,V -> P(new), Q(new), R(new), S(new), T(new), U(new), V(new)
                # And the next iteration uses these new values.
                # So we need to assign v_new to v.
                #v[:] = v_new[:]

            #
            # plot
            #
            plt.gca().clear()
            plt.title('ST6')
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
            plt.savefig("btst6.png")
                
    except Exception as e:
        print(f"An error occurred during file operations: {e}")

    print("Program finished. Results written to tst6.dat")

# Call the main program function
if __name__ == '__main__':
    nroom_program()
