import numpy as np
import math
import matplotlib.pyplot as plt

# --- Parameters ---
PI = 3.14159265358979323 # Using full precision PI for consistency
NMAX = 100  # Maximum number of rooms/points
MTIME = 3000 # Total number of time steps for evolution loop

# --- Main Program: ONED ---
def oned_program():
    """
    Main program to simulate 1D quantum mechanics with fixed boundary conditions
    using a 3-step time evolution method.
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # --- I/O File ---
    try:
        f_p304_out = open('p304.dat', 'w')
    except Exception as e:
        print(f"Error opening p304.dat: {e}")
        return

    plotx = range(NMAX)
    plt.gca().clear()
    plt.title('one-d')
    plt.xlabel(r'n')
    plt.ylabel(r'prob')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlim([0,NMAX])
    plt.ylim([0,20/NMAX])

    # --- Initial Data ---
    dx = 1.0
    a_const = 0.5 / dx / dx
    alpha = 0.03

    # V0 is commented out for user input and hardcoded to 1
    v0 = 1.0
    nl = 50
    nr = 60

    # --- Calculation of Hamiltonian ---
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX - 1, dtype=np.float64) # HB is NMAX-1 elements

    # DO 50 N=1,NMAX
    # IF(NL.LE.N .AND. N.LE.NR) THEN
    #   HA(N) = V0
    # ELSE
    #   HA(N) = 0.0
    # ENDIF
    # 50 CONTINUE
    for n_idx in range(NMAX):
        # Fortran N is (n_idx + 1) in Python
        if nl <= (n_idx + 1) <= nr:
            ha[n_idx] = v0
        else:
            ha[n_idx] = 0.0

    # DO 51 N=1,NMAX-1
    # HB(N) = -A
    # 51 CONTINUE
    for n_idx in range(NMAX - 1):
        hb[n_idx] = -a_const

    # --- Setting Time Step ---
    vmax = 0.0
    vmin = 0.0
    # DO 10 N=1,NMAX
    # IF(HA(N).GT.VMAX) VMAX=HA(N)
    # IF(HA(N).LT.VMIN) VMIN=HA(N)
    # 10 CONTINUE
    vmax = np.max(ha)
    vmin = np.min(ha)

    emax = a_const + vmax
    emin = -a_const + vmin
    ec = max(abs(emax), abs(emin))
    dt = alpha / ec

    # --- Setting Initial Wave Function ---
    # P, Q, R are the main state arrays
    p_psi = np.zeros(NMAX, dtype=np.complex128)
    q_psi = np.zeros(NMAX, dtype=np.complex128)
    r_psi = np.zeros(NMAX, dtype=np.complex128) # R will be used as a state variable in the loop

    x0 = 0.4 * dx * NMAX
    sg = 0.05 * dx * NMAX
    pa_momentum = 0.1 * 2 * PI / dx

    # DO 60 N=1,NMAX
    # X=N*DX
    # P(N)= EXP(-0.5*(X-X0)**2/SG**2) * EXP(IU*PA*(X-X0))
    # 60 CONTINUE
    for n_idx in range(NMAX):
        x = (n_idx + 1) * dx # Fortran N*DX
        p_psi[n_idx] = np.exp(-0.5 * (x - x0)**2 / sg**2) * np.exp(iu * pa_momentum * (x - x0))

    # --- Normalization of Wave Function ---
    pa_norm_factor = 0.0
    # DO 100 N=1,NMAX
    # PA = PA + CABS(P(N))**2
    # 100 CONTINUE
    for n_idx in range(NMAX):
        pa_norm_factor += abs(p_psi[n_idx])**2
    pa_norm_factor = np.sqrt(pa_norm_factor)

    # DO 110 N=1,NMAX
    # P(N)=P(N)/PA
    # 110 CONTINUE
    for n_idx in range(NMAX):
        p_psi[n_idx] = p_psi[n_idx] / pa_norm_factor

    # --- Second OAWave Function (Initial step for 3-step method) ---
    # R(1)= -IU*DT*(HA(1)*P(1)+HB(1)*P(2))
    r_temp_deriv = np.zeros(NMAX, dtype=np.complex128) # Temporary for derivative calculation
    r_temp_deriv[0] = -iu * dt * (ha[0] * p_psi[0] + hb[0] * p_psi[1])
    # DO 61 N=2,NMAX-1
    # R(N)= -IU*DT*( HB(N-1)*P(N-1) + HA(N)*P(N) + HB(N)*P(N+1) )
    # 61 CONTINUE
    for n_idx in range(1, NMAX - 1): # Python index for N=2 to NMAX-1
        r_temp_deriv[n_idx] = -iu * dt * (hb[n_idx-1] * p_psi[n_idx-1] + ha[n_idx] * p_psi[n_idx] + hb[n_idx] * p_psi[n_idx+1])
    # R(NMAX)=-IU*DT*( HB(NMAX-1)*P(NMAX-1) + HA(NMAX)*P(NMAX) )
    r_temp_deriv[NMAX-1] = -iu * dt * (hb[NMAX-2] * p_psi[NMAX-2] + ha[NMAX-1] * p_psi[NMAX-1])

    # DO 62 N=1,NMAX
    # Q(N)=P(N)+R(N)
    # 62 CONTINUE
    q_psi[:] = p_psi[:] + r_temp_deriv[:] # Q_psi is psi(t=dt) approx.

    # Refine Q_psi (corrector step)
    # Q(1)= Q(1)-0.5*IU*DT*(HA(1)*R(1)+HB(1)*R(2))
    r_temp_deriv_2 = np.zeros(NMAX, dtype=np.complex128) # Another temporary for derivative
    r_temp_deriv_2[0] = -0.5 * iu * dt * (ha[0] * r_temp_deriv[0] + hb[0] * r_temp_deriv[1])
    # DO 63 N=2,NMAX-1
    # Q(N)= Q(N)-0.5*IU*DT*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))
    # 63 CONTINUE
    for n_idx in range(1, NMAX - 1): # Python index for N=2 to NMAX-1
        r_temp_deriv_2[n_idx] = -0.5 * iu * dt * (hb[n_idx-1] * r_temp_deriv[n_idx-1] + ha[n_idx] * r_temp_deriv[n_idx] + hb[n_idx] * r_temp_deriv[n_idx+1])
    # Q(NMAX)=Q(NMAX)-0.5*IU*DT*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))
    r_temp_deriv_2[NMAX-1] = -0.5 * iu * dt * (hb[NMAX-2] * r_temp_deriv[NMAX-2] + ha[NMAX-1] * r_temp_deriv[NMAX-1])

    q_psi[:] = q_psi[:] - r_temp_deriv_2[:] # Final Q_psi (psi(t=dt) refined)

    # At this point:
    # p_psi holds psi(t=0)
    # q_psi holds psi(t=dt)
    # r_psi is currently uninitialized as a state variable, it will be calculated in the first loop iteration.

    # --- Time Evolution Loop (3-step method) ---
    # DO 1000 ITIME=1,MTIME,3
    for itime in range(1, MTIME + 1, 3): # Loop from 1 to MTIME, step by 3
        t_val = dt * itime

        # --- Output Wave Function ---
        # IF(MOD(ITIME,10).EQ. 0) THEN
        if (itime % 100) == 0:
            # DO 150 N=1,NMAX
            # WRITE(20,*) CABS(Q(N))**2
            # 150 CONTINUE
            plot_y = [[] for _ in range(1) ]
            for n_idx in range(NMAX):
                f_p304_out.write(f"{abs(q_psi[n_idx])**2}\n") # Output Q (psi(t+dt))
                plot_y[0].append(abs(q_psi[n_idx])**2)
            plt.plot(plot_y[0],linestyle='solid')

        # --- Time Evolution (3-step method) - No In-Place Update ---
        # Capture current states for calculation
        p_curr = p_psi.copy()
        q_curr = q_psi.copy()
        r_curr = r_psi.copy() # This will be zero for ITIME=1, then the calculated R from previous step

        # Calculate the next states
        # R_next_calc: Corresponds to R(N) in Fortran's first block
        # R(N) = -2*IU*DT*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
        r_next_calc = np.zeros(NMAX, dtype=np.complex128)
        r_next_calc[0] = -2 * iu * dt * (ha[0] * q_curr[0] + hb[0] * q_curr[1]) + p_curr[0]
        for n_idx in range(1, NMAX - 1):
            r_next_calc[n_idx] = -2 * iu * dt * (hb[n_idx-1] * q_curr[n_idx-1] + ha[n_idx] * q_curr[n_idx] + hb[n_idx] * q_curr[n_idx+1]) + p_curr[n_idx]
        r_next_calc[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * q_curr[NMAX-2] + ha[NMAX-1] * q_curr[NMAX-1]) + p_curr[NMAX-1]

        # P_next_calc: Corresponds to P(N) in Fortran's second block
        # P(N) = -2*IU*DT*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
        # Uses r_next_calc (the newly calculated R) and q_curr (the old Q)
        p_next_calc = np.zeros(NMAX, dtype=np.complex128)
        p_next_calc[0] = -2 * iu * dt * (ha[0] * r_next_calc[0] + hb[0] * r_next_calc[1]) + q_curr[0]
        for n_idx in range(1, NMAX - 1):
            p_next_calc[n_idx] = -2 * iu * dt * (hb[n_idx-1] * r_next_calc[n_idx-1] + ha[n_idx] * r_next_calc[n_idx] + hb[n_idx] * r_next_calc[n_idx+1]) + q_curr[n_idx]
        p_next_calc[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * r_next_calc[NMAX-2] + ha[NMAX-1] * r_next_calc[NMAX-1]) + q_curr[NMAX-1]

        # Q_next_calc: Corresponds to Q(N) in Fortran's third block
        # Q(N) = -2*IU*DT*(HA(1)*P(1)+HB(1)*P(2))+R(1)
        # Uses p_next_calc (the newly calculated P) and r_next_calc (the newly calculated R)
        q_next_calc = np.zeros(NMAX, dtype=np.complex128)
        q_next_calc[0] = -2 * iu * dt * (ha[0] * p_next_calc[0] + hb[0] * p_next_calc[1]) + r_next_calc[0]
        for n_idx in range(1, NMAX - 1):
            q_next_calc[n_idx] = -2 * iu * dt * (hb[n_idx-1] * p_next_calc[n_idx-1] + ha[n_idx] * p_next_calc[n_idx] + hb[n_idx] * p_next_calc[n_idx+1]) + r_next_calc[n_idx]
        q_next_calc[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * p_next_calc[NMAX-2] + ha[NMAX-1] * p_next_calc[NMAX-1]) + r_next_calc[NMAX-1]

        # Update the main state variables for the next iteration
        # The states shift:
        # P_main (old psi(t)) becomes Q_curr (old psi(t+dt))
        # Q_main (old psi(t+dt)) becomes R_next_calc (new psi(t+2dt))
        # R_main (old psi(t+2dt)) becomes P_next_calc (new psi(t+3dt))
        p_psi[:] = q_curr[:]
        q_psi[:] = r_next_calc[:]
        r_psi[:] = p_next_calc[:]

    # Close the output file
    plt.show()
    plt.savefig("one-d.png")
    
    f_p304_out.close()
    print("Program finished. Results written to p304.dat")

# Call the main program function
if __name__ == '__main__':
    oned_program()
