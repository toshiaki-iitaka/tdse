import numpy as np
import math
import matplotlib.pyplot as plt

# --- Parameters ---
PI = 3.141592 # Using the exact PI from the Fortran code
NMAX = 100  # Maximum number of rooms/points
MTIME = 1500 # Total number of time steps for evolution loop

# --- Main Program: NROOM ---
def nroom_program():
    """
    Simulates N-room quantum mechanics using a 3-step time evolution method.
    This program models the time evolution of a quantum system with N states
    under a given Hamiltonian, specifically for a 1D chain of rooms.
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # --- I/O Files ---
    try:
        # Open the output file for writing. This will overwrite if it exists.
        f_p202_out = open('p202.dat', 'w')
    except Exception as e:
        print(f"Error opening p202.dat: {e}")
        return

    plotx = range(NMAX)
    plt.gca().clear()
    plt.title('boxn')
    plt.xlabel(r'n')
    plt.ylabel(r'prob')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlim([0,NMAX])
    plt.ylim([0,20/NMAX])

    # --- Initial Data ---
    dx = 1.0 # Spatial step size

    # --- Calculation of Hamiltonian ---
    a_coupling = 1.0 # This 'A' defines the coupling strength
    dt = (1.0 / a_coupling) * 1e-1 # Time step

    # HA (diagonal elements) and HB (off-diagonal elements) of the Hamiltonian
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX - 1, dtype=np.float64) # HB has NMAX-1 elements

    # DO 50 N=1,NMAX
    # HA(N) = 0
    # 50 CONTINUE
    # All diagonal elements are set to 0 (no on-site potential)
    for n_idx in range(NMAX):
        ha[n_idx] = 0.0

    # DO 51 N=1,NMAX-1
    # HB(N) = -A
    # 51 CONTINUE
    # All off-diagonal elements are set to -A (constant nearest-neighbor coupling)
    for n_idx in range(NMAX - 1):
        hb[n_idx] = -a_coupling

    # --- Set Initial Wave Function (Gaussian Wave Packet) ---
    # P, Q, R are the main state arrays (complex numbers)
    p_psi = np.zeros(NMAX, dtype=np.complex128)
    q_psi = np.zeros(NMAX, dtype=np.complex128)
    r_psi = np.zeros(NMAX, dtype=np.complex128) # R will be used as a state variable in the loop

    # Initial parameters for the Gaussian wave packet
    # X0, SG, PA are hardcoded in the Fortran, so we use these values directly.
    x0 = 0.3 * dx * NMAX
    sg = 0.05 * dx * NMAX
    pa_momentum = 0.1 * 2 * PI / dx

    # DO 60 N=1,NMAX
    # X=N*DX
    # P(N)= EXP(-0.5*(X-X0)**2/SG**2) * EXP(IU*PA*(X-X0))
    # 60 CONTINUE
    for n_idx in range(NMAX):
        x = (n_idx + 1) * dx # Fortran N*DX corresponds to (n_idx + 1) * dx
        # Calculate the initial wave function (Gaussian wave packet)
        p_psi[n_idx] = np.exp(-0.5 * (x - x0)**2 / sg**2) * np.exp(iu * pa_momentum * (x - x0))

    # --- Normalization of Wave Function ---
    pa_norm_factor = 0.0
    # DO 100 N=1,NMAX
    # PA = PA + CABS(P(N))**2
    # 100 CONTINUE
    # Calculate the sum of absolute squares to find the normalization factor
    for n_idx in range(NMAX):
        pa_norm_factor += np.abs(p_psi[n_idx])**2
    pa_norm_factor = np.sqrt(pa_norm_factor)

    # DO 110 N=1,NMAX
    # P(N)=P(N)/PA
    # 110 CONTINUE
    # Normalize the initial wave function
    for n_idx in range(NMAX):
        p_psi[n_idx] = p_psi[n_idx] / pa_norm_factor

    # --- Second Wave Function (Initialization for the 3-step method) ---
    # This block calculates the first time-step derivative and then uses it
    # to get the wave function at t=dt (stored in Q) using a predictor-corrector approach.

    # R(1)= -IU*DT*(HA(1)*P(1)+HB(1)*P(2))
    # Calculate the derivative at the first point
    r_psi[0] = -iu * dt * (ha[0] * p_psi[0] + hb[0] * p_psi[1])
    # DO 61 N=2,NMAX-1
    # R(N)= -IU*DT*( HB(N-1)*P(N-1) + HA(N)*P(N) + HB(N)*P(N+1) )
    # 61 CONTINUE
    # Calculate derivatives for interior points
    for n_idx in range(1, NMAX - 1): # Python index for N=2 to NMAX-1
        r_psi[n_idx] = -iu * dt * (hb[n_idx-1] * p_psi[n_idx-1] + ha[n_idx] * p_psi[n_idx] + hb[n_idx] * p_psi[n_idx+1])
    # R(NMAX)=-IU*DT*( HB(NMAX-1)*P(NMAX-1) + HA(NMAX)*P(NMAX) )
    # Calculate the derivative at the last point
    r_psi[NMAX-1] = -iu * dt * (hb[NMAX-2] * p_psi[NMAX-2] + ha[NMAX-1] * p_psi[NMAX-1])

    # DO 62 N=1,NMAX
    # Q(N)=P(N)+R(N)
    # 62 CONTINUE
    # Predictor step: Q_psi is approximately psi(t=dt)
    q_psi[:] = p_psi[:] + r_psi[:]

    # Corrector step for Q_psi
    # Q(1)= Q(1)-0.5*IU*DT*(HA(1)*R(1)+HB(1)*R(2))
    q_psi[0] = q_psi[0] - 0.5 * iu * dt * (ha[0] * r_psi[0] + hb[0] * r_psi[1])
    # DO 63 N=2,NMAX-1
    # Q(N)= Q(N)-0.5*IU*DT*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))
    # 63 CONTINUE
    for n_idx in range(1, NMAX - 1): # Python index for N=2 to NMAX-1
        q_psi[n_idx] = q_psi[n_idx] - 0.5 * iu * dt * (hb[n_idx-1] * r_psi[n_idx-1] + ha[n_idx] * r_psi[n_idx] + hb[n_idx] * r_psi[n_idx+1])
    # Q(NMAX)=Q(NMAX)-0.5*IU*DT*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))
    q_psi[NMAX-1] = q_psi[NMAX-1] - 0.5 * iu * dt * (hb[NMAX-2] * r_psi[NMAX-2] + ha[NMAX-1] * r_psi[NMAX-1])

    # At this point:
    # p_psi holds psi(t=0)
    # q_psi holds psi(t=dt)
    # r_psi holds the initial derivative (which will be overwritten in the main loop)

    # --- Time Evolution Loop ---
    # DO 1000 ITIME=1,MTIME,3
    # The loop iterates from 1 up to MTIME, incrementing by 3 in each step.
    for itime in range(1, MTIME + 1, 3):
        t_val = dt * itime # Current simulation time

        # --- Output Wave Function ---
        # IF(MOD(ITIME,30).EQ. 1) THEN
        # Output is performed when (ITIME % 30) is 1.
        if (itime % 30) == 1:
            # DO 150 N=1,NMAX
            # WRITE(*,*) N,CABS(Q(N))**2
            # WRITE(20,*) CABS(Q(N))**2
            # 150 CONTINUE
            plot_y = [[] for _ in range(1) ]
            for n_idx in range(NMAX):
                # Print to console (N and probability density)
                print(f"{n_idx + 1:4d} {np.abs(q_psi[n_idx])**2:15.7E}")
                # Write probability density to file
                f_p202_out.write(f"{np.abs(q_psi[n_idx])**2}\n")
                plot_y[0].append(abs(q_psi[n_idx])**2)
            # WRITE(20,*) ' ' - Write a blank line to the file after each output block
            f_p202_out.write("\n")
            plt.plot(plot_y[0],linestyle='solid')
            
        # --- Time Evolution (3-step method) ---
        # This section updates the wave functions P, Q, and R for the next time step.
        # The updates are "in-place" in the sense that the new R is calculated using the current Q and P,
        # then the new P is calculated using the new R and current Q, etc.

        # Step 1: Calculate new R based on current Q and P
        # R(1) = -2*IU*DT*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
        r_psi[0] = -2 * iu * dt * (ha[0] * q_psi[0] + hb[0] * q_psi[1]) + p_psi[0]
        # DO 70 N=2,NMAX-1
        # R(N)=-2*IU*DT*(HB(N-1)*Q(N-1)+HA(N)*Q(N)+HB(N)*Q(N+1))+P(N)
        # 70 CONTINUE
        for n_idx in range(1, NMAX - 1): # Python index for N=2 to NMAX-1
            r_psi[n_idx] = -2 * iu * dt * (hb[n_idx-1] * q_psi[n_idx-1] + ha[n_idx] * q_psi[n_idx] + hb[n_idx] * q_psi[n_idx+1]) + p_psi[n_idx]
        # R(NMAX)=-2*IU*DT*(HB(NMAX-1)*Q(NMAX-1)+HA(NMAX)*Q(NMAX))+P(NMAX)
        r_psi[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * q_psi[NMAX-2] + ha[NMAX-1] * q_psi[NMAX-1]) + p_psi[NMAX-1]

        # Step 2: Calculate new P based on the newly calculated R and current Q
        # P(1) = -2*IU*DT*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
        p_psi[0] = -2 * iu * dt * (ha[0] * r_psi[0] + hb[0] * r_psi[1]) + q_psi[0]
        # DO 71 N=2,NMAX-1
        # P(N)=-2*IU*DT*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))+Q(N)
        # 71 CONTINUE
        for n_idx in range(1, NMAX - 1): # Python index for N=2 to NMAX-1
            p_psi[n_idx] = -2 * iu * dt * (hb[n_idx-1] * r_psi[n_idx-1] + ha[n_idx] * r_psi[n_idx] + hb[n_idx] * r_psi[n_idx+1]) + q_psi[n_idx]
        # P(NMAX)=-2*IU*DT*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))+Q(NMAX)
        p_psi[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * r_psi[NMAX-2] + ha[NMAX-1] * r_psi[NMAX-1]) + q_psi[NMAX-1]

        # Step 3: Calculate new Q based on the newly calculated P and newly calculated R
        # Q(1) = -2*IU*DT*(HA(1)*P(1)+HB(1)*P(2))+R(1)
        q_psi[0] = -2 * iu * dt * (ha[0] * p_psi[0] + hb[0] * p_psi[1]) + r_psi[0]
        # DO 72 N=2,NMAX-1
        # Q(N)=-2*IU*DT*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))+R(N)
        # 72 CONTINUE
        for n_idx in range(1, NMAX - 1): # Python index for N=2 to NMAX-1
            q_psi[n_idx] = -2 * iu * dt * (hb[n_idx-1] * p_psi[n_idx-1] + ha[n_idx] * p_psi[n_idx] + hb[n_idx] * p_psi[n_idx+1]) + r_psi[n_idx]
        # Q(NMAX)=-2*IU*DT*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))+R(NMAX)
        q_psi[NMAX-1] = -2 * iu * dt * (hb[NMAX-2] * p_psi[NMAX-2] + ha[NMAX-1] * p_psi[NMAX-1]) + r_psi[NMAX-1]

    # Close the output file after the loop finishes
    plt.show()
    plt.savefig("boxn.png")
    
    f_p202_out.close()
    print("Program finished. Results written to p202.dat")

# Call the main program function
if __name__ == '__main__':
    nroom_program()
