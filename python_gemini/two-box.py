import numpy as np
import math
import matplotlib.pyplot as plt

# --- Parameters ---
PI = 3.14159265358979 # Using the exact PI from the Fortran code
MTIME = 10000 # Total number of time steps for evolution loop

# --- Main Program: TWOROOM ---
def tworoom_program():
    """
    Simulates 2-room quantum mechanics using a 3-step time evolution method.
    This program models the time evolution of a quantum system with two states
    (representing two rooms or sites) under a given Hamiltonian.
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # --- I/O File ---
    try:
        # Open the output file for writing. This will overwrite if it exists.
        f_p1061_out = open('p1061.dat', 'w')
    except Exception as e:
        print(f"Error opening p1061.dat: {e}")
        return

    # --- Initial Data Input ---
    print('INPUT PROBABILITY OF RIGHT AND LEFT ROOM')
    try:
        # Read initial probabilities for the two rooms
        pr_initial = float(input("Probability of Right Room (PR): "))
        pl_initial = float(input("Probability of Left Room (PL): "))
    except ValueError:
        print("Invalid input. Please enter numerical values for probabilities.")
        f_p1061_out.close()
        return

    # Initialize wave function array P (complex numbers)
    # P(1) and P(2) in Fortran correspond to P[0] and P[1] in Python
    p_psi = np.zeros(2, dtype=np.complex128)
    
    # Initialize P based on square root of probabilities
    # DSQRT in Fortran is np.sqrt in Python
    p_psi[0] = np.sqrt(pr_initial)
    p_psi[1] = np.sqrt(pl_initial)

    # Initialize Q and R arrays for the 3-step time evolution
    q_psi = np.zeros(2, dtype=np.complex128)
    r_psi = np.zeros(2, dtype=np.complex128)

    # Initialize Hamiltonian arrays HA (diagonal elements) and HB (off-diagonal)
    # HA(1) and HA(2) in Fortran correspond to HA[0] and HA[1] in Python
    ha = np.zeros(2, dtype=np.float64)
    # HB(1) in Fortran corresponds to HB[0] in Python
    hb = np.zeros(1, dtype=np.float64) 

    print('INPUT ENERGY OF RIGHT AND LEFT ROOM')
    try:
        # Read diagonal Hamiltonian elements (energies of each room)
        ha[0] = float(input("Energy of Right Room (HA(1)): "))
        ha[1] = float(input("Energy of Left Room (HA(2)): "))
    except ValueError:
        print("Invalid input. Please enter numerical values for energies.")
        f_p1061_out.close()
        return

    print('INPUT TRANSFER ENERGY A')
    try:
        # Read off-diagonal Hamiltonian element (coupling between rooms)
        hb[0] = float(input("Transfer Energy (HB(1)): "))
    except ValueError:
        print("Invalid input. Please enter a numerical value for transfer energy.")
        f_p1061_out.close()
        return

    # Define time step
    dt = 1e-3

    # --- Normalization of Initial Wave Function ---
    # PA = SQRT( CABS(P(1))**2 + CABS(P(2))**2)
    # CABS in Fortran is np.abs() in Python for complex numbers
    pa_norm_factor = np.sqrt(np.abs(p_psi[0])**2 + np.abs(p_psi[1])**2)
    
    # Normalize P array
    p_psi[0] = p_psi[0] / pa_norm_factor
    p_psi[1] = p_psi[1] / pa_norm_factor

    # --- Second Wave Function (Initialization for the 3-step method) ---
    # This block performs a predictor-corrector like step to get the
    # wave function at t=dt (stored in Q) and a temporary R derivative.

    # R(1)= -IU*DT*( HA(1)*P(1) + HB(1)*P(2))
    r_psi[0] = -iu * dt * (ha[0] * p_psi[0] + hb[0] * p_psi[1])
    # R(2)= -IU*DT*( HB(1)*P(1) + HA(2)*P(2) )
    r_psi[1] = -iu * dt * (hb[0] * p_psi[0] + ha[1] * p_psi[1])
    
    # Q(1)= P(1)+R(1)
    q_psi[0] = p_psi[0] + r_psi[0]
    # Q(2)= P(2)+R(2)
    q_psi[1] = p_psi[1] + r_psi[1]

    # Corrector step for Q
    # Q(1)= Q(1)-0.5*IU*DT*(HA(1)*R(1)+HB(1)*R(2))
    q_psi[0] = q_psi[0] - 0.5 * iu * dt * (ha[0] * r_psi[0] + hb[0] * r_psi[1])
    # Q(2)= Q(2)-0.5*IU*DT*(HB(1)*R(1)+HA(2)*R(2))
    q_psi[1] = q_psi[1] - 0.5 * iu * dt * (hb[0] * r_psi[0] + ha[1] * r_psi[1])

    # At this point:
    # p_psi holds psi(t=0)
    # q_psi holds psi(t=dt)
    # r_psi holds the initial derivative (which will be overwritten in the loop)

    # --- Time Evolution Loop ---
    plot_y = [[] for _ in range(3) ]
    # DO 1000 ITIME=1,MTIME,3
    for itime in range(1, MTIME + 1, 3): # Loop from 1 to MTIME, stepping by 3
        t_val = dt * itime # Current time

        # --- Output ---
        # PR=CABS(Q(1))**2
        # PL=CABS(Q(2))**2
        pr_output = np.abs(q_psi[0])**2 # Probability of right room
        pl_output = np.abs(q_psi[1])**2 # Probability of left room (not explicitly printed in Fortran's format)

        # WRITE(*,100) T, Q(1), PR, -PR
        # WRITE(20,100) T, Q(1), PR, -PR
        # Fortran FORMAT(1H ,10E15.7 ) suggests 10 floating point numbers,
        # but only 4 are provided in the WRITE statement.
        # We'll match the values provided: T, Q(1) real, Q(1) imag, PR, -PR
        # The Fortran code has a typo in the FORMAT, it should likely be 4E15.7 or similar.
        # Assuming 10E15.7 is a general format and it prints 4 values.
        
        # Format for console output
        console_output = (
            f"{t_val:15.7E} {q_psi[0].real:15.7E} {q_psi[0].imag:15.7E} "
            f"{pr_output:15.7E} {pl_output:15.7E}" # error of the fortran program is corrected
        )
        # Format for file output (same as console)
        file_output = console_output

        print(console_output)
        f_p1061_out.write(file_output + "\n")

        plot_y[0].append(t_val)
        plot_y[1].append(pr_output)
        plot_y[2].append(pl_output)

        # --- Time Evolution (3-step method) ---
        # These calculations update P, Q, and R in sequence, using the
        # most recently updated values within the same iteration.
        # This is essentially an in-place update for the three state variables.

        # Step 1: Calculate new R based on current Q and P
        # R(1) = -2*IU*DT*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
        r_psi[0] = -2 * iu * dt * (ha[0] * q_psi[0] + hb[0] * q_psi[1]) + p_psi[0]
        # R(2) = -2*IU*DT*(HB(1)*Q(1)+HA(2)*Q(2))+P(2)
        r_psi[1] = -2 * iu * dt * (hb[0] * q_psi[0] + ha[1] * q_psi[1]) + p_psi[1]

        # Step 2: Calculate new P based on current R (just calculated) and Q
        # P(1) = -2*IU*DT*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
        p_psi[0] = -2 * iu * dt * (ha[0] * r_psi[0] + hb[0] * r_psi[1]) + q_psi[0]
        # P(2) = -2*IU*DT*(HB(1)*R(1)+HA(2)*R(2))+Q(2)
        p_psi[1] = -2 * iu * dt * (hb[0] * r_psi[0] + ha[1] * r_psi[1]) + q_psi[1]

        # Step 3: Calculate new Q based on current P (just calculated) and R (just calculated)
        # Q(1) = -2*IU*DT*(HA(1)*P(1)+HB(1)*P(2))+R(1)
        q_psi[0] = -2 * iu * dt * (ha[0] * p_psi[0] + hb[0] * p_psi[1]) + r_psi[0]
        # Q(2) = -2*IU*DT*(HB(1)*P(1)+HA(2)*P(2))+R(2)
        q_psi[1] = -2 * iu * dt * (hb[0] * p_psi[0] + ha[1] * p_psi[1]) + r_psi[1]

    # Close the output file
    f_p1061_out.close()

    print("Program finished. Results written to p1061.dat")

    #
    # plot
    #
    plt.gca().clear()
    plt.title('two-box')
    plt.xlabel(r't')
    plt.ylabel(r'prob')
    plt.xscale('linear')
    plt.yscale('linear')
    #plt.xlim([0.1,100])
    plt.ylim([0,1])
    plt.plot(plot_y[0],plot_y[1], linestyle='solid')
    plt.plot(plot_y[0],plot_y[2], linestyle='dashed')
    #plt.plot(plot_y[1],plot_y[5], marker='None', linestyle='dashed')
    plt.show()
    plt.savefig("two-box.png")

# Call the main program function
if __name__ == '__main__':
    tworoom_program()
