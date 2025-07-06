import numpy as np
import math
from scipy.special import jv # For Bessel function of the first kind
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

# --- Subroutine CH: Implements Chebyshev propagation ---
def ch_subroutine(ha_norm, hb_norm, egrid, emin, ph0, t, nmax_val):
    """
    Performs time evolution using the Chebyshev propagation method.

    Args:
        ha_norm (np.ndarray): Normalized diagonal Hamiltonian elements.
        hb_norm (np.ndarray): Normalized off-diagonal Hamiltonian elements.
        egrid (float): Energy grid span (EMAX - EMIN).
        emin (float): Minimum energy.
        ph0 (np.ndarray): Initial wave function (complex array).
        t (float): Total time for propagation.
        nmax_val (int): Maximum number of points (NMAX).

    Returns:
        np.ndarray: The propagated wave function (PH).
    """
    iu = 0.0 + 1.0j
    eps = 1e-16 # Convergence criterion for Bessel function series
    imax = 1000 # Maximum iterations for Chebyshev series (from Fortran PARAMETER)

    # Initialize P and Q for the Chebyshev series
    p = np.zeros(nmax_val, dtype=np.complex128)
    q = np.zeros(nmax_val, dtype=np.complex128)
    r = np.zeros(nmax_val, dtype=np.complex128) # Temporary array for calculations
    ph = np.zeros(nmax_val, dtype=np.complex128) # Resulting propagated wave function

    # P(N)=PH0(N)
    for n_idx in range(nmax_val):
        p[n_idx] = ph0[n_idx]

    # Q(1) = -IU*(HA(1)*P(1)+HB(1)*P(2))
    q[0] = -iu * (ha_norm[0] * p[0] + hb_norm[0] * p[1])
    # DO 30 N=2,NMAX-1
    # Q(N) = -IU*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))
    # 30 CONTINUE
    for n_idx in range(1, nmax_val - 1):
        q[n_idx] = -iu * (hb_norm[n_idx-1] * p[n_idx-1] + ha_norm[n_idx] * p[n_idx] + hb_norm[n_idx] * p[n_idx+1])
    # Q(NMAX)= -IU*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))
    q[nmax_val-1] = -iu * (hb_norm[nmax_val-2] * p[nmax_val-2] + ha_norm[nmax_val-1] * p[nmax_val-1])

    # Initial terms of the Chebyshev expansion
    # AN=DJBES(0,EGRID*T/2)
    an_0 = jv(0, egrid * t / 2.0)
    # DO 40 N=1,NMAX
    # PH(N)=AN*P(N)
    # 40 CONTINUE
    for n_idx in range(nmax_val):
        ph[n_idx] = an_0 * p[n_idx]

    # AN=2*DJBES(1,EGRID*T/2)
    an_1 = 2 * jv(1, egrid * t / 2.0)
    # DO 45 N=1,NMAX
    # PH(N)=PH(N)+AN*Q(N)
    # 45 CONTINUE
    for n_idx in range(nmax_val):
        ph[n_idx] = ph[n_idx] + an_1 * q[n_idx]

    # Main Chebyshev series loop
    # DO 1000 I=2, IMAX, 3
    # The loop structure is a bit tricky due to Fortran's GOTO.
    # We will break the loop if convergence is met.
    for i in range(2, imax + 1, 3):
        # AN=2*DJBES(I,EGRID*T/2)
        an_i = 2 * jv(i, egrid * t / 2.0)
        if abs(an_i) < eps:
            print(f"Convergence reached at I={i}, AN={an_i:15.7E}")
            break # Exit loop if Bessel function term is negligible

        # R(1) = -2*IU*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
        r[0] = -2 * iu * (ha_norm[0] * q[0] + hb_norm[0] * q[1]) + p[0]
        # PH(1)=PH(1) +AN*R(1)
        ph[0] = ph[0] + an_i * r[0]
        # DO 70 N=2,NMAX-1
        # R(N)=-2*IU*(HB(N-1)*Q(N-1)+HA(N)*Q(N)+HB(N)*Q(N+1))+P(N)
        # PH(N)=PH(N)+AN*R(N)
        # 70 CONTINUE
        for n_idx in range(1, nmax_val - 1):
            r[n_idx] = -2 * iu * (hb_norm[n_idx-1] * q[n_idx-1] + ha_norm[n_idx] * q[n_idx] + hb_norm[n_idx] * q[n_idx+1]) + p[n_idx]
            ph[n_idx] = ph[n_idx] + an_i * r[n_idx]
        # R(NMAX)=-2*IU*(HB(NMAX-1)*Q(NMAX-1)+HA(NMAX)*Q(NMAX))+P(NMAX)
        r[nmax_val-1] = -2 * iu * (hb_norm[nmax_val-2] * q[nmax_val-2] + ha_norm[nmax_val-1] * q[nmax_val-1]) + p[nmax_val-1]
        # PH(NMAX)=PH(NMAX)+AN*R(NMAX)
        ph[nmax_val-1] = ph[nmax_val-1] + an_i * r[nmax_val-1]

        # AN=2*DJBES(I+1,EGRID*T/2)
        an_i_plus_1 = 2 * jv(i + 1, egrid * t / 2.0)
        if abs(an_i_plus_1) < eps:
            print(f"Convergence reached at I+1={i+1}, AN={an_i_plus_1:15.7E}")
            break

        # P(1) = -2*IU*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
        p_temp = -2 * iu * (ha_norm[0] * r[0] + hb_norm[0] * r[1]) + q[0]
        # PH(1)=PH(1) +AN*P(1)
        ph[0] = ph[0] + an_i_plus_1 * p_temp
        # DO 71 N=2,NMAX-1
        # P(N)=-2*IU*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))+Q(N)
        # PH(N)=PH(N)+AN*P(N)
        # 71 CONTINUE
        for n_idx in range(1, nmax_val - 1):
            p_temp_n = -2 * iu * (hb_norm[n_idx-1] * r[n_idx-1] + ha_norm[n_idx] * r[n_idx] + hb_norm[n_idx] * r[n_idx+1]) + q[n_idx]
            ph[n_idx] = ph[n_idx] + an_i_plus_1 * p_temp_n
            p[n_idx] = p_temp_n # Update p for next iteration
        # P(NMAX)=-2*IU*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))+Q(NMAX)
        p_temp_nmax = -2 * iu * (hb_norm[nmax_val-2] * r[nmax_val-2] + ha_norm[nmax_val-1] * r[nmax_val-1]) + q[nmax_val-1]
        ph[nmax_val-1] = ph[nmax_val-1] + an_i_plus_1 * p_temp_nmax
        p[nmax_val-1] = p_temp_nmax # Update p for next iteration
        p[0] = p_temp # Update p[0]

        # AN=2*DJBES(I+2,EGRID*T/2)
        an_i_plus_2 = 2 * jv(i + 2, egrid * t / 2.0)
        if abs(an_i_plus_2) < eps:
            print(f"Convergence reached at I+2={i+2}, AN={an_i_plus_2:15.7E}")
            break

        # Q(1) = -2*IU*(HA(1)*P(1)+HB(1)*P(2))+R(1)
        q_temp = -2 * iu * (ha_norm[0] * p[0] + hb_norm[0] * p[1]) + r[0]
        # PH(1)=PH(1) +AN*Q(1)
        ph[0] = ph[0] + an_i_plus_2 * q_temp
        # DO 72 N=2,NMAX-1
        # Q(N)=-2*IU*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))+R(N)
        # PH(N)=PH(N)+AN*Q(N)
        # 72 CONTINUE
        for n_idx in range(1, nmax_val - 1):
            q_temp_n = -2 * iu * (hb_norm[n_idx-1] * p[n_idx-1] + ha_norm[n_idx] * p[n_idx] + hb_norm[n_idx] * p[n_idx+1]) + r[n_idx]
            ph[n_idx] = ph[n_idx] + an_i_plus_2 * q_temp_n
            q[n_idx] = q_temp_n # Update q for next iteration
        # Q(NMAX)=-2*IU*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))+R(NMAX)
        q_temp_nmax = -2 * iu * (hb_norm[nmax_val-2] * p[nmax_val-2] + ha_norm[nmax_val-1] * p[nmax_val-1]) + r[nmax_val-1]
        ph[nmax_val-1] = ph[nmax_val-1] + an_i_plus_2 * q_temp_nmax
        q[nmax_val-1] = q_temp_nmax # Update q for next iteration
        q[0] = q_temp # Update q[0]

        print(f"I={i}, AN={an_i:15.7E}") # Fortran also prints this
        # If the loop finishes without breaking, it means it didn't converge within IMAX steps.
    else:
        print('CHEBYSHEV NOT CONVERGE')

    # Final transformation
    # DO 80 N=1,NMAX
    # PH(N)=PH(N)*EXP(-IU*(EGRID/2+EMIN)*T)
    # 80 CONTINUE
    for n_idx in range(nmax_val):
        ph[n_idx] = ph[n_idx] * np.exp(-iu * (egrid / 2.0 + emin) * t)

    return ph

# --- Main Program: TEST ---
def test_program():
    """
    Main program to simulate 1D quantum mechanics using the Chebyshev scheme.
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # Get input from user
    print('M')
    try:
        m = int(input())
    except ValueError:
        print("Invalid input. Please enter an integer for M.")
        return

    # --- Initial Data Calculations ---
    dx = 1.0
    dx2 = dx * dx
    emax = 1.0 / dx2 # EMAX from Fortran's initial calculation

    # --- Calculation of Hamiltonian ---
    # HA and HB arrays
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX - 1, dtype=np.float64)
    v0 = np.zeros(NMAX, dtype=np.float64) # V0 from Fortran, initialized to 0

    # Fortran code had V(N) commented out and HA(N) set to 0.
    # DO 50 N=1,NMAX
    # C        V(N)=0
    # C        HA(N) = V(N) + 1/DX2
    # HA(N)=0
    # 50 CONTINUE
    # -> ha remains all zeros as initialized

    # DO 51 N=1,NMAX-1
    #     HB(N) = -0.5/DX2
    # 51 CONTINUE
    for n_idx in range(NMAX - 1):
        hb[n_idx] = -0.5 / dx2

    # --- Normalization of Hamiltonian ---
    # Note: The Fortran code calculates VMAX/VMIN based on V(N),
    # but V(N) is explicitly set to 0 in the Fortran code.
    # This implies VMAX=0, VMIN=0.
    # The EMAX/EMIN are then hardcoded to 1/DX2 and -1/DX2.
    # I will follow the Fortran code's explicit assignments.
    emax_ham = 1.0 / dx2
    emin_ham = -1.0 / dx2
    egrid = emax_ham - emin_ham

    han = np.zeros(NMAX, dtype=np.float64)
    hbn = np.zeros(NMAX - 1, dtype=np.float64)

    # DO 20 N=1,NMAX
    # HAN(N)=(HA(N)-EGRID/2-EMIN)*2/EGRID
    # IF(N.LT.NMAX) HBN(N)=HB(N)*2/EGRID
    # 20 CONTINUE
    for n_idx in range(NMAX):
        han[n_idx] = (ha[n_idx] - egrid / 2.0 - emin_ham) * 2.0 / egrid
        if n_idx < NMAX - 1:
            hbn[n_idx] = hb[n_idx] * 2.0 / egrid

    print(f'EGRID={egrid:15.7E}')
    print('TIME,M') # This line is from Fortran, but M is an integer, so it's a bit odd.

    # --- Time Evolution Loop ---
    # Open the output file
    try:
        with open('tstch.dat', 'w') as f_out:
            plot_x = range(1, 41)
            plot_y = [[] for _ in range(5) ]
            # DO 1000 ITIME=1,40
            for itime in range(1, 41): # Loop from 1 to 40 inclusive
                t = 10**(itime / 10.0)
                print(f'T={t:15.7E}')
                print(f'ESTIMATED STEP={egrid * t / 2.0:15.7E}')

                # Call EIGEN subroutine
                ph0, em = eigen_subroutine(m, dx, NMAX)
                print(f'PERIOD={2 * PI / em:15.7E}')

                # Call CH subroutine for time evolution
                ph = ch_subroutine(han, hbn, egrid, emin_ham, ph0, t, NMAX)

                # --- Output Wave Function ---
                anorm = 0.0
                ovr = 0.0 + 0.0j
                # DO 150 N=1,NMAX
                # ANORM=ANORM+ABS(PH(N))**2
                # OVR=OVR+DCONJG(PH(N))*PH0(N)*EXP(-IU*EM*T)
                # 150   CONTINUE
                for n_idx in range(NMAX):
                    anorm += abs(ph[n_idx])**2
                    ovr += np.conjugate(ph[n_idx]) * ph0[n_idx] * np.exp(-iu * em * t)

                # Fortran console output:
                # WRITE(*,100) ITIME,EMAX*T,ANORM-1,
                # &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
                # Format: I4, 30E15.7 (for 5 values, 30E15.7 is generous)
                console_output = (
                    f"{itime:4d} {emax*t:15.7E} {anorm-1:15.7E} "
                    f"{1.0 - abs(ovr):15.7E} {np.arctan2(ovr.imag, ovr.real):15.7E}"
                )
                print(console_output)

                # Fortran file output:
                # WRITE(20,100) ITIME,EMAX*T,ABS(ANORM-1),
                # &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
                file_output = (
                    f"{itime:4d} {emax*t:15.7E} {abs(anorm-1):15.7E} "
                    f"{abs(1.0 - abs(ovr)):15.7E} {abs(np.arctan2(ovr.imag, ovr.real)):15.7E}"
                )
                f_out.write(file_output + "\n")
                    
                plot_y[0].append(itime)
                plot_y[1].append(emax*t)
                plot_y[2].append(abs(anorm-1))
                plot_y[3].append(abs(1.0 - abs(ovr)))
                plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))

            #
            # plot
            #
            plt.gca().clear()
            plt.title('Chevshev')
            plt.xlabel(r'Emax*t')
            plt.ylabel(r'error')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlim([1,100])
            plt.ylim([1e-20,1])
            plt.plot(plot_y[1],plot_y[2], marker='o',    linestyle='None',  label='ε_norm')
            plt.plot(plot_y[1],plot_y[4], marker='o',    linestyle='None',  label='ε_phase')
            #plt.plot(plot_y[1],plot_y[5], marker='None', linestyle='dashed',label='ε_theory')
            plt.legend()
            plt.show()
            plt.savefig("btstch.png")
                
    except Exception as e:
        print(f"An error occurred during file operations: {e}")

    print("Program finished. Results written to tstch.dat")

# Call the main program function
if __name__ == '__main__':
    test_program()
