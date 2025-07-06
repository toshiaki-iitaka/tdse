import numpy as np
import math
from scipy.special import jv # For Bessel function of the first kind
import matplotlib.pyplot as plt

# --- Parameters ---
PI = math.pi
NMAX = 100  # Maximum number of rooms/points

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
    p_cheb = np.zeros(nmax_val, dtype=np.complex128)
    q_cheb = np.zeros(nmax_val, dtype=np.complex128)
    r_cheb = np.zeros(nmax_val, dtype=np.complex128) # Temporary array for calculations
    ph_cheb = np.zeros(nmax_val, dtype=np.complex128) # Resulting propagated wave function

    # P(N)=PH0(N)
    for n_idx in range(nmax_val):
        p_cheb[n_idx] = ph0[n_idx]

    # Q(1) = -IU*(HA(1)*P(1)+HB(1)*P(2))
    q_cheb[0] = -iu * (ha_norm[0] * p_cheb[0] + hb_norm[0] * p_cheb[1])
    # DO 30 N=2,NMAX-1
    # Q(N) = -IU*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))
    # 30 CONTINUE
    for n_idx in range(1, nmax_val - 1):
        q_cheb[n_idx] = -iu * (hb_norm[n_idx-1] * p_cheb[n_idx-1] + ha_norm[n_idx] * p_cheb[n_idx] + hb_norm[n_idx] * p_cheb[n_idx+1])
    # Q(NMAX)= -IU*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))
    q_cheb[nmax_val-1] = -iu * (hb_norm[nmax_val-2] * p_cheb[nmax_val-2] + ha_norm[nmax_val-1] * p_cheb[nmax_val-1])

    # Initial terms of the Chebyshev expansion
    # AN=DJBES(0,EGRID*T/2)
    an_0 = jv(0, egrid * t / 2.0)
    # DO 40 N=1,NMAX
    # PH(N)=AN*P(N)
    # 40 CONTINUE
    for n_idx in range(nmax_val):
        ph_cheb[n_idx] = an_0 * p_cheb[n_idx]

    # AN=2*DJBES(1,EGRID*T/2)
    an_1 = 2 * jv(1, egrid * t / 2.0)
    # DO 45 N=1,NMAX
    # PH(N)=PH(N)+AN*Q(N)
    # 45 CONTINUE
    for n_idx in range(nmax_val):
        ph_cheb[n_idx] = ph_cheb[n_idx] + an_1 * q_cheb[n_idx]

    # Main Chebyshev series loop
    # DO 1000 I=2, IMAX, 3
    # The loop structure is a bit tricky due to Fortran's GOTO.
    # We will break the loop if convergence is met.
    for i in range(2, imax + 1, 3):
        # AN=2*DJBES(I,EGRID*T/2)
        an_i = 2 * jv(i, egrid * t / 2.0)
        # Check for convergence based on the current Bessel term
        if abs(an_i) < eps:
            break # Exit loop if Bessel function term is negligible

        # R(1) = -2*IU*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
        r_cheb[0] = -2 * iu * (ha_norm[0] * q_cheb[0] + hb_norm[0] * q_cheb[1]) + p_cheb[0]
        # PH(1)=PH(1) +AN*R(1)
        ph_cheb[0] = ph_cheb[0] + an_i * r_cheb[0]
        # DO 70 N=2,NMAX-1
        # R(N)=-2*IU*(HB(N-1)*Q(N-1)+HA(N)*Q(N)+HB(N)*Q(N+1))+P(N)
        # PH(N)=PH(N)+AN*R(N)
        # 70 CONTINUE
        for n_idx in range(1, nmax_val - 1):
            r_cheb[n_idx] = -2 * iu * (hb_norm[n_idx-1] * q_cheb[n_idx-1] + ha_norm[n_idx] * q_cheb[n_idx] + hb_norm[n_idx] * q_cheb[n_idx+1]) + p_cheb[n_idx]
            ph_cheb[n_idx] = ph_cheb[n_idx] + an_i * r_cheb[n_idx]
        # R(NMAX)=-2*IU*(HB(NMAX-1)*Q(NMAX-1)+HA(NMAX)*Q(NMAX))+P(NMAX)
        r_cheb[nmax_val-1] = -2 * iu * (hb_norm[nmax_val-2] * q_cheb[nmax_val-2] + ha_norm[nmax_val-1] * q_cheb[nmax_val-1]) + p_cheb[nmax_val-1]
        # PH(NMAX)=PH(NMAX)+AN*R(NMAX)
        ph_cheb[nmax_val-1] = ph_cheb[nmax_val-1] + an_i * r_cheb[nmax_val-1]

        # AN=2*DJBES(I+1,EGRID*T/2)
        an_i_plus_1 = 2 * jv(i + 1, egrid * t / 2.0)
        if abs(an_i_plus_1) < eps:
            break

        # P(1) = -2*IU*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
        p_temp_0 = -2 * iu * (ha_norm[0] * r_cheb[0] + hb_norm[0] * r_cheb[1]) + q_cheb[0]
        # PH(1)=PH(1) +AN*P(1)
        ph_cheb[0] = ph_cheb[0] + an_i_plus_1 * p_temp_0
        # DO 71 N=2,NMAX-1
        # P(N)=-2*IU*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))+Q(N)
        # PH(N)=PH(N)+AN*P(N)
        # 71 CONTINUE
        for n_idx in range(1, nmax_val - 1):
            p_temp_n = -2 * iu * (hb_norm[n_idx-1] * r_cheb[n_idx-1] + ha_norm[n_idx] * r_cheb[n_idx] + hb_norm[n_idx] * r_cheb[n_idx+1]) + q_cheb[n_idx]
            ph_cheb[n_idx] = ph_cheb[n_idx] + an_i_plus_1 * p_temp_n
            p_cheb[n_idx] = p_temp_n # Update p_cheb for next iteration
        # P(NMAX)=-2*IU*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))+Q(NMAX)
        p_temp_nmax = -2 * iu * (hb_norm[nmax_val-2] * r_cheb[nmax_val-2] + ha_norm[nmax_val-1] * r_cheb[nmax_val-1]) + q_cheb[nmax_val-1]
        ph_cheb[nmax_val-1] = ph_cheb[nmax_val-1] + an_i_plus_1 * p_temp_nmax
        p_cheb[nmax_val-1] = p_temp_nmax # Update p_cheb for next iteration
        p_cheb[0] = p_temp_0 # Update p_cheb[0]

        # AN=2*DJBES(I+2,EGRID*T/2)
        an_i_plus_2 = 2 * jv(i + 2, egrid * t / 2.0)
        if abs(an_i_plus_2) < eps:
            break

        # Q(1) = -2*IU*(HA(1)*P(1)+HB(1)*P(2))+R(1)
        q_temp_0 = -2 * iu * (ha_norm[0] * p_cheb[0] + hb_norm[0] * p_cheb[1]) + r_cheb[0]
        # PH(1)=PH(1) +AN*Q(1)
        ph_cheb[0] = ph_cheb[0] + an_i_plus_2 * q_temp_0
        # DO 72 N=2,NMAX-1
        # Q(N)=-2*IU*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))+R(N)
        # PH(N)=PH(N)+AN*Q(N)
        # 72 CONTINUE
        for n_idx in range(1, nmax_val - 1):
            q_temp_n = -2 * iu * (hb_norm[n_idx-1] * p_cheb[n_idx-1] + ha_norm[n_idx] * p_cheb[n_idx] + hb_norm[n_idx] * p_cheb[n_idx+1]) + r_cheb[n_idx]
            ph_cheb[n_idx] = ph_cheb[n_idx] + an_i_plus_2 * q_temp_n
            q_cheb[n_idx] = q_temp_n # Update q_cheb for next iteration
        # Q(NMAX)=-2*IU*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))+R(NMAX)
        q_temp_nmax = -2 * iu * (hb_norm[nmax_val-2] * p_cheb[nmax_val-2] + ha_norm[nmax_val-1] * p_cheb[nmax_val-1]) + r_cheb[nmax_val-1]
        ph_cheb[nmax_val-1] = ph_cheb[nmax_val-1] + an_i_plus_2 * q_temp_nmax
        q_cheb[nmax_val-1] = q_temp_nmax # Update q_cheb for next iteration
        q_cheb[0] = q_temp_0 # Update q_cheb[0]

    else: # This block executes if the loop completes without a 'break'
        print('CHEBYSHEV NOT CONVERGE')

    # Final transformation
    # DO 80 N=1,NMAX
    # PH(N)=PH(N)*EXP(-IU*(EGRID/2+EMIN)*T)
    # 80 CONTINUE
    for n_idx in range(nmax_val):
        ph_cheb[n_idx] = ph_cheb[n_idx] * np.exp(-iu * (egrid / 2.0 + emin) * t)

    return ph_cheb

# --- Main Program: TEST ---
def test_program():
    """
    Main program to simulate 1D quantum mechanics using the Chebyshev scheme.
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # Get input from user for M (number of time steps to iterate, though only 50 are used)
    print('M')
    try:
        # The Fortran code reads 'M' but then iterates 50 times regardless.
        # So, 'M' is not actually used in the loop limit.
        # We'll still read it for fidelity to the original.
        m_input = int(input())
    except ValueError:
        print("Invalid input. Please enter an integer for M.")
        return

    # --- Initial Data Calculations ---
    dx = 1.0
    dx2 = dx * dx

    # --- Calculation of Hamiltonian ---
    ha = np.zeros(NMAX, dtype=np.float64)
    hb = np.zeros(NMAX - 1, dtype=np.float64)
    v_potential = np.zeros(NMAX, dtype=np.float64) # V from Fortran

    # DO 50 N=1,NMAX
    # IF(45.LE.N .AND. N.LE.55) THEN
    #   V(N)=0.05*COS((I-50)/5.0 * PI/2)**2  <- Assuming 'I' should be 'N'
    # ELSE
    #   V(N)=0
    # ENDIF
    # HA(N) = V(N)
    # 50 CONTINUE
    for n_idx in range(NMAX):
        # Fortran N is (n_idx + 1) in Python
        if 45 <= (n_idx + 1) <= 55:
            # (I-50) in Fortran is ((n_idx+1)-50) in Python
            v_potential[n_idx] = 0.05 * (np.cos(((n_idx + 1) - 50) / 5.0 * PI / 2.0))**2
        else:
            v_potential[n_idx] = 0.0
        ha[n_idx] = v_potential[n_idx] # HA(N) = V(N)

    # DO 51 N=1,NMAX-1
    # HB(N) = -0.5/DX2
    # 51 CONTINUE
    for n_idx in range(NMAX - 1):
        hb[n_idx] = -0.5 / dx2

    # --- Normalization of Hamiltonian ---
    # VMAX, VMIN calculation based on V(N)
    # The Fortran code then hardcodes VMAX=1.0, overriding the calculation.
    # vmax_calc = np.max(v_potential) if v_potential.size > 0 else 0.0
    # vmin_calc = np.min(v_potential) if v_potential.size > 0 else 0.0

    # Fortran's explicit override
    vmax_actual = 1.0
    print(f'VMAX={vmax_actual:15.7E}')

    # EMAX, EMIN, EGRID for Chebyshev normalization
    # Fortran: EMAX=1/DX2 +VMAX
    # Fortran: EMIN=-1/DX2
    emax_cheb = 1.0 / dx2 + vmax_actual
    emin_cheb = -1.0 / dx2
    egrid = emax_cheb - emin_cheb

    han_norm = np.zeros(NMAX, dtype=np.float64)
    hbn_norm = np.zeros(NMAX - 1, dtype=np.float64)

    # DO 20 N=1,NMAX
    # HAN(N)=(HA(N)-EGRID/2-EMIN)*2/EGRID
    # IF(N.LT.NMAX) HBN(N)=HB(N)*2/EGRID
    # 20 CONTINUE
    for n_idx in range(NMAX):
        han_norm[n_idx] = (ha[n_idx] - egrid / 2.0 - emin_cheb) * 2.0 / egrid
        if n_idx < NMAX - 1:
            hbn_norm[n_idx] = hb[n_idx] * 2.0 / egrid

    print(f'EGRID={egrid:15.7E}')
    print('TIME,M') # 'M' is just a label in Fortran output, not a variable

    # Open the output files
    try:
        with open('wave.dat', 'w') as f_wave_out, \
             open('scch.dat', 'w') as f_scch_out:
            plot_x = range(1, 51)
            plot_y = [[] for _ in range(6) ]
            # DO 1000 ITIME=1,50
            for itime in range(1, 51): # Loop from 1 to 50
                # T=2*10**(ITIME/25.0)
                t_val = 2 * 10**(itime / 25.0)
                print(f'T={t_val:15.7E}')
                print(f'ESTIMATED STEP={egrid*t_val/2:15.7E}')

                # --- Setting Initial Wave Function (Gaussian Wave Packet) ---
                ph0 = np.zeros(NMAX, dtype=np.complex128)
                
                # Fortran:
                # X0=0.25*DX*NMAX
                # SG=0.1*DX*NMAX
                # PA=0.1*2*PI/DX
                x0 = 0.25 * dx * NMAX
                sg = 0.1 * dx * NMAX
                pa_momentum = 0.1 * 2 * PI / dx # Renamed to avoid conflict with normalization PA

                # DO 60 N=1,NMAX
                # X=N*DX
                # PH0(N)= EXP(-0.5*(X-X0)**2/SG**2) * EXP(IU*PA*(X-X0))
                # 60 CONTINUE
                for n_idx in range(NMAX):
                    x = (n_idx + 1) * dx # Fortran N*DX
                    ph0[n_idx] = np.exp(-0.5 * (x - x0)**2 / sg**2) * np.exp(iu * pa_momentum * (x - x0))

                # --- Normalization of Wave Function ---
                pa_norm_factor = 0.0
                # DO 101 N=1,NMAX
                # PA = PA + ABS(PH0(N))**2
                # 101 CONTINUE
                for n_idx in range(NMAX):
                    pa_norm_factor += abs(ph0[n_idx])**2
                pa_norm_factor = np.sqrt(pa_norm_factor)

                # DO 110 N=1,NMAX
                # PH0(N)=PH0(N)/PA
                # 110 CONTINUE
                for n_idx in range(NMAX):
                    ph0[n_idx] = ph0[n_idx] / pa_norm_factor

                # --- Call Chebyshev Subroutine ---
                # CALL CH(HAN,HBN,EGRID,EMIN,PH0,T,PH)
                ph_result = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, t_val, NMAX)

                # --- Output Wave Function and Normalization Check ---
                anorm = 0.0
                # OVR is commented out in Fortran, so we won't calculate it here.
                # OVR=0.0

                # DO 150 N=1,NMAX
                # ANORM=ANORM+ABS(PH(N))**2
                # WRITE(21,100) N,ABS(PH(N))**2+ 0.01*ITIME
                # 150 CONTINUE
                for n_idx in range(NMAX):
                    anorm += abs(ph_result[n_idx])**2
                    # Fortran FORMAT(1H ,I4,30E15.7)
                    f_wave_out.write(f"{n_idx + 1:4d} {abs(ph_result[n_idx])**2 + 0.01 * itime:15.7E}\n")

                # WRITE(*,100) ITIME,EMAX*T,ANORM-1
                # WRITE(20,100) ITIME,EMAX*T,ABS(ANORM-1)
                console_output = f"{itime:4d} {emax_cheb*t_val:15.7E} {anorm-1:15.7E}"
                f_scch_output = f"{itime:4d} {emax_cheb*t_val:15.7E} {abs(anorm-1):15.7E}"
                
                print(console_output)
                f_scch_out.write(f_scch_output + "\n")

                plot_y[0].append(itime)
                plot_y[1].append(emax_cheb*t_val)
                plot_y[2].append(abs(anorm-1))
                #plot_y[3].append(abs(1.0 - abs(ovr)))
                #plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))
                #plot_y[5].append(abs(em*dt)**3 / 6.0 * itime)

            #
            # plot
            #
            plt.gca().clear()
            plt.title('SCCH')
            plt.xlabel(r'Emax*t')
            plt.ylabel(r'error')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlim([1,100])
            plt.ylim([1e-20,1])
            plt.plot(plot_y[1],plot_y[2], marker='o', linestyle='None')
            #plt.plot(plot_y[1],plot_y[4], marker='o', linestyle='None')
            #plt.plot(plot_y[1],plot_y[5], marker='None', linestyle='dashed')
            plt.show()
            plt.savefig("scch.png")

    except Exception as e:
        print(f"An error occurred during file operations: {e}")

    print("Program finished. Results written to wave.dat and scch.dat")

# Call the main program function
if __name__ == '__main__':
    test_program()
    
