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

# --- Main Program: SCCN ---
def sccn_program():
    """
    Main program to simulate 1D quantum mechanics using the Crank-Nicolson scheme.
    It also uses the Chebyshev scheme for comparison (exact solution).
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # --- Initial Data ---
    # Open files
    try:
        f_sccn_out = open('sccn.dat', 'w')
    except Exception as e:
        print(f"Error opening sccn.dat: {e}")
        return

    print('ALPHA')
    try:
        alpha = float(input())
    except ValueError:
        print("Invalid input. Please enter a number for ALPHA.")
        f_sccn_out.close()
        return

    dx = 1.0
    dx2 = dx * dx
    emax_initial = 1.0 / dx2 + 1.0
    dt = alpha / emax_initial
    
    mtime = (2 * PI / alpha) * 40
    mtime = min(int(mtime), 100000) # Ensure MTIME is an integer

    # --- Calculation of Hamiltonian ---
    ha = np.zeros(NMAX, dtype=np.float64)
    # HB is NMAX elements in SCCN main program for A1, B1, A2, B2 calculations
    hb = np.zeros(NMAX, dtype=np.float64) 
    v_potential = np.zeros(NMAX, dtype=np.float64) # V from Fortran

    # A1, B1, A2, B2 arrays for Crank-Nicolson coefficients
    a1_cn = np.zeros(NMAX, dtype=np.complex128)
    b1_cn = np.zeros(NMAX, dtype=np.complex128)
    a2_cn = np.zeros(NMAX, dtype=np.complex128)
    b2_cn = np.zeros(NMAX, dtype=np.complex128)

    # DO 50 N=1,NMAX
    # IF(45.LE.N .AND. N.LE.55) THEN
    #   V(N)=0.05*COS((I-50)/5.0 * PI/2)**2  <- Assuming 'I' should be 'N'
    # ELSE
    #   V(N)=0
    # ENDIF
    # HA(N) = V(N)
    # HB(N) = -0.5/DX2
    # A1(N) = 1.-0.5*IU*DT*HA(N)
    # A2(N) = 1.+0.5*IU*DT*HA(N)
    # B1(N) = -0.5*IU*DT*HB(N)
    # B2(N) = +0.5*IU*DT*HB(N)
    # 50 CONTINUE
    for n_idx in range(NMAX):
        # Fortran N is (n_idx + 1) in Python
        if 45 <= (n_idx + 1) <= 55:
            # (I-50) in Fortran is ((n_idx+1)-50) in Python
            v_potential[n_idx] = 0.05 * (np.cos(((n_idx + 1) - 50) / 5.0 * PI / 2.0))**2
        else:
            v_potential[n_idx] = 0.0
        ha[n_idx] = v_potential[n_idx]
        hb[n_idx] = -0.5 / dx2 # HB is NMAX elements here

        a1_cn[n_idx] = 1.0 - 0.5 * iu * dt * ha[n_idx]
        a2_cn[n_idx] = 1.0 + 0.5 * iu * dt * ha[n_idx]
        b1_cn[n_idx] = -0.5 * iu * dt * hb[n_idx]
        b2_cn[n_idx] = +0.5 * iu * dt * hb[n_idx]

    # --- Normalization of Hamiltonian (for Chebyshev subroutine) ---
    vmax_actual = 1.0 # Fortran's explicit override
    print(f'VMAX={vmax_actual:15.7E}')

    emax_cheb = 1.0 / dx2 + vmax_actual
    emin_cheb = -1.0 / dx2
    egrid = emax_cheb - emin_cheb

    han_norm = np.zeros(NMAX, dtype=np.float64)
    # HBN for CH subroutine is NMAX-1 elements
    hbn_norm = np.zeros(NMAX - 1, dtype=np.float64) 

    # DO 20 N=1,NMAX
    # HAN(N)=(HA(N)-EGRID/2-EMIN)*2/EGRID
    # IF(N.LT.NMAX) HBN(N)=HB(N)*2/EGRID
    # 20 CONTINUE
    for n_idx in range(NMAX):
        han_norm[n_idx] = (ha[n_idx] - egrid / 2.0 - emin_cheb) * 2.0 / egrid
        if n_idx < NMAX - 1: # HBN is NMAX-1 elements
            hbn_norm[n_idx] = hb[n_idx] * 2.0 / egrid

    print(f'EGRID={egrid:15.7E}')
    print('TIME,M')

    # --- Setting Initial Wave Function (Gaussian Wave Packet) ---
    ph0 = np.zeros(NMAX, dtype=np.complex128)
    
    x0 = 0.25 * dx * NMAX
    sg = 0.1 * dx * NMAX
    pa_momentum = 0.1 * 2 * PI / dx

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
    
    # P array for Crank-Nicolson time evolution
    p_cn = np.zeros(NMAX, dtype=np.complex128)
    # DO 81 N=1,NMAX
    # P(N)=PH0(N)
    # 81 CONTINUE
    for n_idx in range(NMAX):
        p_cn[n_idx] = ph0[n_idx]

    # Arrays for tridiagonal solver
    o_cn = np.zeros(NMAX, dtype=np.complex128)
    e_cn = np.zeros(NMAX, dtype=np.complex128)
    f_cn = np.zeros(NMAX, dtype=np.complex128)

    # --- Time Evolution Loop (Crank-Nicolson) ---
    # DO 1000 ITIME=0,MTIME
    plot_x = range(mtime + 1)
    plot_y = [[] for _ in range(6) ]
    for itime in range(mtime + 1): # Loop from 0 to MTIME inclusive
        t_val = itime * dt

        # --- Output Wave Function (every 30 steps) ---
        # IF(MOD(ITIME,30).EQ. 0) THEN
        if (itime % 30) == 0:
            # CALL CH(HAN,HBN,EGRID,EMIN,PH0,T,PH)
            ph_exact = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, t_val, NMAX)

            anorm = 0.0
            ovr = 0.0 + 0.0j
            # DO 150 N=1,NMAX
            # ANORM=ANORM+ABS(P(N))**2
            # OVR=OVR+DCONJG(P(N))*PH(N)
            # 150 CONTINUE
            for n_idx in range(NMAX):
                anorm += abs(p_cn[n_idx])**2
                ovr += np.conjugate(p_cn[n_idx]) * ph_exact[n_idx]

            # WRITE(*,100) ITIME,EMAX*T,ANORM-1,
            # &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
            # WRITE(20,100) ITIME,EMAX*T,ABS(ANORM-1),
            # &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
            
            # Format string for Fortran's 30E15.7 (adjusted for Python's f-strings)
            # For 5 values: I7, 4x E15.7
            console_output = (
                f"{itime:7d} {emax_initial*t_val:15.7E} {anorm-1:15.7E} "
                f"{1.0 - abs(ovr):15.7E} {math.atan2(ovr.imag, ovr.real):15.7E}"
            )
            file_output = (
                f"{itime:7d} {emax_initial*t_val:15.7E} {abs(anorm-1):15.7E} "
                f"{abs(1.0 - abs(ovr)):15.7E} {abs(math.atan2(ovr.imag, ovr.real)):15.7E}"
            )
            print(console_output)
            f_sccn_out.write(file_output + "\n")
            
            plot_y[0].append(itime)
            plot_y[1].append(emax_initial*t_val)
            plot_y[2].append(abs(anorm-1))
            plot_y[3].append(abs(1.0 - abs(ovr)))
            plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))
            #plot_y[5].append(abs(em*dt)**3 / 6.0 * itime)

        # --- Crank-Nicolson Time Evolution Step ---
        # Multiplication of Symmetric Tridiagonal Matrix: O = (AB)*P
        # O(1) = A1(1)*P(1)+B1(1)*P(2)
        o_cn[0] = a1_cn[0] * p_cn[0] + b1_cn[0] * p_cn[1]
        # DO 70 N=2,NMAX-1
        # O(N) = B1(N-1)*P(N-1)+A1(N)*P(N)+B1(N)*P(N+1)
        # 70 CONTINUE
        for n_idx in range(1, NMAX - 1): # Python index n_idx corresponds to Fortran N
            o_cn[n_idx] = b1_cn[n_idx-1] * p_cn[n_idx-1] + a1_cn[n_idx] * p_cn[n_idx] + b1_cn[n_idx] * p_cn[n_idx+1]
        # O(NMAX) = B1(NMAX-1)*P(NMAX-1)+A1(NMAX)*P(NMAX)
        o_cn[NMAX-1] = b1_cn[NMAX-2] * p_cn[NMAX-2] + a1_cn[NMAX-1] * p_cn[NMAX-1]

        # Inversion of Symmetric Tridiagonal Matrix (Thomas Algorithm): P = O / (AB)
        # E(1) = -A2(1)/B2(1)
        e_cn[0] = -a2_cn[0] / b2_cn[0]
        # F(1) = O(1)/B2(1)
        f_cn[0] = o_cn[0] / b2_cn[0]

        # DO 80 N=2,NMAX
        # E(N)= -(B2(N-1)/E(N-1)+A2(N)) / B2(N)
        # F(N)= O(N)/B2(N) + (B2(N-1)/B2(N))*(F(N-1)/E(N-1))
        # 80 CONTINUE
        for n_idx in range(1, NMAX): # Python index n_idx corresponds to Fortran N
            # Fortran's B2(N) is b2_cn[n_idx]
            # Fortran's B2(N-1) is b2_cn[n_idx-1]
            # Fortran's E(N-1) is e_cn[n_idx-1]
            # Fortran's A2(N) is a2_cn[n_idx]
            # Fortran's O(N) is o_cn[n_idx]
            # Fortran's F(N-1) is f_cn[n_idx-1]

            # Check for division by zero for E(N-1)
            if abs(e_cn[n_idx-1]) < 1e-20: # Small epsilon to avoid near-zero division
                print(f"Warning: Division by near-zero E({n_idx}) at ITIME={itime}. Numerical instability possible.")
                # Handle this case, e.g., by breaking or returning an error,
                # or by using a very small number instead of 0.
                # For direct translation, we'll let it proceed as Fortran might.
                # In real applications, this would need robust handling.

            e_cn[n_idx] = -(b2_cn[n_idx-1] / e_cn[n_idx-1] + a2_cn[n_idx]) / b2_cn[n_idx]
            f_cn[n_idx] = o_cn[n_idx] / b2_cn[n_idx] + (b2_cn[n_idx-1] / b2_cn[n_idx]) * (f_cn[n_idx-1] / e_cn[n_idx-1])

        # P(NMAX)=-F(NMAX)/E(NMAX)
        p_cn[NMAX-1] = -f_cn[NMAX-1] / e_cn[NMAX-1]
        
        # DO 90 N=NMAX-1,1,-1
        # P(N)=(P(N+1)-F(N))/E(N)
        # 90 CONTINUE
        for n_idx in range(NMAX - 2, -1, -1): # Loop from NMAX-2 down to 0
            p_cn[n_idx] = (p_cn[n_idx+1] - f_cn[n_idx]) / e_cn[n_idx]

    #
    # plot
    #
    plt.gca().clear()
    plt.title('SCCN')
    plt.xlabel(r'Emax*t')
    plt.ylabel(r'error')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1,100])
    plt.ylim([1e-20,1])
    plt.plot(plot_y[1],plot_y[2], marker='o', linestyle='None')
    plt.plot(plot_y[1],plot_y[4], marker='o', linestyle='None')
    #plt.plot(plot_y[1],plot_y[5], marker='None', linestyle='dashed')
    plt.show()
    plt.savefig("sccn.png")

    # Close the output file
    f_sccn_out.close()
    print("Program finished. Results written to sccn.dat")

# Call the main program function
if __name__ == '__main__':
    sccn_program()
