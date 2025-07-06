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

# --- Main Program: NROOM ---
def nroom_program():
    """
    Main program to simulate N-room quantum mechanics with a potential well
    and a combination of Chebyshev and 7-step time evolution.
    """
    # Define complex imaginary unit
    iu = 0.0 + 1.0j

    # Get input from user
    print('ALPHA')
    try:
        alpha = float(input())
    except ValueError:
        print("Invalid input. Please enter a number for ALPHA.")
        return

    # --- Initial Data Calculations ---
    dx = 1.0
    dx2 = dx * dx
    # EMAX from Fortran: EMAX=1/DX2+1
    emax_initial = 1.0 / dx2 + 1.0
    dt = alpha / emax_initial
    
    # Fortran: MTIME=(2*PI/ALPHA) * 40
    # Fortran: MTIME=MIN(MTIME,100000)
    mtime = (2 * PI / alpha) * 40
    mtime = min(int(mtime), 100000) # Ensure MTIME is an integer

    # --- Calculation of Hamiltonian ---
    ha = np.zeros(NMAX, dtype=np.float64)
    # HB is NMAX-1 elements for off-diagonal terms
    hb = np.zeros(NMAX - 1, dtype=np.float64)
    v0_potential = np.zeros(NMAX, dtype=np.float64) # V0 from Fortran

    # DO 50 N=1,NMAX
    # IF(45.LE.N .AND. N.LE.55) THEN
    #   V0(N)=0.05*COS((I-50)/5.0 * PI/2)**2
    # ELSE
    #   V0(N)=0
    # ENDIF
    # HA(N) = V0(N)
    # 50 CONTINUE
    for n_idx in range(NMAX):
        # Fortran I is N here, so (n_idx+1)
        if 45 <= (n_idx + 1) <= 55:
            # (I-50) in Fortran is (n_idx+1-50) in Python
            v0_potential[n_idx] = 0.05 * (np.cos(((n_idx + 1) - 50) / 5.0 * PI / 2.0))**2
        else:
            v0_potential[n_idx] = 0.0
        ha[n_idx] = v0_potential[n_idx] # HA(N) = V0(N)

    # DO 51 N=1,NMAX-1
    # HB(N) = -0.5/DX2
    # 51 CONTINUE
    for n_idx in range(NMAX - 1):
        hb[n_idx] = -0.5 / dx2

    # --- Normalization of Hamiltonian ---
    # VMAX, VMIN calculation based on V0(N)
    # The Fortran code then hardcodes VMAX=1.0, overriding the calculation.
    # vmax_calc = np.max(v0_potential) if v0_potential.size > 0 else 0.0
    # vmin_calc = np.min(v0_potential) if v0_potential.size > 0 else 0.0

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
    # HBN is NMAX-1 elements to match HB in CH subroutine
    hbn_norm = np.zeros(NMAX - 1, dtype=np.float64)

    # DO 20 N=1,NMAX
    # HAN(N)=(HA(N)-EGRID/2-EMIN)*2/EGRID
    # IF(N.LT.NMAX) HBN(N)=HB(N)*2/EGRID
    # 20 CONTINUE
    for n_idx in range(NMAX):
        han_norm[n_idx] = (ha[n_idx] - egrid / 2.0 - emin_cheb) * 2.0 / egrid
        # Only assign HBN up to NMAX-1, as HB is NMAX-1 elements
        if n_idx < NMAX - 1:
            hbn_norm[n_idx] = hb[n_idx] * 2.0 / egrid

    print(f'EGRID={egrid:15.7E}')
    #print('TIME,M') # M is not read in this program, it's just a print from Fortran

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

    # --- Normalization of Initial Wave Function ---
    pa_norm_factor = 0.0 # Renamed to avoid conflict with momentum PA
    # DO 101 N=1,NMAX
    # PA = PA + ABS(PH0(N))**2
    # 101 CONTINUE
    for n_idx in range(NMAX):
        pa_norm_factor += abs(ph0[n_idx])**2
    pa_norm_factor = np.sqrt(pa_norm_factor)

    # Initialize P, Q, R, S, T, U, V for the 7-step time evolution
    # These will represent psi(t), psi(t+dt), psi(t+2dt), ..., psi(t+6dt)
    p_psi = np.zeros(NMAX, dtype=np.complex128)
    q_psi = np.zeros(NMAX, dtype=np.complex128)
    r_psi = np.zeros(NMAX, dtype=np.complex128)
    s_psi = np.zeros(NMAX, dtype=np.complex128)
    t_psi = np.zeros(NMAX, dtype=np.complex128)
    u_psi = np.zeros(NMAX, dtype=np.complex128)
    v_psi = np.zeros(NMAX, dtype=np.complex128) 

    # DO 110 N=1,NMAX
    # PH0(N)=PH0(N)/PA
    # P(N)=PH0(N)
    # 110 CONTINUE
    for n_idx in range(NMAX):
        ph0[n_idx] = ph0[n_idx] / pa_norm_factor
        p_psi[n_idx] = ph0[n_idx] # P is psi(t=0)

    # Initial calls to CH for the first few time steps
    q_psi = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, dt, NMAX)
    r_psi = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, 2 * dt, NMAX)
    s_psi = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, 3 * dt, NMAX)
    t_psi = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, 4 * dt, NMAX)
    u_psi = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, 5 * dt, NMAX)
    v_psi = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, 6 * dt, NMAX)


    # --- Time Evolution Loop (7-step method) ---
    # Open the output file
    try:
        with open('sc6.dat', 'w') as f_out:
            plot_x = range(5, mtime + 1, 7)
            plot_y = [[] for _ in range(6) ]
            # DO 1000 ITIME=5,MTIME,7
            for itime in range(5, mtime + 1, 7): # Loop starts at 5, increments by 7
                tt = dt * itime

                # --- Output Wave Function (every 28th step, starting from 5) ---
                # IF(MOD(ITIME,28).EQ. 5) THEN
                if (itime % 28) == 5:
                    # CALL CH(HAN,HBN,EGRID,EMIN,PH0,TT,PH)
                    # PH is the exact solution at time TT, calculated by Chebyshev
                    ph_exact = ch_subroutine(han_norm, hbn_norm, egrid, emin_cheb, ph0, tt, NMAX)

                    anorm = 0.0
                    ovr = 0.0 + 0.0j
                    # DO 150 N=1,NMAX
                    # ANORM=ANORM+ABS(U(N))**2
                    # OVR=OVR+DCONJG(U(N))*PH(N)
                    # 150   CONTINUE
                    for n_idx in range(NMAX):
                        anorm += abs(u_psi[n_idx])**2 # U is the current state for output
                        ovr += np.conjugate(u_psi[n_idx]) * ph_exact[n_idx]

                    # Fortran console output is commented out.
                    # c WRITE(*,100) ITIME,EMAX*TT,ANORM-1,
                    # c&1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
                    # WRITE(20,100) ITIME,EMAX*TT,ABS(ANORM-1),
                    # &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
                    file_output = (
                        f"{itime:7d} {emax_initial*tt:15.7E} {abs(anorm-1):15.7E} "
                        f"{abs(1.0 - abs(ovr)):15.7E} {abs(np.arctan2(ovr.imag, ovr.real)):15.7E}"
                    )
                    f_out.write(file_output + "\n")

                    plot_y[0].append(itime)
                    plot_y[1].append(emax_initial*tt)
                    plot_y[2].append(abs(anorm-1))
                    plot_y[3].append(abs(1.0 - abs(ovr)))
                    plot_y[4].append(abs(abs(np.arctan2(ovr.imag, ovr.real))))
                    #plot_y[5].append(abs(em*dt)**3 / 6.0 * itime)

                # --- Time Evolution (7-step method) ---
                # These are in-place updates. The variables P, Q, R, S, T, U, V
                # are updated sequentially, with each step using the most recently
                # updated values of the variables.

                # Temporary array for W calculation
                w_temp_calc = np.zeros(NMAX, dtype=np.complex128)

                # Step 1: Calculate new V(N)
                # W(N) = (11*Q(N)-14*R(N)+26*S(N)-14*T(N)+11*U(N))/20.0
                # Uses current Q, R, S, T, U
                w_temp_calc[:] = (11 * q_psi - 14 * r_psi + 26 * s_psi - 14 * t_psi + 11 * u_psi) / 20.0
                
                # V(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + P(1)
                v_psi[0] = -6 * iu * dt * (ha[0] * w_temp_calc[0] + hb[0] * w_temp_calc[1]) + p_psi[0]
                # DO 70 N=2,NMAX-1
                # V(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + P(N)
                for n_idx in range(1, NMAX - 1):
                    v_psi[n_idx] = -6 * iu * dt * (hb[n_idx-1] * w_temp_calc[n_idx-1] + ha[n_idx] * w_temp_calc[n_idx] + hb[n_idx] * w_temp_calc[n_idx+1]) + p_psi[n_idx]
                # V(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + P(NMAX)
                v_psi[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w_temp_calc[NMAX-2] + ha[NMAX-1] * w_temp_calc[NMAX-1]) + p_psi[NMAX-1]

                # Step 2: Calculate new P(N)
                # W(N) = (11*R(N)-14*S(N)+26*T(N)-14*U(N)+11*V(N))/20.0
                # Uses current R, S, T, U, and the NEW V
                w_temp_calc[:] = (11 * r_psi - 14 * s_psi + 26 * t_psi - 14 * u_psi + 11 * v_psi) / 20.0
                
                # P(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + Q(1)
                p_psi[0] = -6 * iu * dt * (ha[0] * w_temp_calc[0] + hb[0] * w_temp_calc[1]) + q_psi[0]
                # DO 71 N=2,NMAX-1
                # P(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + Q(N)
                for n_idx in range(1, NMAX - 1):
                    p_psi[n_idx] = -6 * iu * dt * (hb[n_idx-1] * w_temp_calc[n_idx-1] + ha[n_idx] * w_temp_calc[n_idx] + hb[n_idx] * w_temp_calc[n_idx+1]) + q_psi[n_idx]
                # P(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + Q(NMAX)
                p_psi[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w_temp_calc[NMAX-2] + ha[NMAX-1] * w_temp_calc[NMAX-1]) + q_psi[NMAX-1]

                # Step 3: Calculate new Q(N)
                # W(N) = (11*S(N)-14*T(N)+26*U(N)-14*V(N)+11*P(N))/20.0
                # Uses current S, T, U, and the NEW V, P
                w_temp_calc[:] = (11 * s_psi - 14 * t_psi + 26 * u_psi - 14 * v_psi + 11 * p_psi) / 20.0
                
                # Q(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + R(1)
                q_psi[0] = -6 * iu * dt * (ha[0] * w_temp_calc[0] + hb[0] * w_temp_calc[1]) + r_psi[0]
                # DO 72 N=2,NMAX-1
                # Q(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + R(N)
                for n_idx in range(1, NMAX - 1):
                    q_psi[n_idx] = -6 * iu * dt * (hb[n_idx-1] * w_temp_calc[n_idx-1] + ha[n_idx] * w_temp_calc[n_idx] + hb[n_idx] * w_temp_calc[n_idx+1]) + r_psi[n_idx]
                # Q(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + R(NMAX)
                q_psi[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w_temp_calc[NMAX-2] + ha[NMAX-1] * w_temp_calc[NMAX-1]) + r_psi[NMAX-1]

                # Step 4: Calculate new R(N)
                # W(N) = (11*T(N)-14*U(N)+26*V(N)-14*P(N)+11*Q(N))/20.0
                # Uses current T, U, and the NEW V, P, Q
                w_temp_calc[:] = (11 * t_psi - 14 * u_psi + 26 * v_psi - 14 * p_psi + 11 * q_psi) / 20.0
                
                # R(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + S(1)
                r_psi[0] = -6 * iu * dt * (ha[0] * w_temp_calc[0] + hb[0] * w_temp_calc[1]) + s_psi[0]
                # DO 73 N=2,NMAX-1
                # R(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + S(N)
                for n_idx in range(1, NMAX - 1):
                    r_psi[n_idx] = -6 * iu * dt * (hb[n_idx-1] * w_temp_calc[n_idx-1] + ha[n_idx] * w_temp_calc[n_idx] + hb[n_idx] * w_temp_calc[n_idx+1]) + s_psi[n_idx]
                # R(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + S(NMAX)
                r_psi[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w_temp_calc[NMAX-2] + ha[NMAX-1] * w_temp_calc[NMAX-1]) + s_psi[NMAX-1]

                # Step 5: Calculate new S(N)
                # W(N) = (11*U(N)-14*V(N)+26*P(N)-14*Q(N)+11*R(N))/20.0
                # Uses current U, and the NEW V, P, Q, R
                w_temp_calc[:] = (11 * u_psi - 14 * v_psi + 26 * p_psi - 14 * q_psi + 11 * r_psi) / 20.0
                
                # S(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + T(1)
                s_psi[0] = -6 * iu * dt * (ha[0] * w_temp_calc[0] + hb[0] * w_temp_calc[1]) + t_psi[0]
                # DO 74 N=2,NMAX-1
                # S(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + T(N)
                for n_idx in range(1, NMAX - 1):
                    s_psi[n_idx] = -6 * iu * dt * (hb[n_idx-1] * w_temp_calc[n_idx-1] + ha[n_idx] * w_temp_calc[n_idx] + hb[n_idx] * w_temp_calc[n_idx+1]) + t_psi[n_idx]
                # S(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + T(NMAX)
                s_psi[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w_temp_calc[NMAX-2] + ha[NMAX-1] * w_temp_calc[NMAX-1]) + t_psi[NMAX-1]

                # Step 6: Calculate new T(N)
                # W(N) = (11*V(N)-14*P(N)+26*Q(N)-14*R(N)+11*S(N))/20.0
                # Uses the NEW V, P, Q, R, S
                w_temp_calc[:] = (11 * v_psi - 14 * p_psi + 26 * q_psi - 14 * r_psi + 11 * s_psi) / 20.0
                
                # T(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + U(1)
                t_psi[0] = -6 * iu * dt * (ha[0] * w_temp_calc[0] + hb[0] * w_temp_calc[1]) + u_psi[0]
                # DO 75 N=2,NMAX-1
                # T(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + U(N)
                for n_idx in range(1, NMAX - 1):
                    t_psi[n_idx] = -6 * iu * dt * (hb[n_idx-1] * w_temp_calc[n_idx-1] + ha[n_idx] * w_temp_calc[n_idx] + hb[n_idx] * w_temp_calc[n_idx+1]) + u_psi[n_idx]
                # T(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + U(NMAX)
                t_psi[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w_temp_calc[NMAX-2] + ha[NMAX-1] * w_temp_calc[NMAX-1]) + u_psi[NMAX-1]

                # Step 7: Calculate new U(N)
                # W(N) = (11*P(N)-14*Q(N)+26*R(N)-14*S(N)+11*T(N))/20.0
                # Uses the NEW P, Q, R, S, T
                w_temp_calc[:] = (11 * p_psi - 14 * q_psi + 26 * r_psi - 14 * s_psi + 11 * t_psi) / 20.0
                
                # U(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + V(1)
                u_psi[0] = -6 * iu * dt * (ha[0] * w_temp_calc[0] + hb[0] * w_temp_calc[1]) + v_psi[0]
                # DO 76 N=2,NMAX-1
                # U(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + V(N)
                for n_idx in range(1, NMAX - 1):
                    u_psi[n_idx] = -6 * iu * dt * (hb[n_idx-1] * w_temp_calc[n_idx-1] + ha[n_idx] * w_temp_calc[n_idx] + hb[n_idx] * w_temp_calc[n_idx+1]) + v_psi[n_idx]
                # U(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + V(NMAX)
                u_psi[NMAX-1] = -6 * iu * dt * (hb[NMAX-2] * w_temp_calc[NMAX-2] + ha[NMAX-1] * w_temp_calc[NMAX-1]) + v_psi[NMAX-1]

            #
            # plot
            #
            plt.gca().clear()
            plt.title('SC6')
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
            plt.savefig("sc6.png")
                
    except Exception as e:
        print(f"An error occurred during file operations: {e}")

    print("Program finished. Results written to sc6.dat")

# Call the main program function
if __name__ == '__main__':
    nroom_program()
