import numpy as np
import matplotlib.pyplot as plt

# Define the degree of the polynomial

for NDEG in [1,2,4,6]:
    # Open the output file
    # Using 'with open' ensures the file is properly closed even if errors occur.
    try:
        with open('stability'+str(NDEG)+'.dat', 'w') as f_out:
            # Loop for parameter A from 0.0 to 1.0 with a step of 0.01
            # np.arange handles floating-point steps more reliably than a simple while loop.
            # Adding a small epsilon to the stop value to ensure 1.0 is included due to
            # floating point precision issues.
            plot_x = np.arange(0.0, 2.0001, 0.01)
            plot_y = [[] for _ in range(NDEG) ]
            for a_val in plot_x:
                # Define the imaginary unit
                iu = 0.0 + 1.0j

                # Define the coefficients of the polynomial z^2 + (2*i*A)*z - 1 = 0
                # numpy.roots expects coefficients in descending order of power: [c_n, c_{n-1}, ..., c_0]
                # For a quadratic equation c2*z^2 + c1*z + c0 = 0
                # Here, c2 = 1, c1 = 2*iu*a_val, c0 = -1
                coefficients    = {}
                coefficients[1] = [1.0, iu * a_val -1.0]
                coefficients[2] = [1.0, 2.0 * iu * a_val, -1.0]
                coefficients[4] = [1.0, +8./3.*iu*a_val, -4./3.*iu*a_val, +8./3.*iu*a_val, -1.0]
                coefficients[6] = [1.0, +33./10.*iu*a_val, -21./5.*iu*a_val, +39./5.*iu*a_val, -21./5.*iu*a_val, +33./10.*iu*a_val, -1.0]

                # Calculate the roots of the polynomial
                # np.roots returns an array of roots
                roots = np.roots(coefficients[NDEG])

                # Iterate through the calculated roots
                for i in range(NDEG):
                    root = roots[i]
                    # Calculate the stability check: |root| - 1
                    stability_check = abs(root) - 1.0
                    plot_y[i].append(stability_check)
                    # Print to standard output (console)
                    # Mimicking Fortran's E15.7 format for scientific notation
                    print(f"{a_val:15.7E} {root.real:15.7E} {root.imag:15.7E} {stability_check:15.7E}")

                # Write to the output file
                # The Fortran code writes all roots on one line for each A value.
                # We will format the output similarly.
                # The Fortran format '30E15.7' suggests 30 fields, each 15 chars wide, 7 decimal places.
                # For NDEG=2, we have A, (REAL(R1), IMAG(R1), ABS(R1)), (REAL(R2), IMAG(R2), ABS(R2))
                # This makes 1 + 3*NDEG = 1 + 3*2 = 7 fields. The Fortran format is quite generous.
                # We'll stick to the required number of fields for our output.
                file_output_line = f"{a_val:15.7E}"
                for i in range(NDEG):
                    root = roots[i]
                    file_output_line += f" {root.real:15.7E} {root.imag:15.7E} {abs(root):15.7E}"
                f_out.write(file_output_line + "\n")

            print("A")
            #
            # plot
            #
            plt.gca().clear()
            plt.title("stability"+str(NDEG))
            plt.xlabel(r'α')
            plt.ylabel(r'|λ|-1')
            for i in range(NDEG):
                #print(f'i={i}')
                #print(plot_x)
                #print(plot_y)
                plt.plot(plot_x,plot_y[i], marker='o', linestyle='None')
            plt.show()
            plt.savefig("stability"+str(NDEG)+".png")
            
    except Exception as e:
        print(f"An error occurred: {e}")

print("Program finished. Results written to stability*.dat")
