import numpy as np
import matplotlib.pyplot as plt

# Define the degree of the polynomial
for NDEG in [1,2,4,6]:
    try:
        with open('stability'+str(NDEG)+'.dat', 'w') as f_out:
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
