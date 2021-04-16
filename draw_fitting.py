import numpy as np
import matplotlib.pyplot as plt

def read_measurement(file_name):
    '''
    Read txt file to generate a list of strokes

    Args:
        fileName(str): The name of the txt file, including .txt

    Returns:
        list_of_strokes: The list of strokes of the input character
    '''

    with open(file_name) as f:
        lines = f.read().splitlines()

    coordinates = []
    # coordinates.append(np.array([6.5, 0, 0]))
    for line in lines:
        if(line[0] == "#"):
            continue
        line = line.split(' ')
        line = line[0].split('\t')
        current_coordinate = []
        for num in line:
            if(num):
                current_coordinate.append(float(num))
        coordinates.append(np.array(current_coordinate))

    return np.array(coordinates)


def fit_and_plot(z, data, label, deg=1):
    coeff_drag = np.polyfit(z, data, deg)
    var_sim = np.polyval(coeff_drag, z)
    print(coeff_drag)
    plt.plot(z, data)
    plt.plot(z, var_sim)
    plt.xlabel("z /m")
    plt.ylabel(label + " /m")
    print("The average fitting error" + "of "+label+": ",
          np.linalg.norm(data-var_sim)/np.mean(data)/len(data))
    plt.show()


def fit(z, data, deg=1):
    coeff_drag = np.polyfit(z, data, deg)
    var_sim = np.polyval(coeff_drag, z)
    print(coeff_drag)
    print("The average fitting error" + ": ",
          np.linalg.norm(data-var_sim)/np.mean(data)/len(data))
    return var_sim


if __name__ == "__main__":
    tipLen = 0.03

    coordinates = read_measurement("line_measurement_2020_02_23.txt")
    z = (coordinates[:, 0])
    width = coordinates[:, 1]/1000
    drag = coordinates[:, 2]/1000
    offset = coordinates[:, 3]/1000

    # fit_and_plot(z,drag,"drag")
    # fit_and_plot(z,width,"width")
    # fit_and_plot(z,offset,"offset",2)

    new_function = drag+np.sqrt(np.power(tipLen-z, 2)+np.power(offset, 2))
    sim_offset2 = np.sqrt(np.power(tipLen-drag, 2)-np.power(tipLen-z, 2))
    # fit_and_plot(z,new_function,"tipLen")

    # fit_and_plot(z, drag+offset-z, "drag plus offset plus z")
    sim_width = fit(z, width)
    sim_drag = fit(z, drag)
    sim_offset = fit(z, offset, 2)
    sim_tipLen = fit(z, new_function)

    fig, (ax0,ax1,ax2) = plt.subplots(1, 3)
    ax0.plot(z*100, sim_width*100,c='#ff7f0e')
    ax0.scatter(z*100, width*100, s=20, c='#1f77b4')
    # axs[0, 0].set_title('Width')
    ax0.set(title='Width /cm')
    ax0.set(xlabel='z /cm')
    ax1.plot(z*100, sim_drag*100,c='#ff7f0e')
    ax1.scatter(z*100, drag*100, s=20, c='#1f77b4')
    ax1.set(title='Drag /cm')
    ax1.set(xlabel='z /cm')
    ax2.plot(z*100, sim_offset*100,c='#ff7f0e')
    ax2.scatter(z*100, offset*100, s=20, c='#1f77b4')
    ax2.set(title='Offset /cm')
    ax2.set(xlabel='z /cm')

    # for ax in axs.flat:
    #     ax.set(xlabel='z /m')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #     ax.label_outer()

    plt.show()
