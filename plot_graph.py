from scipy import optimize
import scipy
import numpy
import matplotlib.pyplot as plt

def func(x, a, b):
        return a + b * numpy.log2(x)

def create_graph():
    with open('./results_n.txt') as f:
        y_pts_n = f.read().splitlines()
    with open('./results_n_2.txt') as f:
        y_pts_n_2 = f.read().splitlines()
    with open('./results_log_p.txt') as f:
        y_pts_log_p = f.read().splitlines()

    line_n, _ = scipy.optimize.curve_fit(lambda t,a,b: a+b*numpy.log2(t), range(2, 25, 2),  y_pts_n)
    line_n_2, _ = scipy.optimize.curve_fit(lambda t,a,b: a+b*numpy.log2(t), range(2, 49, 2),  y_pts_n_2)
    line_log_p, _ = scipy.optimize.curve_fit(lambda t,a,b: a+b*numpy.log2(t), range(2, 111, 2),  y_pts_log_p)

    fig, ax = plt.subplots()
    # ax.plot(range(2, 25, 2), y_pts_n, linestyle='None', marker='o', color='b', label=r"$p(n) = n$")
    ax.plot(range(2, 111, 2), func(range(2, 111, 2), *line_n),
            linewidth=2.0, linestyle='-', color='b', label=r"$Fitted Curve: n$"
    )
    # ax.plot(range(2, 49, 2), y_pts_n_2,
    #         linestyle='None', marker='o', color='r',  label=r"$p(n) = \frac{n}{2}$"
    # )
    ax.plot(range(2, 111, 2), func(range(2, 111, 2), *line_n_2),
            linewidth=2.0, linestyle='-', color='r', label=r"$Fitted Curve: \frac{n}{2}$"
    )
    # ax.plot(range(2, 111, 2), y_pts_log_p,
    #         linestyle='None', marker='o', color='g', label=r"$p(n) = \frac{n}{p} \geq \log \, p$"
    # )
    ax.plot(range(2, 111, 2), func(range(2, 111, 2), *line_log_p),
        linewidth=2.0, linestyle='-', color='g', label=r"$Fitted Curve: \frac{n}{p} \geq \log \, p$"
    )
    ax.set(xlabel='n - points count', ylabel=r'time $(\mu s)$',
           title='Line-of-Sight')
    ax.grid()
    ax.set_ylim(ymin=0)
    ax.legend(loc="upper left")
    fig.savefig("common_graph.pdf", format="pdf")
    plt.show()

if __name__ == '__main__':
    create_graph()