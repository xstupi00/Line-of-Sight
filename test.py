import numpy as np
import pandas as pd
import random
import subprocess
from termcolor import colored
import matplotlib.pyplot as plt

MIN_NUMBERS = 8
MAX_NUMBERS = 24
NUMBER_STEP = 1
TEST_REPEAT = 5
MIN_RANGE = 0
MAX_RANGE = 100
FLOAT_MIN = -3.40282e+38

def compute_angles(numbers):
    angles = [round(np.arctan( (number - numbers[0]) / (idx+1)), 5)  for idx, number in enumerate(numbers[1:])]
    angles.insert(0, FLOAT_MIN)
    angles_series = pd.Series(angles)
    # max-scan
    min_previous_max = angles_series.cummax().tolist()
    # remove overall maximum (transform to max-prescan)
    del min_previous_max[-1]
    # add neutral item I (transform to max-prescan)
    min_previous_max.insert(0, FLOAT_MIN)
    return min_previous_max

def run_test_check():
    for numbers_count in range (MIN_NUMBERS, MAX_NUMBERS, NUMBER_STEP):
        for _ in range (0, TEST_REPEAT):
            input_num = [random.randint(MIN_RANGE, MAX_RANGE) for _ in range (0, numbers_count)]
            input_str = ','.join([str(number) for number in input_num])
            reference_output = compute_angles(input_num)
            for option in range(1, 4):
                out = subprocess.check_output(["./test.sh", input_str, str(option)])
                out_num = out.decode("utf-8").split("\n")[0].split(",")
                out_num = \
                    [round(float(number), 5) if float(number) != FLOAT_MIN else float(number) for number in out_num]
                if out_num == reference_output:
                    print(colored("Test (" + str(numbers_count) + " - " + str(option) + ") successful.", 'green'))
                else:
                    print(colored("Test (" + str(numbers_count) + " - " + str(option) + ") unsuccessful.", 'red'))
                    print("-----------------------------------------------------")
                    print(input_num)
                    print(reference_output)
                    print(out_num)
                    print("-----------------------------------------------------")

def create_graph(elapsed_time, samples):
    fig, ax = plt.subplots()
    print(samples, elapsed_time[0])
    ax.plot(samples, elapsed_time[0], linestyle='-', marker='o', color='b')
    ax.plot(samples, elapsed_time[1], linestyle='-', marker='o', color='r')
    ax.plot(samples, elapsed_time[2], linestyle='-', marker='o', color='g')
    ax.set(xlabel='n - points count', ylabel='time (us)',
           title='Line-of-Sight')
    ax.grid()
    fig.savefig("plot.png")
    plt.show()

def rewrite_results(results, filename):
    with open(filename, 'r+') as f:
        _ = f.read()
        f.seek(0, 0)
        for result in results:
            f.write(str(result) + '\n')

def run_test_measure():
    elapsed_times = [[],[],[]]
    for numbers_count in range (MIN_NUMBERS, MAX_NUMBERS, NUMBER_STEP):
        print(numbers_count)
        sub_times_log = []
        sub_times_n_2 = []
        sub_times_n = []
        for _ in range (0, TEST_REPEAT):
            input_num = [random.randint(MIN_RANGE, MAX_RANGE) for _ in range (0, numbers_count)]
            input_str = ','.join([str(number) for number in input_num])
            for option in range(1, 4):
                out = subprocess.check_output(["./test.sh", input_str, str(option)])
                if option == 1:
                    sub_times_log.append(float(out.decode("utf-8").split("\n")[0]))
                elif option == 2:
                    sub_times_n_2.append(float(out.decode("utf-8").split("\n")[0]))
                elif option == 3:
                    sub_times_n.append(float(out.decode("utf-8").split("\n")[0]))
        elapsed_times[0].append(min(sub_times_log))
        elapsed_times[1].append(min(sub_times_n_2))
        elapsed_times[2].append(min(sub_times_n))
        rewrite_results(elapsed_times[0], "results_1.txt")
        rewrite_results(elapsed_times[1], "results_2.txt")
        rewrite_results(elapsed_times[2], "results_3.txt")
    create_graph(elapsed_times, range(MIN_NUMBERS, MAX_NUMBERS, NUMBER_STEP))

if __name__ == '__main__':
    run_test_check()
    # run_test_measure()