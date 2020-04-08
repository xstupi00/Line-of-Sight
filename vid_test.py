from subprocess import Popen, PIPE
from argparse import ArgumentParser
import math
import random
import itertools
from termcolor import colored

FAIL = '\033[91m\033[1m'
GREEN = '\033[32m\033[1m'
ENDC = '\033[0m'
BOLD = '\033[1m'


def main():
    parser = ArgumentParser()
    parser.add_argument("-m", "--mpi", type=str, default="mpirun")
    parser.add_argument("-e", "--executable", type=str, default="vid")
    parser.add_argument("-l", "--limit", type=int, default=15)

    argv = parser.parse_args()
    if not argv.executable:
        return

    altitudes = [random.randint(1, 1024)]

    args = [argv.mpi, '--hostfile', 'hostfile', '-np', '', argv.executable, '']
    for inputSize in range(1, 31):
        print("[Input size]: ", inputSize)
        # permutations = itertools.permutations(altitudes)

        for i in range(0, min(argv.limit, math.factorial(inputSize))):
            if i == 0:
                altitudes.sort()
            elif i == 1:
                altitudes.sort(reverse=True)
            else:
                random.shuffle(altitudes)

            args[6] = ','.join(map(str, altitudes))
            print("[Input string]: ", {args[6]})
            angles = [-math.inf] + [math.atan((x - altitudes[0]) / float(i))
                    for i, x in enumerate(altitudes[1:], 1)]

            max_prescan = [-math.inf] + list(itertools.accumulate(angles, max))
            max_prescan.pop()
            visibilities = ["v" if (angle > max_prev_angle) else "u"
                            for (angle, max_prev_angle) in zip(angles, max_prescan)]
            visibilities[0] = "_"

            ref_output = ",".join(visibilities)
            for procCount in range(1, inputSize + 1):
                args[4] = str(procCount)
                process = Popen(args, stdout=PIPE, stderr=PIPE)
                output = process.communicate()[0].decode('utf-8').rstrip('\n')
                if output != ref_output:
                    print(colored(ref_output, "red"))
                    print(colored(output, "red"))
                    print(procCount)
                    print("------------------------------------")

        altitudes.append(random.randint(1, 1024))


if __name__ == "__main__":
    main()
