#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def load_files(quantity):
    path_to_data = os.path.join("data", quantity)
    file_names = os.listdir(path_to_data)
    nr_of_iterations = len(file_names)
    print(nr_of_iterations)

    out = []
    for file in sorted(file_names):
        path_to_file = os.path.join(path_to_data, file)
        with open(path_to_file) as fp:
            content = fp.readlines()
        values = [float(i[:-1]) for i in content]
        out.append(values)
    N = 10+2
    out = np.array(out)
    # print(out.shape)
    out = np.reshape(out, (nr_of_iterations, N, N))
    return out

def zero_pad(num, length=4):
    while len(str(num)) < length:
        num = f"0{num}"
    return num

if __name__ == "__main__":
    d = load_files("d")
    u = load_files("u")
    v = load_files("v")

    N = 10+2

    for iter_idx in tqdm(range(0, d.shape[0], 20)):
        plt.imshow(d[iter_idx])
        iter_idx = zero_pad(iter_idx)
        plt.title(iter_idx)
        plt.savefig(f"plots/d/{iter_idx}.png")

    # print(d)
