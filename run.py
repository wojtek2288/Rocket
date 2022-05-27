import os
import subprocess
import shutil
import time
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeClassifierCV

f = open("result_scores.txt", "a+")
f2 = open("result_coefs.txt", "a+")
for k in [("BasicMotions", 6), ("Epilepsy", 3), ("Handwriting", 3), ("NATOPS", 24)]:
    print(f"Data for {k[0]} {k[1]}")
    f.write(f"Data for {k[0]} {k[1]}\n")
    f2.write(f"Data for {k[0]} {k[1]}\n")
    DATASET_NAME = k[0]
    DIMENSIONS = k[1]
    KERNELS = 100
    SEED = 124
    USE_MEAN = "FALSE"
    ITERATION_COUNT = 10

    print("Deleting previous results")
    print()

    shutil.rmtree(".results", ignore_errors=True)

    print("Starting rocket algorithm")
    print()
    start = time.time()
    p = subprocess.Popen(f"Rscript rocket.R {DATASET_NAME} {DIMENSIONS} {KERNELS} {SEED} {USE_MEAN}", stdout=subprocess.PIPE)
    for line in iter(p.stdout.readline, b''):
        print(line.decode("utf-8"))
    p.stdout.close()
    p.wait()
    end = time.time()
    estimated_time = end - start

    print(f"Finished rocket algorithm in {estimated_time}, analyzing data")

    os.chdir(".results")

    y_train = np.loadtxt('Y_TRAIN.txt', comments="#", delimiter=",", skiprows=1, dtype='str')
    y_test = np.loadtxt('Y_TEST.txt', comments="#", delimiter=",", skiprows=1, dtype='str')

    scores = []
    coefs = []
    classes = []

    for i in range(ITERATION_COUNT):
        X_train_transform = np.loadtxt(f'{i + 1}_Result_TRAIN.txt', comments="#", delimiter=",", skiprows=1)
        X_test_transform = np.loadtxt(f'{i + 1}_Result_TEST.txt', comments="#", delimiter=",", skiprows=1)
        classifier = RidgeClassifierCV(alphas=np.logspace(-3, 3, 10))
        var = classifier.fit(X_train_transform, y_train)
        score = classifier.score(X_test_transform, y_test)
        coefs.append(var.coef_)
        classes = var.classes_
        print(f'Score for {i+1}: {score}')
        f.write(f'Score for {i+1}: {score}\n')
        scores.append(score)

    print(f"Result mean: {np.mean(scores)}")
    f.write(f"Result mean: {np.mean(scores)}\n")
    print(f"Result standard deviation: {np.std(scores)}")
    f.write(f"Result standard deviation: {np.std(scores)}\n")
    f.write(f"Algorithm took {estimated_time} seconds \n")
    f.write(f"Average time per iteration: {estimated_time / ITERATION_COUNT} seconds \n")

    shape = coefs[0].shape
    coefs2 = np.empty(shape)
    for i in range(coefs[0].shape[0]):
        for j in range(coefs[0].shape[1]):
            sum = 0
            for k in range(ITERATION_COUNT):
                sum += coefs[k][i,j]
            coefs2[i, j] = sum / ITERATION_COUNT

    df = pd.DataFrame(coefs2, index=classes)

    df.to_csv(f2, mode="a", sep="\t", header=False)

    os.chdir("..")
    # -- odchylenie standardowe
    # -- alpha
    # logowanie z ridge CV

f.close()
f2.close()
