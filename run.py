import os
import subprocess
import shutil
import numpy as np
from sklearn.linear_model import RidgeClassifierCV

f = open("result.txt", "a")
for k in [("BasicMotions", 6), ("Epilepsy", 3), ("Handwriting", 3), ("NATOPS", 24)]:
    print(f"Data for {k[0]} k{1}")
    f.write(f"Data for {k[0]} k{1}\n")
    DATASET_NAME = k[0]
    DIMENSIONS = k[1]
    KERNELS = 100
    SEED = 124

    print("Deleting previous results")
    print()

    shutil.rmtree("results", ignore_errors=True)

    print("Starting rocket algorithm")
    print()

    p = subprocess.Popen(f"Rscript rocket.R {DATASET_NAME} {DIMENSIONS} {KERNELS} {SEED}", stdout=subprocess.PIPE)
    for line in iter(p.stdout.readline, b''):
        print(line.decode("utf-8"))
    p.stdout.close()
    p.wait()

    print("Finished rocket algorithm, analyzing data")

    os.chdir("results")

    y_train = np.loadtxt('Y_TRAIN.txt', comments="#", delimiter=",", skiprows=1, dtype='str')
    y_test = np.loadtxt('Y_TEST.txt', comments="#", delimiter=",", skiprows=1, dtype='str')

    scores = []

    for i in range(10):
        X_train_transform = np.loadtxt(f'{i + 1}_Result_TRAIN.txt', comments="#", delimiter=",", skiprows=1)
        X_test_transform = np.loadtxt(f'{i + 1}_Result_TEST.txt', comments="#", delimiter=",", skiprows=1)
        classifier = RidgeClassifierCV(alphas=np.logspace(-3, 3, 10))
        classifier.fit(X_train_transform, y_train)
        score = classifier.score(X_test_transform, y_test)
        print(f'Score for {i+1}: {score}')
        f.write(f'Score for {i+1}: {score}\n')
        scores.append(score)

    print(f"Result mean: {sum(scores) / len(scores)}")
    f.write(f"Result mean: {sum(scores) / len(scores)}\n")

    os.chdir("..")

f.close()
