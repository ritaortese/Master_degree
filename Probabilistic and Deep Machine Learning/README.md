# PML
# Fake Image Detection with Uncertainty Estimation

In this project, we explored two approaches for detecting fake images (fake vs. real) using Bayesian Convolutional Neural Networks (BCNNs), with a particular focus on estimating the model's uncertainty. The ability of a model to quantify its confidence in a prediction is crucial, especially in sensitive tasks like detecting deepfakes or manipulated content.

Both notebooks utilize an image dataset (presumably CIFAKE, given the "fake" vs. "real" classification context and relative data paths) and train BCNN models for binary classification. Images are resized to 64x64 pixels before processing.

## Project Files

1.  **`MonteCarloDropout.ipynb`**:
    * **Main Concept**: This notebook implements the Monte Carlo Dropout (MC-Dropout) technique to approximate Bayesian behavior and estimate uncertainty.
    * **Methodology**:
        * We used a standard CNN architecture (named `BayesianCNN` in the notebook) with convolutional layers, max pooling, and fully connected layers, including a Dropout layer (p=0.5).
        * During **training**, the model was trained conventionally using CrossEntropyLoss and the Adam optimizer. The notebook shows a decrease in loss over 5 epochs.
        * During **inference**, Dropout was kept active. We performed multiple forward passes ($T=30$) on the same input image, each with a different Dropout mask.
        * From these $T$ predictions, we calculated the **mean of probabilities** (for the final prediction), the **variance** (as a measure of epistemic uncertainty), and the **predictive entropy**.
    * **Key Metrics and Results**:
        * The model was evaluated on a test set, generating a confusion matrix and metrics like precision, recall, F1-score, and accuracy.
        * The reported accuracy is approximately **0.9283**.

2.  **`Variational_Inference.ipynb`**:
    * **Main Concept**: This notebook implements a BCNN using Variational Inference (VI) to learn distributions over the model weights, allowing for a more formal estimation of uncertainty.
    * **Methodology**:
        * The `BayesianCNN` architecture uses specific Bayesian layers from the `bayesian_torch` library (`Conv2dReparameterization`, `LinearReparameterization`). These layers learn a distribution (mean $\mu$ and standard deviation $\sigma$) for each weight, using the reparameterization trick.
        * **Training** minimizes the Evidence Lower Bound (ELBO), which combines the CrossEntropy loss with a KL divergence term between the learned variational distribution $q(\theta)$ and the prior $p(\theta)$ over the weights. The notebook shows a decrease in the total loss (ELBO) over 5 epochs.
        * During **inference**, similar to MC-Dropout but by sampling from the learned weight distributions, we performed multiple forward passes ($T=30$) to obtain a distribution of predictions.
        * We calculated the same uncertainty metrics: mean of probabilities, variance (epistemic uncertainty), and predictive entropy.
    * **Key Metrics and Results**:
        * The model was evaluated, generating a confusion matrix and classification metrics.
        * The reported accuracy is approximately **0.9057**.
        * An example visual of a prediction with its original and transformed image is included.
    * **Comparison with MC-Dropout (as noted in this notebook)**:
        * MC-Dropout achieved a slightly higher accuracy and more balanced F1-scores in its respective notebook.
        * The BCNN with VI appears to provide more nuanced and sharper uncertainty values, potentially being better at identifying borderline cases.
        * MC-Dropout's training was noted as slightly faster and more stable; VI is more computationally intensive but offers a more principled Bayesian treatment of uncertainty.

## Dataset

The project uses the CIFAKE dataset, which contains synthetically generated images (FAKE) and real-world images (REAL) for image classification tasks.

* You can find and download the dataset from Kaggle: [CIFAKE: Real and AI-Generated Synthetic Images](https://www.kaggle.com/datasets/birdy654/cifake-real-and-ai-generated-synthetic-images)

## General Conclusion

Both notebooks demonstrate how Bayesian techniques can be incorporated into CNNs for fake image detection, providing not only classifications but also valuable uncertainty metrics. We found that MC-Dropout is simpler to implement and offers good predictive performance. Variational Inference, while more complex, provides a more theoretically grounded and potentially more informative uncertainty estimation. The choice between them would depend on the specific application's requirements regarding accuracy, the interpretability of uncertainty, and computational resources.

Copyright © 2025 | Arahi Fernandez Monagas & Rita Ortese.
Project for the Probabilistic Machine Learning course (2024/2025), University of Trieste.

# Deep Learning 

# Project: Physics-Informed Deep Learning for Gene Expression Modeling

This project explores the use of Physics-Informed Neural Networks (PINNs) to model gene expression dynamics by coupling neural networks with biologically inspired differential equations describing mRNA and protein interactions. The PINN is trained using sparse and noisy time-series observations while enforcing the underlying physical constraints via differential equation residuals and initial conditions.

The framework addresses both forward modeling, where system parameters are known, and inverse modeling, where kinetic parameters are learned directly from data. The approach demonstrates how physics-based regularization improves generalization and parameter identifiability compared to purely data-driven models.

## Software & Tools: Python, PyTorch, NumPy, Matplotlib
## Models & Techniques Used: Physics-Informed Neural Networks; deep learning for dynamical systems; forward and inverse PINNs; automatic differentiation; ODE residual minimization; parameter estimation; sparse and noisy data learning.
