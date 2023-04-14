Updated on 12/04/2023. 

This repository hosts the code for the deep inverse Rosenblatt transport methods.

Load both `src` and `external` folders to use the code. The `MultiIndices` class in the `external` folder is used by `SparseFun` for sparse polynomial approximations. 

See the repository `https://github.com/DeepTransport/deep-tensor-examples` for more examples and `https://github.com/dolgov/TT-IRT` for a different implementation.

Important change of interfaces:
* All polynomial classes have name started with Capital letters. 
* All diagonal map classes also have name started with Capital letters. This affects `EmpiricalMap`, `UniformMap`, and the abstract `DiagonalMap`.
* The DIRT (including TTDIRT and SparseDIRT) class now has a simplified interface. 
    * In the simplest case, one only needs to passin the density function, the parameter dimension and the approximation domain. This automatically override the default FTT option (in DIRT class), use a default 2nd order Lagrange polynimal basis, and use the Gaussian reference. See `example_dirt_default.m`.
    * Alternatively, one can specify the density function, the parameter dimension and the approximation polynomial basis (with a given domain). This automatically override the default FTT option (in DIRT class) and use the Gaussian reference. See `example_conditional_dirt.m`.
    * One can also have the full specification, with polynomial basis for each layer, reference measure, FTT options, etc.  See `example_dirt.m`.
* All the help files are due to updates.

References: 
* For SIRT, DIRT, and IRT: 
    * Cui, Dolgov and Zahm (2023). Self-reinforced polynomial approximation methods for concentrated probability densities. arXiv preprint: 2303.02554
    * Cui, Dolgov and Scheichl (2022). Deep importance sampling using tensor-trains with application to a priori and a posteriori rare event estimation. arXiv preprint: 2209.01941
    * Cui, Dolgov and Zahm (2023). Scalable conditional deep inverse Rosenblatt transports using tensor-trains and gradient-based dimension reduction. Journal of Computational Physics, 112103
    * Cui and Dolgov (2022). Deep composition of tensor trains using squared inverse Rosenblatt transports. Foundation of Computational Mathematics, Foundations of Computational Mathematics 22 (6), 1863-1922
    * Dolgov, Anaya-Izquierdo, Fox and Scheichl (2020). Approximation and sampling of multivariate probability distributions in the tensor train decomposition. Statistics and Computing 30(3), 603-625.
* For building the tensor train using AMEN:
    * Dolgov, Savostyanov (2014). Alternating minimal energy methods for linear systems in higher dimensions. SIAM Journal on Scientific Computing 36(5), A2248-A2271.
* For functional tensor train:
    * Gorodetsky, Karaman and Marzouk (2018). A continuous analogue of the tensor-train decomposition. Computer Methods in Applied Mechanics and Engineering 347, 59-84.
    * Bigoni, Engsig-Karup, Marzouk (2016). Spectral tensor-train decomposition. SIAM Journal on Scientific Computing 38(4), A2405-A2439.

The `MultiIndices` class is from the `ApproximationToolbox` package of Nouy:
    * https://github.com/anthony-nouy/ApproximationToolbox
    * https://anthony-nouy.github.io/ApproximationToolbox/
